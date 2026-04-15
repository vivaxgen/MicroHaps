"""
gff3_csq.py
-----------
GFF3 parser and genomic <-> amino acid position converter.
Designed for Plasmodium drug-resistance gene analysis (bcftools csq context).
"""

from __future__ import annotations

import dataclasses
import warnings
from collections import defaultdict
from typing import Callable, Optional


# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

# A callable that returns the reference base (uppercase, single char) at a
# 1-based genomic position on a given sequence/contig.
# Compatible with pyfaidx:  lambda seqid, pos: str(fasta[seqid][pos-1:pos])
# Compatible with pysam:    lambda seqid, pos: fasta.fetch(seqid, pos-1, pos).upper()
# Compatible with a plain dict (for tests): lookup_dict.get
RefLookup = Callable[[str, int], str]


from typing import NamedTuple

class GenomicPosition(NamedTuple):
    """A fully-qualified genomic coordinate: chromosome name + 1-based position."""
    seqid: str
    pos:   int

    def __str__(self) -> str:
        return f"{self.seqid}:{self.pos}"


@dataclasses.dataclass
class CodonResult:
    """
    Result of translating one codon affected by one or more variants.

    Attributes
    ----------
    aa_pos      : 1-based amino acid position in the protein
    ref_codon   : 3-base reference codon string (on the + strand of the CDS)
    alt_codon   : 3-base alternate codon string (with variant bases applied)
    ref_aa      : single-letter reference amino acid
    alt_aa      : single-letter alternate amino acid
    ref_filled  : genomic positions that were absent from the input and were
                  filled from the reference sequence (empty if all 3 were provided)
    is_synonymous : True when ref_aa == alt_aa
    is_stop_gained : True when alt_aa == '*' and ref_aa != '*'
    """
    aa_pos       : int
    ref_codon    : str
    alt_codon    : str
    ref_aa       : str
    alt_aa       : str
    ref_filled   : list[GenomicPosition]  # codon positions filled from reference

    @property
    def is_synonymous(self) -> bool:
        return self.ref_aa == self.alt_aa

    @property
    def is_stop_gained(self) -> bool:
        return self.alt_aa == '*' and self.ref_aa != '*'

    def csq_string(self) -> str:
        """bcftools csq-style amino acid change string, e.g. '86N>86Y'."""
        if self.is_synonymous:
            return f"{self.aa_pos}{self.ref_aa}"
        return f"{self.aa_pos}{self.ref_aa}>{self.aa_pos}{self.alt_aa}"

    def __str__(self) -> str:
        filled_note = (
            f" [ref-filled: {self.ref_filled}]" if self.ref_filled else ""
        )
        return (
            f"AA {self.aa_pos}: {self.ref_codon}>{self.alt_codon} "
            f"({self.ref_aa}>{self.alt_aa}){filled_note}"
        )


# ---------------------------------------------------------------------------
# Codon / amino acid tables
# ---------------------------------------------------------------------------

CODON_TABLE: dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

UNKNOWN_BASE = 'N'
UNKNOWN_AA   = 'X'
COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def codon2aa(codon: str) -> str:
    """Translate a 3-character codon to a single-letter amino acid."""
    codon = codon.upper()
    return CODON_TABLE.get(codon, UNKNOWN_AA)


def reverse_complement_codon(codon: str) -> str:
    """Return the reverse complement of a 3-character codon string."""
    return codon.translate(COMP)[::-1]

def complement_codon(codon: str) -> str:
    """Return the complement of a codon string (no reversal)."""
    return codon.translate(COMP)


# ---------------------------------------------------------------------------
# GFF3 data structures
# ---------------------------------------------------------------------------

class GFF3Feature:
    """Represents a single GFF3 feature (typically a CDS)."""

    def __init__(self, line: str):
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            raise ValueError(f"Invalid GFF3 line: {line!r}")
        self.seqid  = parts[0]
        self.source = parts[1]
        self.type   = parts[2]
        self.start  = int(parts[3])                           # 1-based inclusive
        self.end    = int(parts[4])                           # 1-based inclusive
        self.score  = float(parts[5]) if parts[5] != '.' else None
        self.strand = parts[6]
        self.phase  = int(parts[7]) if parts[7] != '.' else None  # 0, 1, or 2
        self.attributes: dict[str, str] = {}
        for item in parts[8].split(';'):
            if '=' in item:
                k, v = item.split('=', 1)
                self.attributes[k.strip()] = v.strip()

    def get_attribute(self, key: str) -> str:
        return self.attributes.get(key, '')

    @property
    def length(self) -> int:
        """Length of this feature in bases (inclusive)."""
        return self.end - self.start + 1

    def contains(self, pos: int) -> bool:
        """True if 1-based genomic *pos* falls within this feature."""
        return self.start <= pos <= self.end

    def __repr__(self) -> str:
        return (
            f"GFF3Feature({self.type} {self.seqid}:{self.start}-{self.end} "
            f"{self.strand} phase={self.phase} "
            f"ID={self.get_attribute('ID')!r})"
        )


class GFF3:
    """Minimal GFF3 reader – loads features for codon lookups."""

    def __init__(self, path: str = ''):
        self.features: list[GFF3Feature] = []
        if path:
            self.read_file(path)

    def read_file(self, path: str) -> None:
        from smart_open import open as smart_open
        with smart_open(path) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                self.features.append(GFF3Feature(line))

    def query_features(
        self,
        query_func: Callable[[GFF3Feature], bool],
        ftype: Optional[str] = 'CDS',
    ) -> list[GFF3Feature]:
        """
        Return features matching *query_func*.
        If *ftype* is given, only features of that type are considered.
        Pass ftype=None to search across all feature types.
        """
        return [
            f for f in self.features
            if (ftype is None or f.type == ftype) and query_func(f)
        ]

    def get_transcripts(self, gene_id: str) -> dict[str, list[GFF3Feature]]:
        """
        Return CDS features grouped by transcript ID for *gene_id*.

        Matches on the Parent= attribute (normal GFF3) and falls back to
        ID= so that single-feature flat GFFs are also handled.

        Returns
        -------
        dict mapping transcript_id -> list[GFF3Feature]
            Exons within each transcript are sorted in transcription order:
            ascending start for + strand, descending start for - strand.

        Raises
        ------
        ValueError
            If inconsistent strands are found across a transcript's CDS
            features, which would indicate a malformed GFF3.
        """
        raw: dict[str, list[GFF3Feature]] = defaultdict(list)

        cds_features = self.query_features(
            lambda f: (
                gene_id in f.get_attribute('Parent') or
                gene_id in f.get_attribute('ID')
            ),
            ftype='CDS',
        )

        for f in cds_features:
            parent = f.get_attribute('Parent') or f.get_attribute('ID') or gene_id
            raw[parent].append(f)

        transcripts: dict[str, list[GFF3Feature]] = {}
        for tid, exons in raw.items():
            # Sanity-check: all CDS in a transcript must share one strand
            strands = {f.strand for f in exons}
            if len(strands) > 1:
                raise ValueError(
                    f"Transcript {tid!r} has mixed strands {strands} — "
                    "check your GFF3 file."
                )
            strand = next(iter(strands))
            exons.sort(key=lambda f: f.start, reverse=(strand == '-'))
            transcripts[tid] = exons

        return transcripts


# ---------------------------------------------------------------------------
# Position conversion
# ---------------------------------------------------------------------------

class TranscriptCoords:
    """
    Genomic <-> CDS <-> amino acid coordinate converter for one transcript.

    Parameters
    ----------
    exons : list[GFF3Feature]
        CDS features in transcription order (as returned by
        GFF3.get_transcripts). Must all share the same strand.
    """

    def __init__(self, transcript_id: str, exons: list[GFF3Feature]):
        if not exons:
            raise ValueError("exons list is empty")
        self.transcript_id = transcript_id
        self.exons  = exons
        self.strand = exons[0].strand

        # Phase offset: bases to skip at the very start of the first CDS exon.
        # Usually 0 for Plasmodium but handled for correctness.
        self.phase_offset: int = exons[0].phase or 0

        # Total CDS length after phase trimming
        self.cds_length: int = (
            sum(f.length for f in exons) - self.phase_offset
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _cds_offset_of(self, genomic_pos: int) -> Optional[int]:
        """
        Convert a 1-based genomic position to a 1-based CDS offset.
        Returns None if the position is not within any CDS exon.

        CDS offset counts from the first translated base (after phase trim),
        so offset 1 = first base of start codon.
        """
        accumulated = 0
        for i, exon in enumerate(self.exons):
            if exon.contains(genomic_pos):
                if self.strand == '+':
                    raw_offset = accumulated + (genomic_pos - exon.start + 1)
                else:
                    raw_offset = accumulated + (exon.end - genomic_pos + 1)

                # Subtract phase from the first exon only
                cds_offset = raw_offset - self.phase_offset
                return cds_offset if cds_offset > 0 else None
            accumulated += exon.length

        return None  # position not in any CDS exon

    def _genomic_of_cds_offset(self, cds_offset: int) -> Optional[int]:
        """
        Convert a 1-based CDS offset to a genomic position.
        Returns None if the offset is out of range.
        """
        if cds_offset < 1 or cds_offset > self.cds_length:
            return None

        # Re-introduce the phase trim when walking exons
        raw_target = cds_offset + self.phase_offset
        accumulated = 0
        for exon in self.exons:
            if accumulated + exon.length >= raw_target:
                offset_in_exon = raw_target - accumulated - 1  # 0-based into exon
                if self.strand == '+':
                    return exon.start + offset_in_exon
                else:
                    return exon.end - offset_in_exon
            accumulated += exon.length

        return None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def genomic_to_aa(self, genomic_pos: int) -> Optional[tuple[int, int]]:
        """
        Convert a 1-based genomic position to (aa_position, codon_base).

        codon_base is 1, 2, or 3 — which base within the codon this is.
        Returns None if the position is intronic or outside the CDS.

        Example
        -------
        >>> tc.genomic_to_aa(423456)
        (101, 2)   # AA position 101, second base of its codon
        """
        cds_offset = self._cds_offset_of(genomic_pos)
        if cds_offset is None:
            return None
        aa_pos     = (cds_offset - 1) // 3 + 1
        codon_base = (cds_offset - 1) %  3 + 1
        return aa_pos, codon_base

    def aa_to_genomic(self, aa_pos: int) -> Optional[list[GenomicPosition]]:
        """
        Convert a 1-based amino acid position to the genomic coordinates of
        its three codon bases.

        Returns a list of 3 GenomicPosition(seqid, pos) named tuples in
        transcription order, or None if the AA position is out of range.

        Note: the three positions may span two exons (split codon) and on
        the - strand they will be in descending coordinate order.

        Example
        -------
        >>> tc.aa_to_genomic(86)
        [GenomicPosition(seqid='Pf3D7_05_v3', pos=423201),
         GenomicPosition(seqid='Pf3D7_05_v3', pos=423202),
         GenomicPosition(seqid='Pf3D7_05_v3', pos=423203)]
        """
        max_aa = self.cds_length // 3
        if aa_pos < 1 or aa_pos > max_aa:
            return None

        seqid = self.exons[0].seqid
        codon_cds_offsets = [(aa_pos - 1) * 3 + i for i in range(1, 4)]
        positions = [self._genomic_of_cds_offset(o) for o in codon_cds_offsets]

        if any(p is None for p in positions):
            return None  # shouldn't happen if aa_pos is in range
        return [GenomicPosition(seqid, p) for p in positions]  # type: ignore[arg-type]

    def codon_genomic_range(self, aa_pos: int) -> Optional[tuple[str, int, int]]:
        """
        Return (seqid, min_pos, max_pos) spanning the codon for *aa_pos*.
        Useful for interval lookups regardless of strand.
        """
        positions = self.aa_to_genomic(aa_pos)
        if positions is None:
            return None
        seqid = positions[0].seqid
        coords = [p.pos for p in positions]
        return seqid, min(coords), max(coords)

    def translate_variant(
        self,
        aa_pos: int,
        ref_lookup: RefLookup,
        variant_bases: Optional[dict[int, str]] = None,
    ) -> CodonResult:
        """
        Translate a specific amino acid position, substituting any provided
        variant bases and filling the rest from the reference sequence.

        Parameters
        ----------
        aa_pos : int
            1-based amino acid position to translate.
        ref_lookup : RefLookup
            Callable(seqid, genomic_pos) -> uppercase reference base.
            seqid is taken from the first exon's GFF3Feature.seqid.
        variant_bases : dict[int, str], optional
            Mapping of 1-based genomic position -> alternate base (A/C/G/T).
            Only positions belonging to *aa_pos*'s codon are used; others
            are silently ignored.
            If None or empty, all three codon positions are filled from the
            reference — useful for looking up the reference amino acid.

        Returns
        -------
        CodonResult
            Single result for *aa_pos*. CodonResult.ref_filled lists any
            codon positions that were absent from *variant_bases* and were
            filled from the reference.

        Raises
        ------
        ValueError
            If *aa_pos* is out of range for this transcript.

        Warnings
        --------
        UserWarning
            Issued if one or two codon positions were filled from the
            reference rather than supplied in *variant_bases*.

        Examples
        --------
        # Full codon supplied — no warning
        >>> tc.translate_variant(86, ref_lookup, {423201: 'A', 423202: 'C', 423203: 'T'})

        # Only one base supplied — warning issued, other two filled from ref
        >>> tc.translate_variant(86, ref_lookup, {423202: 'G'})

        # No variants — returns reference amino acid, no warning
        >>> tc.translate_variant(86, ref_lookup)
        """
        codon_gpos = self.aa_to_genomic(aa_pos)
        if codon_gpos is None:
            raise ValueError(
                f"AA position {aa_pos} is out of range for transcript "
                f"{self.transcript_id!r} (max {self.cds_length // 3})."
            )

        provided = variant_bases or {}

        ref_bases:   list[str]             = []
        alt_bases:   list[str]             = []
        filled_here: list[GenomicPosition] = []

        # codon_gpos is already in CDS/transcription order (3 GenomicPositions)
        for gp in codon_gpos:
            ref_base = ref_lookup(gp.seqid, gp.pos).upper()
            ref_bases.append(ref_base)

            if gp.pos in provided:
                alt_bases.append(provided[gp.pos].upper())
            else:
                alt_bases.append(ref_base)
                # Only flag as filled if the caller actually supplied some
                # variants (otherwise all-ref is the intended behaviour).
                if provided:
                    filled_here.append(gp)

        # On - strand: bases are collected in transcription order but are
        # minus-strand bases — reverse-complement the triplet to get the
        # mRNA codon before translation.
        ref_codon_raw = ''.join(ref_bases)
        alt_codon_raw = ''.join(alt_bases)

        if self.strand == '-':
            ref_codon = complement_codon(ref_codon_raw)
            alt_codon = complement_codon(alt_codon_raw)
            #ref_codon = reverse_complement_codon(ref_codon_raw)
            #alt_codon = reverse_complement_codon(alt_codon_raw)
        else:
            ref_codon = ref_codon_raw
            alt_codon = alt_codon_raw

        if filled_here:
            filled_str = ', '.join(str(gp) for gp in filled_here)
            warnings.warn(
                f"[{self.transcript_id}] AA {aa_pos}: "
                f"{len(filled_here)} codon position(s) filled from reference: "
                f"{filled_str}. "
                "Result may not reflect a compound variant.",
                UserWarning,
                stacklevel=2,
            )

        return CodonResult(
            aa_pos    = aa_pos,
            ref_codon = ref_codon,
            alt_codon = alt_codon,
            ref_aa    = codon2aa(ref_codon),
            alt_aa    = codon2aa(alt_codon),
            ref_filled = filled_here,
        )

    def __repr__(self) -> str:
        return (
            f"TranscriptCoords({self.transcript_id!r} "
            f"strand={self.strand!r} "
            f"exons={len(self.exons)} "
            f"cds_length={self.cds_length} "
            f"max_aa={self.cds_length // 3})"
        )


# ---------------------------------------------------------------------------
# Convenience: gene-level entry point
# ---------------------------------------------------------------------------

class GeneCoords:
    """
    Coordinate converter for all transcripts of a single gene.

    Wraps GFF3.get_transcripts + TranscriptCoords into one object.

    Parameters
    ----------
    gff     : GFF3 instance (already loaded)
    gene_id : e.g. "PF3D7_0523000"

    Usage
    -----
    >>> gc = GeneCoords(gff, "PF3D7_0523000")
    >>> gc.genomic_to_aa(423456)
    {'PF3D7_0523000.1': (101, 2)}
    >>> gc.aa_to_genomic(86)
    {'PF3D7_0523000.1': [423201, 423202, 423203]}
    """

    def __init__(self, gff: GFF3, gene_id: str):
        self.gene_id = gene_id
        transcripts  = gff.get_transcripts(gene_id)
        if not transcripts:
            raise ValueError(f"No CDS features found for gene {gene_id!r}")
        if len(transcripts) > 1:
            warnings.warn(
                f"{gene_id!r} has {len(transcripts)} transcripts: "
                f"{list(transcripts)}. Results returned per transcript.",
                stacklevel=2,
            )
        self.transcripts: dict[str, TranscriptCoords] = {
            tid: TranscriptCoords(tid, exons)
            for tid, exons in transcripts.items()
        }

    def genomic_to_aa(
        self, genomic_pos: int
    ) -> dict[str, Optional[tuple[int, int]]]:
        """
        Return {transcript_id: (aa_pos, codon_base)} for all transcripts.
        Value is None for transcripts where the position is intronic/absent.
        """
        return {
            tid: tc.genomic_to_aa(genomic_pos)
            for tid, tc in self.transcripts.items()
        }

    def aa_to_genomic(
        self, aa_pos: int
    ) -> dict[str, Optional[list[GenomicPosition]]]:
        """
        Return {transcript_id: [GenomicPosition, ...]} for all transcripts.
        Value is None where the AA position is out of range.
        """
        return {
            tid: tc.aa_to_genomic(aa_pos)
            for tid, tc in self.transcripts.items()
        }

    def translate_variant(
        self,
        aa_pos: int,
        ref_lookup: RefLookup,
        variant_bases: Optional[dict[int, str]] = None,
    ) -> dict[str, CodonResult]:
        """
        Translate a specific amino acid position across all transcripts.

        Parameters
        ----------
        aa_pos : int
            1-based amino acid position to translate.
        ref_lookup : RefLookup
            Callable(seqid, genomic_pos) -> uppercase reference base.
        variant_bases : dict[int, str], optional
            {genomic_pos: alt_base} — passed through to each transcript.
            If None, returns the reference amino acid for each transcript.

        Returns
        -------
        dict mapping transcript_id -> CodonResult
            Transcripts where *aa_pos* is out of range are omitted.

        Example
        -------
        >>> gc = GeneCoords(gff, "PF3D7_0523000")
        >>> results = gc.translate_variant(86, ref_lookup, {423202: 'G'})
        >>> for tid, r in results.items():
        ...     print(f"{tid}  {r.csq_string()}")
        PF3D7_0523000.1  86N>86Y
        """
        out: dict[str, CodonResult] = {}
        for tid, tc in self.transcripts.items():
            try:
                out[tid] = tc.translate_variant(aa_pos, ref_lookup, variant_bases)
            except ValueError:
                pass  # aa_pos out of range for this transcript — omit
        return out

    def __repr__(self) -> str:
        return (
            f"GeneCoords({self.gene_id!r} "
            f"transcripts={list(self.transcripts)})"
        )
