Implementing microhaplotypes calling for different panel
=========================================================

This document describes the implementation of microhaplotypes calling for different panel. 
The pipeline is designed to be flexible and can be adapted to different panels by modifying the configuration files and the reference files.
The following sections outline the steps involved in implementing microhaplotypes calling for different panel.

Configuration Files
----------------------

The configuration files are located in the ``path-to-vvg-Microhaps/envs/MicroHaps/configs`` directory.
To create a new panel, create a new directory in the ``configs/refs`` to store the reference files and configuration.
E.g., ``path-to-vvg-Microhaps/envs/MicroHaps/configs/refs/new_panel``.
Additionally, the ``configs/refs/`` directory should contain a yaml file referencing the relative path to the actual configuration file:
E.g., given that the name of the new panel is ``new_panel``, the file should be named ``new_panel.yaml`` and created in the ``configs/refs`` directory containing the following content:

.. code-block:: yaml

    refs/new_panel/<actual configuration>.yaml

The actual configuration file should be created in the ``configs/refs/new_panel`` directory.
The configuration file should contain the following sections:

.. code-block:: yaml

    #
    # new_panel.yaml - configuration file for the new microhaps panel
    #
    # --- NGS-Pipeline settings ---

    optical_dedup: False
    instrument: generic
    platform: generic
    libprep: generic

    min_read_qual: 15
    read_length: 100

    refseq_file: configs/refs/Pv/PvP01_v2/PvP01_v2.fasta  # To be modified, see, Reference Files section

    # perform variant calling on all chromosomes as single workflow
    complete_region: PvP01_all_v2  # To be modified

    # perform targeted variant calling
    targetregion_file: configs/refs/pv-vvg/Microhaps_inserts.bed # To be modified, see, Reference Files section

    regions: # To be modified
    - PvP01_01_v2
    - PvP01_02_v2
    - PvP01_03_v2
    - PvP01_04_v2
    - PvP01_05_v2
    - PvP01_06_v2
    - PvP01_07_v2
    - PvP01_08_v2
    - PvP01_09_v2
    - PvP01_10_v2
    - PvP01_11_v2
    - PvP01_12_v2
    - PvP01_13_v2
    - PvP01_14_v2
    - PvP01_API_v2
    - PvP01_MIT_v2

    # this is amplicon sequencing, so turn off deduplication
    deduplicate: False

    # set below to True for new samples (for the purpose of submitting to SRA database)
    keep_paired_bam: False

    # set below to True for structural variant (SV) analysis such as CNV
    keep_final_bam: True

    # filter to apply for deduplicated BAM
    # read_filters: --remove_FR --remove_RF --remove_FF --remove_RR --remove_trans --remove_unmapped --remove_secondary --remove_supplementary
    # for below, we keep FR, RF, FF and RR orientation for future SV analysis
    read_filters: --remove_trans --remove_unmapped --remove_secondary --remove_supplementary


    reads_trimmer_wf: trimmer_null.smk
    reads_mapper_wf: mapper_bwa-mem2.smk

    # multistep-variant-caller parameters
    prepare_sample_directory_flags: --force
    sample_variant_caller_flags: --force --no-config-cascade
    joint_variant_caller_flags: --force --no-config-cascade

    variant_caller_wf: ssf_varcall_gatk_drag.smk

    gatk_calibrate_str: False
    haplotypecaller_flags: --max-reads-per-alignment-start 0 --do-not-run-physical-phasing --pileup-detection --dont-use-soft-clipped-bases

    # --- end of NGS-Pipeline settings ---

    # --- parameters for Microhaplotype settings ---

    class: "parasite"
    maxEE: "5,5"
    trim_right: "10,10"
    min_length: 30 
    truncQ: "5,5"
    max_consist: 10
    omegaA: 1e-120
    justconcat: 0
    platform: "PE"
    trimqv: 15

    #refs

    insertseq: configs/refs/pv-vvg/Microhaps_Inserts_wMito.fasta  # To be modified, see, Reference Files section
    primer_fw: configs/refs/pv-vvg/microhap_pr_fwd.min_overlap.fasta # To be modified, see, Reference Files section
    primer_rev: configs/refs/pv-vvg/microhap_pr_rv.min_overlap.fasta # To be modified, see, Reference Files section

    # EOF


----------

- ``regions``: The chromosomes in the ``refseq_file`` for discovery variant calling.
- ``complete_region``: The complete region for variant calling in the new panel.

An example of the complete config file for a different panel can be found at `configs/refs/pf-spotmal/pfspotmal-mhap.yaml <https://raw.githubusercontent.com/vivaxgen/MicroHaps/refs/heads/main/configs/refs/pf-spotmal/pfspotmal-mhap.yaml>`_.

Reference Files
-----------------
These files are required for the new panel and should be placed in the ``configs/refs/<new_panel>`` directory. All the filepath should be relative to the ``path-to-vvg-Microhaps/envs/MicroHaps/`` directory.

- ``refseq_file``: The reference genome file for the new panel. This file should be in FASTA format and should contain the reference genome for the new panel.
- ``targetregion_file``: The target region file for the new panel. This file should be in BED format and should contain the target regions for the new panel.
- ``insertseq``: The insert sequence file for the new panel. This file should be in FASTA format and should contain the insert sequences for the new panel.
- ``primer_fw``: The forward primer file for the new panel. This file should be in FASTA format and should contain the forward primers for the new panel.
- ``primer_rev``: The reverse primer file for the new panel. This file should be in FASTA format and should contain the reverse primers for the new panel.
