rule qc_haplotypes:
    localrule: True
    input:
        tab = f"{outdir}/malamp/outputHaplotypes.tsv" if new_postprocess != "old" else f"{outdir}/malamp/outputCIGAR.tsv",
    output:
        depths = f"{outdir}/malamp/depths.tsv"
    params:
        mindepth = 10,
        outdir = lambda w, output: output.depths.rsplit('/', 1)[0]
    shell:
        "ngs-pl tab-to-QC -d {params.mindepth} --outdir {params.outdir} {input}"


rule qc_heatmap:
    input:
        depth = "{pfx}/malamp/depths.tsv",
    output:
        depth = "{pfx}/malamp/depths-microhaps.png",
        marker = "{pfx}/malamp/markers-microhaps.png",
    shell:
        "ngs-pl tab-to-plots --index-column Amplicon_name --additional-title 'Microhaps'"
        "  --outheatmap {output.depth} --outmarker {output.marker} {input.depth}"


rule depth_ratio:
    localrule: True
    input:
        depth_mapped = "{pfx}/depths-mapped.tsv",
        depth_microhaps = "{pfx}/malamp/depths.tsv",
    output:
        ratio = "{pfx}/malamp/depth-ratio.tsv"
    run:

        import pandas as pd

        # Read the depth files
        mapped_df = pd.read_csv(input.depth_mapped, sep='\t')
        microhaps_df = pd.read_csv(input.depth_microhaps, sep='\t')

        depth_mapped = mapped_df.iloc[:, 4:].set_index(mapped_df.Amplicon_name)
        depth_microhaps = microhaps_df.iloc[:, 4:].set_index(microhaps_df.Amplicon_name)

        ratio_df = (depth_microhaps / depth_mapped).fillna(0).reindex(depth_mapped.index).reset_index()

        ratio_df.to_csv(output.ratio, sep='\t', index=False)


rule depth_ratio_heatmap:
    localrule: True
    input:
        ratio = "{pfx}/depth-ratio.tsv",
    output:
        heatmap = "{pfx}/depth-ratio-heatmap.png",
        marker = "{pfx}/depth-ratio-markers.png",
    shell:
        "ngs-pl tab-to-plots --value-is-ratio --index-column Amplicon_name --additional-title 'Microhaps Ratio'"
        "  --outheatmap {output.heatmap} --outmarker {output.marker} {input.ratio}"