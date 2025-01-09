# VivaxGEN MicroHaplotypes

This repository contains 2 pipelines for handling microhaplotype amplicon
sequencing data:

* DADA2-based MicroHaplotype Analysis Pipeline; Original software can be found at (https://benjjneb.github.io/dada2/) but is modified for this pipeline. DADA2-based pipeline COMING SOON.

* SNP-based variant calling which output ordinary VCF file, based on
  [vivaxGEN NGS-Pipeline](https://github.com/vivaxgen/ngs-pipeline).


## Installation on local devices and HPCs

Note for Conda-based users!
Be sure you are not in a conda environment or in the (base) conda environment prior to installing. 
To deactivate your conda environment or (base) environment, enter:

	conda deactivate

For vivaxGEN MicroHaps sequencing pipeline, install with the following command:

	"${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/MicroHaps/main/install.sh)

Copy and paste command into the command line in the folder you want the install to be saved to (we recommend you create a specific tools or software folder). 
When prompted for "Pipeline base directory? [./vvg-MicroHaps]" press enter again for the install to proceed.

The installation requires ~ 20-45 minutes as most of R packages need to be recompiled
during installation.

Once the installation finished, it will show the command to activate the
pipeline, as such:

	/path/to/vvg-MicroHaps/bin/activate

This activation command has to be executed before all commands of the pipeline
can be run.
When activated, the terminal will show the **(µhaps)** prompt.

The installation process also performs indexing of the reference files.
However, in case that the indexing fails, please perform manual indexeing
using the command:

	ngs-pl initialize --target wgs

 To test your install, and read about programme specifications / options:

 	ngs-pl run-discovery-variant-caller --help


## Updating the pipeline

To update the pipeline line, assuming that the environment has been activated,
run the following command:

	$VVGBIN/update-pipeline.sh


# Tutorial: How does the Menzies MicroHap pipeline work and what is it doing?

Introduction to SNP-based variant calling
------------

This pipeline contains a wrapper for vivaxGEN NGS-Pipeline, a SNP-based variant
calling pipeline, which has been setup to process P vivax sequence data.
This pipeline will produce ordinary VCF file that can be used for further
downstream analysis.

Currently, this pipeline was setup to use multi-step mode of vivaxGEN
NGS-Pipeline in a single command line with GATK-based workflow to generate VCF
file
There is also Freebayes-based workflow, but it requires modification of the
setting, which will not be covered in this documentation.

Quick Tutorial
--------------

The following instructions show how to run this pipeline.
It is assumed that the pipeline has been installed.

1.  Activate the environment. The terminal will show ``(µhaps)`` prompt.
```
	/path/to/vvg-MicroHaps/bin/activate
```
2.  Go to the directory that will be used to process and analysis::
```
	cd MY_ANALYSIS_DIRECTORY
```
3.  Provide the FASTQ reads from the sequencing result, by either copying the
    FASTQ files or alternatively generate soft link as necessary.
    The soft link approach is preferred since it will prevent duplication of
    the files, if the files are already reside in the local storage.
    The FASTQ files should be in fastq.gz format (gzip-compressed), and the
    filenames should reflect the sample name, eg: my-sample-01_R1.fastq.gz.
    If the FASTQ filenames contains something that are not part of the sample
    names, there is an option ``--underscore`` in the next step that can be
    used::

    	mkdir reads
    	cp /path/to/fastq/files/*.fastq.gz .

    or::

    	ln -s SOME_SOURCE_DIR reads

Run command:
```
	ngs-pl run-discovery-variant-caller -o MY_OUTPUT reads/*.fastq.gz
```
example for paired-end (Illumina) reads::
```
	ngs-pl run-discovery-variant-caller -j 1 -o output_dir -u 4 --paired  reads/*.fastq.gz
```
LAPTOP USERS - It is essential that you specify "-j 1" as this limits the number of jobs running at one time. Without this argument, 
the pipeline will utilise too much system memory and crash. We are working on improving this issue.

The "-u 4" prompt changes depending on how your files are named as it counts the number of underscores "_" in a file name. For a sample called "sample4_date_batch_info_R1.fastq.gz",
the "-u 4" argument will use result in "sample4" as the ID. For a sample called "sample_4_date_batch_pool_country_R1.fastq.gz", the argument "-u 5" will result in "sample_4" as the ID.

4. When the command finishes, examine the content of ``MY_OUTPUT`` directory::
```
	cd MY_OUTPUT
	ls
```
The layout of the output directory is::

    MY_OUTPUT/
              metafile/
                       manifest.tsv
              analysis/
                       SAMPLE-1/
                       SAMPLE-2/
                       ...
              joint/
                    vcfs/

------------
Introduction to DADA2-based MicroHaplotype software
------------
Introduction to Microhaplotype based analysis
COMING SOON.
This pipeline contains a wrapper for the MIT-Broad team DADA2 software. (https://github.com/broadinstitute/malaria-amplicon-pipeline) 

Quick Tutorial
--------------

Prepare pipeline install, analysis directory and raw data as above. 

1.  Execute the ``ngs-pl run-full-analysis`` command with as follow::
```
	ngs-pl run-full-analysis -o outdir test-data/*.fastq.gz
```
2. When the command finishes, examine the content of ``outdir`` directory
```
    outdir/
            alignments/
                marker_1.fasta
                marker_1.msa
                marker_2.fasta
                marker_2.msa
                ...
            trimmed/
                sample_1_R1.trimmed.fastq.gz
                sample_1_R2.trimmed.fastq.gz
                ...
            malamp/
                dada2/
                    ...
                ASVSeqs.fasta
                ASVTable.txt
                asv_to_cigar
                outputCIGAR.tsv
                meta
```

The primary output file of interest is the `outputCIGAR.tsv` which contains the haplotype and their frequencies across the samples.

------------

