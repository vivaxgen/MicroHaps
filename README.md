# Menzies_MicroHaps

This repository contains 2 pipelines for handling microhaplotype amplicon
sequencing data:

* DADA2-based MicroHaplotype Analysis Pipeline [DADA2] (https://benjjneb.github.io/dada2/).

* SNP-based variant calling which output ordinary VCF file, based on
  [vivaxGEN NGS-Pipeline](https://github.com/vivaxgen/ngs-pipeline).

## Installation on local devices and HPCs

For vivaxGEN MicroHaps sequencing pipeline, install with the following command:

	"${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/install/main/MicroHaps-pl.sh)

Copy and paste command into the command line in the folder you want the install to be saved to (we recommend you create a specific tools or software folder). 
When promted for "Pipeline base directory? [./vvg-MicroHaps]" press enter again for the install to proceed.

The installation requires ~20 minutes as most of R packages need to be recompiled
during installation.

Once the installation finished, it will show the command to activate the
pipeline. This activation command has to be executed before all commands of
the pipeline can be run. When activated, the terminal will show the **(µhaps)**
prompt.


## Updating the pipeline

To update the pipeline line, assuming that the environment has been activated,
run the following command:

```
update-pipeline.sh
```

# Tutorial: How does the Menzies MicroHap pipeline work and what is it doing?
Introduction
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

2.  Go to the directory that will be used to process and analysis::

		cd MY_ANALYSIS_DIRECTORY

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
    	cp SOME_WHERE/*.fastq.gz

    or::

    	ln -s SOME_SOURCE_DIR reads

4.  Execute the ``run-discovery-variant-caller`` command with as follow::

		ngs-pl run-discovery-variant-caller -o MY_OUTPUT reads/*.fastq.gz

5. When the command finishes, examine the content of ``MY_OUTPUT`` directory::

		cd MY_OUTPUT
		ls

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

