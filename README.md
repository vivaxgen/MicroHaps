# VivaxGEN MicroHaplotypes

This repository hosts vivaxGEN-Microhaps, an open-source pipeline designed for processing targeted amplicon sequencing data. It offers two distinct sub-pipelines:

* **Microhaplotype Calling:** This sub-pipeline, adapted from the
  [malaria-amplicon-pipeline](https://github.com/broadinstitute/malaria-amplicon-pipeline/) and converted to use Snakemake, generates microhaplotype allele tables.

* **SNP-based Variant Calling:** Leveraging the [vivaxGEN NGS-Pipeline](https://github.com/vivaxgen/ngs-pipeline), this sub-pipeline generates standard VCF files for SNP analysis.


## Documentation

The main documentation, which cover complete installation and tutorials, is [here](https://vivaxgen-microhaps.readthedocs.io/en/latest/).


## Quick installation on local devices and HPCs

Note for Conda-based users!
Be sure you are not in a conda environment or in the (base) conda environment prior to installing. 
To deactivate your conda environment or (base) environment, enter:

	conda deactivate

Install the pipeline with the following command (note: must be run under
relatively current version of bash):

	"${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/MicroHaps/main/install.sh)


The installation requires ~ 20-45 minutes as most of R packages need to be recompiled
during installation.

Once the installation finished, it will show the command to activate the
pipeline, such as:

	/path/to/vvg-MicroHaps/bin/activate

This activation command has to be executed before all commands of the pipeline
can be run.
When activated, the terminal will show the **(µhaps)** prompt.

The installation process also performs indexing of the reference files.
However, in case that the indexing fails, please perform manual indexeing
using the command:

	ngs-pl initialize --target wgs

 To test your install, and read about programme specifications / options:

 	ngs-pl run-microhaplotype-caller --help


## Updating the pipeline

To update the pipeline line, assuming that the environment has been activated,
run the following command:

	$VVGBIN/update-box

------------

