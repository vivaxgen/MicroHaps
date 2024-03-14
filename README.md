# VivaxGEN MicroHaplotypes

This repository contains 2 pipelines for handling microhaplotype amplicon
sequencing data:

* DADA2-based MicroHaplotype Analysis Pipeline; Original software can be found at (https://benjjneb.github.io/dada2/) but is modified for this pipeline. DADA2-based pipeline COMING SOON.

* SNP-based variant calling which output ordinary VCF file, based on
  [vivaxGEN NGS-Pipeline](https://github.com/vivaxgen/ngs-pipeline).

## Installation on local devices and HPCs
Be sure you are not in a conda environemt or in the (base) conda environment prior to installing. 
To deactivate your conda environment or (base) environment, enter:

	conda deactivate

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

Before running, the reference files need to be indexed using the command:

	ngs-pl initialize --target wgs

 To test your install, and read about programme specifications / options:

 	ngs-pl run-discovery-variant-caller --help

## Updating the pipeline

To update the pipeline line, assuming that the environment has been activated,
run the following command:

```
update-pipeline.sh
```

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

		source /path/to/vvg-MicroHaps/bin/activate

3.  Go to the directory that will be used to process and analysis::

		cd MY_ANALYSIS_DIRECTORY

4.  Provide the FASTQ reads from the sequencing result, by either copying the
    FASTQ files or alternatively generate soft link as necessary.
    The soft link approach is preferred since it will prevent duplication of
    the files, if the files are already reside in the local storage.
    The FASTQ files should be in fastq.gz format (gzip-compressed), and the
    filenames should reflect the sample name, eg: my-sample-01_R1.fastq.gz.
    If the FASTQ filenames contains something that are not part of the sample
    names, there is an option ``--underscore`` in the next step that can be
    used::

    	mkdir reads
    	cp /path/to/fastq/files/*.fastq.gz

    or::

    	ln -s SOME_SOURCE_DIR reads

5.  Execute the ``run-discovery-variant-caller`` command with as follow::

		ngs-pl run-discovery-variant-caller -o MY_OUTPUT reads/*.fastq.gz

	example for paired-end (Illumina) reads::

  		ngs-pl run-discovery-variant-caller -o output_dir -u 4 --paired  reads/*.fastq.gz

 	The "-u 4" prompt changes depending on how your files are named as it counts the number of underscores "_" in a file name. For a sample called "sample4_date_batch_info_R1.fastq.gz",
	the "-u 4" argument will use result in "sample4" as the ID. For a sample called "sample_4_date_batch_pool_country_R1.fastq.gz", the argument "-u 5" will result in "sample_4" as the ID.

7. When the command finishes, examine the content of ``MY_OUTPUT`` directory::

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
------------
Introduction to DADA2-based MicroHaplotype software
------------
COMING SOON.

This pipeline contains a wrapper for the MIT-Broad team DADA2 software.

------------

# DEPRECIATED - MicroHaplotype python-wrapped Quality Control step and DADA2 pipeline.

Bare-bones python-wrapped script for small datasets.

------------

## Running only the FASTQ Quality Control step on local devices and HPCs
Create sample list CSV file, using the command below, to run script in the folder containing your FASTQ files. 

Alternatively, you can manually create a CSV sample file with only the samples you require. The CSV needs the sample IDs (which should correspond to your FASTQ file IDs) in a single column with "sample" as the column name. The column name "sample" is case sensitive.
```
ls *_R1.fastq.gz | sed 's/.fastq.gz//' | sed 's/_R1$//' | (echo "sample" && cat -) | sed 's/ \+/,/g' > sample_file.csv
```
The MicroHaplotype Quality Control script needs to be run in the same directory as your FASTQ files. Remember to create a separate copy of your raw FASTQ files elsewhere. 

DO NOT RUN WITHOUT MAKING A BACKUP COPY OF YOUR RAW FASTQ FILES STORED ELSEWHERE.
```
usage: microhap_QC.py [-h] --index-file INDEX_FILE --ref REF --bed BED [--version]

MicroHaplotype Quality Control script

arguments:
  -h, --help            show this help message and exit
  --index-file          INDEX_FILE CSV file containing field "Sample" (default: None)
  --ref REF             Reference fasta (default: None)
  --bed BED             BED file with MicroHaplotype locations (default: None)
  --version             show program's version number and exit

```
Example of usage with your input files stored in a separate directory. The command is run in the directory containing copies of the FASTQ files listed in the CSV file.
```
microhap_QC.py --index-file ~/Documents/microhaps/sample_file.csv --ref ~/Documents/microhaps/PlasmoDB-51_PvivaxP01_Genome.fasta --bed ~/Documents/microhaps/microhap.bed
```
## Running the complete MicroHaplotype pipeline on local devices and HPCs
POSTPROC_DADA2 AND ASV_TO_CIGAR COMING SOON.

Running the MicroHaplotype pipeline carries out both quality control of raw read data, as well as downstream processing, including DADA 2.

Run the below command to generate a sample file for the pipeline in the sirectory containing your raw FASTQ reads. 

Alternatively, you can manually create a CSV sample file with only the samples you require. The CSV needs the sample IDs (which should correspond to your FASTQ file IDs) in a single column with "sample" as the column name. The column name "sample" is case sensitive.
```
ls *_R1.fastq.gz | sed 's/.fastq.gz//' | sed 's/_R1$//' | (echo "sample" && cat -) | sed 's/ \+/,/g' > sample_file.csv
```
Run the pipeline using the sample_file.csv as the INDEX_FILE. You can store the index file elsewhere but the FASTQ files should be located in the directory where you run the pipeline. Remember to create a separate copy of your raw FASTQ files elsewhere. 

DO NOT RUN WITHOUT MAKING A BACKUP COPY OF YOUR RAW FASTQ FILES STORED ELSEWHERE.
```
usage: microhap_pipeline_beta.py [-h] --index-file INDEX_FILE --ref REF --bed BED [--trim]
                                 [--trim-qv TRIM_QV] --output_file OUTPUT_FILE --pattern_fw
                                 PATTERN_FW --pattern_rv PATTERN_RV --pr1 PR1 --pr2 PR2
                                 [--Class CLASS] [--maxEE MAXEE] [--trimRight TRIMRIGHT]
                                 [--minLen MINLEN] [--truncQ TRUNCQ] [--max_consist MAX_CONSIST]
                                 [--omegaA OMEGAA] [--justConcatenate JUSTCONCATENATE]
                                 [--saveRdata SAVERDATA] [--version]

MicroHaplotype Pipeline
  -h, --help            show this help message and exit

required arguments:
  --index-file          INDEX_FILE
                        CSV file containing field "Sample" (default: None)
  --ref REF             Reference fasta (default: None)
  --bed BED             BED file with MicroHaplotype locations (default: None)
  --trim                Perform triming (default: False)
  --output_file         OUTPUT_FILE
                        Output meta file; to be used in path to meta (default: None)
  --pattern_fw          PATTERN_FW
                        Pattern for forward reads, e.g. "*_R1.fastq.gz" (default: None)
  --pattern_rv          PATTERN_RV
                        Pattern for reverse reads, e.g. "*_R2.fastq.gz" (default: None)
  --pr1 PR1             Path to forward primers FASTA file (default: None)
  --pr2 PR2             Path to reverse primers FASTA file (default: None)

optional arguments:
  --trim-qv             TRIM_QV
                        Quality value to use in the sliding window analysis (default: 5)
  --Class CLASS         Specify Analysis class. Accepts one of two: parasite/vector (default:
                        parasite)
  --maxEE MAXEE         Maximum Expected errors (dada2 filtering argument) (default: 5,5)
  --trimRight           TRIMRIGHT
                        Hard trim number of bases at 5` end (dada2 filtering argument) (default:
                        10,10)
  --minLen MINLEN       Minimum length filter (dada2 filtering argument) (default: 30)
  --truncQ TRUNCQ       Soft trim bases based on quality (dada2 filtering argument) (default: 5,5)
  --max_consist         MAX_CONSIST
                        Number of cycles for consistency in error model (dada2 argument) (default:
                        10)
  --omegaA OMEGAA       p-value for the partitioning algorithm (dada2 argument) (default: 1e-120)
  --justConcatenate     JUSTCONCATENATE
                        whether reads should be concatenated with N's during merge (dada2
                        argument) (default: 0)
  --version             show program's version number and exit
```

Example Usage:
```
microhap_pipeline_beta.py --index-file ~/Documents/microhaps/sample_file.csv --ref ~/Documents/microhaps/PlasmoDB-51_PvivaxP01_Genome.fasta --bed ~/Documents/microhaps/microhap.bed --trim --output_file meta_file --pattern_fw "*_R1.trimmed.fastq.gz" --pattern_rv "*_R2.trimmed.fastq.gz" --pr1 ~/Documents/microhaps/microhap_pr_fwd.min_overlap.fasta --pr2 ~/Documents/microhaps/microhap_pr_rv.min_overlap.fasta
```
