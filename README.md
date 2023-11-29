# Menzies_MicroHaps
GitHub adaptation of MicroHaplotype pipeline.

## Tutorial: How does the Menzies MicroHap pipeline work and what is it doing?
COMING SOON

## Installation on local devices and private servers (Skip this if using ADA)
Create conda environment with required packages by downloading the microhapQC.yaml file provided in input files (https://github.com/aosborne13/Menzies_MicroHaps/blob/main/input_files/microhapQC.yaml) and following the commands below. 

Manual installation is found at the bottom of this page if required or if individual packages are missing/incorrectly installed.
```
mkdir tools
cd /tools/
conda env create -f microhapQC.yaml
conda activate microhapQC
```
R packages managed by BiocManager are not currently included in the conda environment and require manual installation.

Install DADA2 into the R client while the microhapQC conda environment is active.
```
# Open R in command line
R

#install DADA2 and pre-requisites using BiocManager
BiocManager::install("GenomeInfoDb")
BiocManager::install("GenomicRanges")
BiocManager::install("Biostrings")
BiocManager::install("Rsamtools")
BiocManager::install("SummarizedExperiment")
BiocManager::install("GenomicAlignments")
BiocManager::install("ShortRead")
BiocManager::install("dada2")
BiocManager::install("limma")

# Quit R and do not save current workspace using 'n'
q()
n
```

Install pre-requisite GitHub repositories; store repositories in easily accessible "tools" folder for quick maintenance.
```
git clone https://github.com/pathogenseq/fastq2matrix.git
cd fastq2matrix
python setup.py install

cd ..
git clone https://github.com/aosborne13/Menzies_MicroHaps
cd Menzies_MicroHaps
python setup.py install
```
## Running only the FASTQ Quality Control step (on local devices/private servers)
Create sample list CSV file, using the command below, to run script in the folder containing your FASTQ files. 

Alternatively, you can manually create a CSV sample file with only the samples you require. The CSV needs the sample IDs (which should correspond to your FASTQ file IDs) in a single column with "sample" as the column name. The column name "sample" is case sensitive.
```
ls *_R1.fastq.gz | sed 's/.fastq.gz//' | sed 's/_R1$//' | (echo "sample" && cat -) | sed 's/ \+/,/g' > sample_file.csv
```
The MicroHaplotype Quality Control script needs to be run in the same directory as your FASTQ files. Remember to create a back-up copy of your raw FASTQ files elsewhere, in case of error.
```
usage: microhap_QC.py [-h] --index-file INDEX_FILE --ref REF --bed BED [--version]

MicroHaplotype Quality Control script

optional arguments:
  -h, --help            show this help message and exit
  --index-file          INDEX_FILE CSV file containing field "Sample" (default: None)
  --ref REF             Reference fasta (default: None)
  --bed BED             BED file with MicroHaplotype locations (default: None)
  --version             show program's version number and exit

```
Example of usage with your input files stored in a separate directory. The command is run in the directory containing the FASTQ files listed in the CSV file.
```
microhap_QC.py --index-file ~/Documents/microhaps/sample_file.csv --ref ~/Documents/microhaps/PlasmoDB-51_PvivaxP01_Genome.fasta --bed ~/Documents/microhaps/microhap.bed
```
## Running the complete MicroHaplotype pipeline on local devices and private servers (not ADA)
POSTPROC_DADA2 AND ASV_TO_CIGAR COMING SOON.

Running the MicroHaplotype pipeline carries out both quality control of raw read data, as well as downstream processing, including DADA 2.

Run the below command to generate a sample file for the pipeline in the sirectory containing your raw FASTQ reads. 

Alternatively, you can manually create a CSV sample file with only the samples you require. The CSV needs the sample IDs (which should correspond to your FASTQ file IDs) in a single column with "sample" as the column name. The column name "sample" is case sensitive.
```
ls *_R1.fastq.gz | sed 's/.fastq.gz//' | sed 's/_R1$//' | (echo "sample" && cat -) | sed 's/ \+/,/g' > sample_file.csv
```
Run the pipeline using the sample_file.csv as the INDEX_FILE. You can store the index file elsewhere but the FASTQ files should be located in the directory where you run the pipeline.
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
## Running the MicroHaplotype pipeline on ADA
Running the MicroHaplotype pipeline carries out both quality control of raw read data, as well as downstream processing, including DADA 2. Running on ADA using JSON inputs to submit a patch job for processing.

COMING SOON

## Manual Conda environment setup and R package installation
Create conda environment to store required packages. Conda channel configuration is shown in instructions for first time users. If your conda is already configured, please skip those steps.
```
conda create -n microhapQC
conda activate microhapQC

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install python=3.8 bwa samtools bcftools freebayes parallel datamash gatk4=4.1.4.1 delly tqdm trimmomatic minimap2 biopython bedtools r-ggplot2 iqtree fastqc mosdepth samclip sambamba multiqc pandas cutadapt r-BiocManager r-RCurl r-argparse r-data.table

```
Install DADA2 into R Client
```
# Open R in command line
R

#install DADA2 and pre-requisites using BiocManager
BiocManager::install("GenomeInfoDb")
BiocManager::install("GenomicRanges")
BiocManager::install("Biostrings")
BiocManager::install("Rsamtools")
BiocManager::install("SummarizedExperiment")
BiocManager::install("GenomicAlignments")
BiocManager::install("ShortRead")
BiocManager::install("dada2")
BiocManager::install("limma")

# Quit R and do not save current workspace using 'n'
q()
n
```
