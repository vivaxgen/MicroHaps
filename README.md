# Menzies_MicroHaps
GitHub adaptation of MicroHaplotype pipeline for collaborators.

## Installation
Create conda environment to store required packages. Conda channel configuration is shown in instructions for first time users. If your conda is already configured, please skip those steps.
```
conda create -n microhapQC
conda activate microhapQC

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install python=3.8 bwa samtools bcftools freebayes parallel datamash gatk4=4.1.4.1 delly tqdm trimmomatic minimap2 biopython bedtools r-ggplot2 iqtree fastqc mosdepth samclip sambamba multiqc pandas
```
Install pre-requisite packages and repositories; storing repositories in easily accessible "tools" folder for quick maintenance.
```
mkdir tools
cd /tools/
git clone https://github.com/pathogenseq/fastq2matrix.git
cd fastq2matrix
python setup.py install

cd ..
git clone https://github.com/aosborne13/Menzies_MicroHaps
cd Menzies_MicroHaps
python setup.py install
```
## Run MicroHap pipeline
Running the MicroHap pipeline carries out both quality control of raw read data, as well as downstream processing, including DADA 2.

DOWNSTREAM PROCESSING NOT ADDED YET - COMING SOON
```
usage: microhap_pipeline.py [-h] --index-file INDEX_FILE --ref REF --bed BED [--trim]
                                    [--trim-qv TRIM_QV] --path_to_fq PATH_TO_FQ --output_file
                                    OUTPUT_FILE --pattern_fw PATTERN_FW --pattern_rv PATTERN_RV
                                    [--version]

MicroHaplotype Pipeline

optional arguments:
  -h, --help            show this help message and exit
  --index-file          INDEX_FILE
                        CSV file containing field "Sample" (default: None)
  --ref REF             Reference fasta (default: None)
  --bed BED             BED file with MicroHaplotype locations (default: None)
  --trim                Perform triming (default: False)
  --trim-qv TRIM_QV     Quality value to use in the sliding window analysis (default: 5)
  --path_to_fq PATH_TO_FQ
                        Path to fastq files (default: None)
  --output_file OUTPUT_FILE
                        Output meta file (default: None)
  --pattern_fw PATTERN_FW
                        Pattern for forward reads, e.g. "*_R1.fastq.gz" (default: None)
  --pattern_rv PATTERN_RV
                        Pattern for reverse reads, e.g. "*_R2.fastq.gz" (default: None)
  --version             show program's version number and exit
```

Example Usage:
```
microhap_pipeline.py --index-file sample_file.csv --ref PlasmoDB-51_PvivaxP01_Genome.fasta --bed microhap.bed --trim --path_to_fq /path/to/fastq/FASTQ_FILES --output_file meta_file --pattern_fw "*_R1.fastq.gz" --pattern_rv "*_R2.fastq.gz"
```
## Run MicroHap Quality Control ONLY
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
