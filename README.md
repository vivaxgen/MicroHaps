# Menzies_MicroHaps
GitHub adaptation of MicroHaplotype pipeline for collaborators.

## Installation
Install pre-requisite packages and repositories; storing repositories in easily accessible tools folder for quick maintenance.
```
conda create -n microhapQC
conda activate microhapQC

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install python=3.7 bwa samtools bcftools freebayes parallel datamash gatk4=4.1.4.1 delly tqdm trimmomatic minimap2 biopython bedtools r-ggplot2 iqtree fastqc mosdepth samclip sambamba

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
## Run MicroHap Quality Control step
Create sample list CSV file to run script.
```
ls *_R1.fastq.gz | sed 's/.fastq.gz//' > samples.txt
sed -i 's/_R1$//' samples.txt
echo -e "sample" | cat - samples.txt > samples_header.txt
sed 's/ \+/,/g' samples_header.txt > sample_file.csv

usage: microhap_QC.py [-h] --index-file INDEX_FILE --ref REF [--version]

MicroHaplotype Quality Control script

optional arguments:
  -h, --help            show this help message and exit
  --index-file INDEX_FILE
                        CSV file containing field "Sample" (default: None)
  --ref REF             Reference fasta (default: None)
  --version             show program's version number and exit

```
