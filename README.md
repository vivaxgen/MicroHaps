# Menzies_MicroHaps
GitHub adaptation of MicroHaplotype pipeline for collaborators.

## Installation
Install pre-requisite packages and repositories
```
conda create -n microhapQC
conda activate microhapQC

conda install python=3.7 bwa samtools bcftools freebayes parallel datamash gatk4=4.1.4.1 delly tqdm trimmomatic minimap2 biopython bedtools r-ggplot2 iqtree fastqc mosdepth samclip sambamba

mkdir tools
cd /tools/
git clone https://github.com/pathogenseq/fastq2matrix.git
python setup.py install

```

Installation of MicroHaps GitHub
```
git clone https://github.com/aosborne13/Menzies_MicroHaps
cd Menzies_MicroHaps
python setup.py install
```
