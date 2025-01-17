
echo "Installing latest bedtools"
micromamba -y install bedtools -c conda-forge -c bioconda -c defaults

echo "Installing datamash"
micromamba -y install datamash -c conda-forge -c bioconda -c defaults

echo "Installing trimmomatic"
micromamba -y install trimmomatic -c conda-forge -c bioconda -c defaults

echo "Installing mosdepth"
micromamba -y install mosdepth -c conda-forge -c bioconda -c defaults

echo "Installing samclip"
micromamba -y install samclip -c conda-forge -c bioconda -c defaults

echo "Installing muscle version v3.8.1551"
micromamba -y install "muscle=3.8.1551" -c conda-forge -c bioconda -c defaults

echo "Installing required R packages"
micromamba -y install r-ggplot2 r-BiocManager r-RCurl r-argparse r-data.table r-seqinr r-doMC -c conda-forge -c bioconda -c defaults

echo "Installing ivar"
micromamba -y install ivar -c bioconda

echo "Installing dada2 amd limma"
micromamba -y install bioconductor-dada2 bioconductor-limma -c conda-forge -c bioconda

echo "Installing bbtools"
micromamba -y install bbmap -c conda-forge -c bioconda -c defaults

echo "installing required Python modules"

# to use latest of all python-related stuff, uncomment below and remove the conda parts
pip3 install biopython
pip3 install cutadapt
pip3 install tqdm

echo "Reloading profiles"
reload_vvg_profiles

echo "Initialize enviroment"
ngs-pl initialize --target wgs

# EOF
