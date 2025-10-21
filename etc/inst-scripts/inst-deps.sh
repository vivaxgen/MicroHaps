
#echo "Installing latest bedtools"
#micromamba -y install bedtools -c conda-forge -c bioconda -c defaults

#echo "Installing datamash"
#micromamba -y install datamash -c conda-forge -c bioconda -c defaults

#echo "Installing trimmomatic"
#micromamba -y install trimmomatic -c conda-forge -c bioconda -c defaults

#echo "Installing mosdepth"
#micromamba -y install mosdepth -c conda-forge -c bioconda -c defaults

#echo "Installing samclip"
#micromamba -y install samclip -c conda-forge -c bioconda -c defaults

echo "Installing other dependencies with micromamba"
retry 5 micromamba -y install -n ${uMAMBA_ENVNAME} -f ${ENVS_DIR}/MicroHaps/etc/inst-scripts/env.yaml
retry 5 micromamba create -y -n ${uMAMBA_ENVNAME}-dada2 

echo "Reloading profiles"
reload_vvg_profiles

echo "Initialize enviroment"
ngs-pl initialize --panel pvvvg-mhap --target wgs
ngs-pl initialize --panel pfspotmal-mhap --target wgs

echo "Installing additional R packages"
retry 5 micromamba install -y -n ${uMAMBA_ENVNAME}-dada2  -f ${ENVS_DIR}/MicroHaps/etc/inst-scripts/dada2.yaml
retry 5 micromamba -y remove -n ${uMAMBA_ENVNAME}-dada2 --no-prune-deps bioconductor-dada2
micromamba run -n ${uMAMBA_ENVNAME}-dada2 Rscript -e 'if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager", repos = "https://cloud.r-project.org")}; BiocManager::install(version = "3.20", ask= FALSE, force=TRUE); BiocManager::install(c("dada2", "limma"));  if (!require("dcifer", quietly = TRUE)){ install.packages("dcifer", repos = "https://cloud.r-project.org")}; if (!require("moire", quietly = TRUE)){install.packages("moire", repos = c("https://eppicenter.r-universe.dev", "https://cloud.r-project.org"))}'