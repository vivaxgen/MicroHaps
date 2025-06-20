
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

echo "Reloading profiles"
reload_vvg_profiles

echo "Initialize enviroment"
ngs-pl initialize --panel pvvvg-mhap --target wgs
ngs-pl initialize --panel pfspotmal-mhap --target wgs

echo "Installing additional R packages"
Rscript -e 'if (!require("dcifer", quietly = TRUE)){ install.packages("dcifer", repos = "https://cloud.r-project.org")}; if (!require("moire", quietly = TRUE)){install.packages("moire", repos = c("https://eppicenter.r-universe.dev", "https://cloud.r-project.org"))}'
# EOF
