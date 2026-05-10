
INST_SCRIPTS_DIR="${ENVS_DIR}/ngs-pipeline/etc/inst-scripts"

echo ">>>Installing generic global dependencies"
pixi-global-install ${INST_SCRIPTS_DIR}/global-generics.spec

echo ">>>Installing generic workspace dependencies (Python and R)"
pixi-add ${INST_SCRIPTS_DIR}/workspace-generics.spec

echo ">>>Reloading profiles"
reload_vvg_profiles

echo ">>>Initialize enviroment"
ngs-pl initialize --panel pvvvg-mhap --target wgs
ngs-pl initialize --panel pf-m4h --target wgs
ngs-pl initialize --panel pfspotmal-mhap --target wgs

echo ">>>Installing additional R packages"
#retry 5 micromamba -y remove -n ${uMAMBA_ENVNAME} --no-prune-deps bioconductor-dada2
Rscript -e 'if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager", repos = "https://cloud.r-project.org")}; BiocManager::install(version = "3.20", ask= FALSE, force=TRUE); BiocManager::install(c("dada2", "limma"));  if (!require("dcifer", quietly = TRUE)){ install.packages("dcifer", repos = "https://cloud.r-project.org")}; if (!require("moire", quietly = TRUE)){install.packages("moire", repos = c("https://eppicenter.r-universe.dev", "https://cloud.r-project.org"))}'
