
echo ">>> Linking resource files"
${VVGBIN}/link-resource-files.sh ${ENVS_DIR}/MicroHaps/etc/bashrc.d

echo ">>> Reloading profiles"
reload_vvg_profiles

ENVS_DIR="${ENVS_DIR:-${VVG_BASEDIR}/envs}"
INST_SCRIPTS_DIR="${ENVS_DIR}/MicroHaps/etc/inst-scripts"

echo ">>> Installing generic global dependencies"
pixi-global-install ${INST_SCRIPTS_DIR}/global-generics.spec

echo ">>> Installing generic workspace dependencies (Python and R)"
set +u
pixi-add ${INST_SCRIPTS_DIR}/workspace-generics.spec
set -u

echo ">>> Installing additional R packages"
#retry 5 micromamba -y remove -n ${uMAMBA_ENVNAME} --no-prune-deps bioconductor-dada2
#Rscript -e 'if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager", repos = "https://cloud.r-project.org")}; BiocManager::install(version = "3.20", ask= FALSE, force=TRUE); BiocManager::install(c("dada2", "limma"));  if (!require("dcifer", quietly = TRUE)){ install.packages("dcifer", repos = "https://cloud.r-project.org")}; if (!require("moire", quietly = TRUE)){install.packages("moire", repos = c("https://eppicenter.r-universe.dev", "https://cloud.r-project.org"))}'
echo ">>> Installing R dcifer"
Rscript -e 'if (!require("dcifer", quietly = TRUE)){ install.packages("dcifer", repos = "https://cloud.r-project.org")};'
echo ">>> Installing R moire"
Rscript -e 'if (!require("moire", quietly = TRUE)){install.packages("moire", repos = c("https://eppicenter.r-universe.dev", "https://cloud.r-project.org"))}'

echo ">>> Initialize panels"
ngs-pl initialize --panel pvvvg-mhap --target wgs
ngs-pl initialize --panel pf-m4h --target wgs
ngs-pl initialize --panel pfspotmal-mhap --target wgs

# EOF

