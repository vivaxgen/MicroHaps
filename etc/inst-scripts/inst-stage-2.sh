
echo -e "\e[32m>>>> Linking resource files\e[0m"
${VVGBIN}/link-resource-files.sh ${ENVS_DIR}/MicroHaps/etc/bashrc.d

echo -e "\e[32m>>>> Reloading profiles\e[0m"
reload_vvg_profiles

ENVS_DIR="${ENVS_DIR:-${VVG_BASEDIR}/envs}"
INST_SCRIPTS_DIR="${ENVS_DIR}/MicroHaps/etc/inst-scripts"

if [[ -z ${VVG_MANIFEST_FILE:-} ]]; then
  echo -e "\e[32m>>> No manifest file provided, installing dependencies with inst-deps.sh\e[0m"
  source ${INST_SCRIPTS_DIR}/inst-deps.sh
fi

echo -e "\e[32m>>>> Installing additional R packages\e[0m"
#retry 5 micromamba -y remove -n ${uMAMBA_ENVNAME} --no-prune-deps bioconductor-dada2
#Rscript -e 'if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager", repos = "https://cloud.r-project.org")}; BiocManager::install(version = "3.20", ask= FALSE, force=TRUE); BiocManager::install(c("dada2", "limma"));  if (!require("dcifer", quietly = TRUE)){ install.packages("dcifer", repos = "https://cloud.r-project.org")}; if (!require("moire", quietly = TRUE)){install.packages("moire", repos = c("https://eppicenter.r-universe.dev", "https://cloud.r-project.org"))}'
echo -e "\e[32m>>>> Installing R dcifer\e[0m"
Rscript -e 'if (!require("dcifer", quietly = TRUE)){ install.packages("dcifer", repos = "https://cloud.r-project.org")};'
echo -e "\e[32m>>>> Installing R moire\e[0m"
Rscript -e 'if (!require("moire", quietly = TRUE)){install.packages("moire", repos = c("https://eppicenter.r-universe.dev", "https://cloud.r-project.org"))}'

echo -e "\e[32m>>>> Initialize panels\e[0m"
ngs-pl initialize --panel pvvvg-mhap --target wgs
ngs-pl initialize --panel pf-m4h --target wgs
ngs-pl initialize --panel pfspotmal-mhap --target wgs

