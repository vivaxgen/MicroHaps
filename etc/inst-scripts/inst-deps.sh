


echo -e "\e[32m>>>> Installing generic global dependencies\e[0m"
pixi-global-install ${INST_SCRIPTS_DIR}/global-generics.spec

echo -e "\e[32m>>>> Installing generic workspace dependencies (Python and R)\e[0m"
set +u
pixi-add ${INST_SCRIPTS_DIR}/workspace-generics.spec
set -u


# EOF

