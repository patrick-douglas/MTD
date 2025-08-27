#!/bin/bash
# Check R package versions in a nice table for MTD

# Cores
w=$(tput sgr0)       # branco
g=$(tput setaf 2)    # verde
r=$(tput setaf 1)    # vermelho

# Versão do R
R_ver=$(R --version | grep version | grep R | awk '{print $3}')

# Pacotes a checar
pkgs=( tidyverse car rstatix ggpubr plyr textshaping ragg BiocManager rlang vctrs BiocGenerics S4Vectors IRanges UCSC.utils GenomeInfoDbData GenomeInfoDb matrixStats formatR lambda.r futile.options futile.logger RColorBrewer )

# Cabeçalho
echo "${g}╔═══════════════════╦════════════════╗"
printf "║ R                 ║ ${w}%-14s${g} ║\n" "$R_ver"
printf "║ Conda Environment ║ ${w}%-14s${g} ║\n" "MTD"
echo "╠═══════════════════╬════════════════╣"

# Loop pelos pacotes
for pkg in "${pkgs[@]}"; do
    ver=$(R --no-restore -e "packageVersion(\"$pkg\")" 2>/dev/null | grep '\[1\]' | awk '{print $2}' | sed -r 's/^.{1}//; s/.$//')
    if [ -z "$ver" ]; then
        ver="not installed"
        color=$r
    else
        color=$w
    fi
    printf "║ %-17s ║ ${color}%-14s${g} ║\n" "$pkg" "$ver"
done

# Rodapé
echo "╚═══════════════════╩════════════════╝"
