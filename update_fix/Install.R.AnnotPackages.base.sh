#!/bin/bash

# Install all required dependencies for running create_annotation_package.R

echo "📦 Checking if R is installed..."
if ! command -v R &> /dev/null
then
    echo "❌ R is not installed. Installing now..."
    sudo apt update
    sudo apt install -y r-base
else
    echo "✅ R is already installed."
fi

echo "📁 Creating temporary R script to install required packages..."

cat << 'EOF' > install_R_dependencies.R
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

cran_packages <- c("optparse", "readr", "httr")
for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

bioc_packages <- c("AnnotationForge", "biomaRt")
for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
}
EOF

echo "📦 Installing R packages..."
Rscript install_R_dependencies.R

echo "🧹 Cleaning up temporary files..."
rm install_R_dependencies.R

echo "✅ Done! You can now safely run create_annotation_package.R."
