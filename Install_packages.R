# Run below commands to install required packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment", version = "3.8")
BiocManager::install("scater", version = "3.8")
BiocManager::install("monocle", version = "3.8")
install.packages("mclust")
