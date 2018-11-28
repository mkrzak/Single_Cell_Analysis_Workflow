# Run below commands to install required packages

#Compulsory for the analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment", version = "3.8")
BiocManager::install("scater", version = "3.8")

#Optional for Import data
BiocManager::install("BiocFileCache", version = "3.8")
BiocManager::install("org.Mm.eg.db", version = "3.8")