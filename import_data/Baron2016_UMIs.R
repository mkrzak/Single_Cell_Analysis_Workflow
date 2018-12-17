main.path="/home/monika/Desktop/Single_Cell_Analysis_Workflow"
setwd(main.path)

library(SingleCellExperiment)

# human1
h1 <- read.csv("import_data/GSM2230757_human1_umifm_counts.csv", header = T)
rownames(h1) <- h1[,1]
cell_type <- as.character(h1$assigned_cluster)
h1 <- h1[,4:ncol(h1)]
h1 <- t(h1)

sce <- SingleCellExperiment(list(counts=h1), colData=data.frame(cell_type))
save(sce, file ="data/Baron2016_UMIs.RData")
