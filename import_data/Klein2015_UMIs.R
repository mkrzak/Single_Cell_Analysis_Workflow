main.path="/home/monika/Desktop/Single_Cell_Analysis_Workflow"
setwd(main.path)

library(SingleCellExperiment)

d0 <- read.csv("import_data/klein/GSM1599494_ES_d0_main.csv", header = FALSE)
d2 <- read.csv("import_data/klein/GSM1599497_ES_d2_LIFminus.csv", header = FALSE)
d4 <- read.csv("import_data/klein/GSM1599498_ES_d4_LIFminus.csv", header = FALSE)
d7 <- read.csv("import_data/klein/GSM1599499_ES_d7_LIFminus.csv", header = FALSE)
d <- cbind(d0, d2[,2:ncol(d2)], d4[,2:ncol(d4)], d7[,2:ncol(d7)])
rownames(d) <- d[,1]
d <- d[,2:ncol(d)]
colnames(d) <- paste0("cell", 1:ncol(d))

### ANNOTATIONS
ann <- data.frame(
  day = c(rep("d0", ncol(d0) - 1), 
                 rep("d2", ncol(d2) - 1),
                 rep("d4", ncol(d4) - 1),
                 rep("d7", ncol(d7) - 1)))
rownames(ann) <- colnames(d)

sce <- SingleCellExperiment(assays = list(counts = as.matrix(d)), colData=ann)
save(sce, file ="data/Klein2015_UMIs.RData")
