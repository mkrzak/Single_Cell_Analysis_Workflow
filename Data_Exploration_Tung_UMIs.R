#Tung2017 dataset:
#Three lines of human induced pluripotent stem cell (iPSC)
#Link to study: https://www.nature.com/articles/srep39921
#Data type: UMI counts
#Spike-ins: yes

main.path="/home/monika/Desktop/Single_Cell_Analysis_Workflow"
setwd(main.path)

#########################     SCATER     ####################################
library(SingleCellExperiment)
library(scater)

#CREATE SINGLE CELL EXPERIENT OBJECT 
#############################################################################
counts <- read.table("data/Tung/molecules.txt", sep = "\t")
anno <- read.table("data/Tung/annotation.txt", sep = "\t", header = TRUE)

head(counts[, 1:3])
head(anno)

#Create SingleCellExperiment object
umi <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData = anno)

table(colData(umi)$individual)
table(colData(umi)$replicate)
table(colData(umi)$batch)

#Filter not expressed genes
keep_feature <- rowSums(counts(umi) > 0) > 0
umi <- umi[keep_feature, ]

rownames(umi[grep("ERCC", rownames(umi)),])
isSpike(umi, "Spike") <- grepl("^ERCC-", rownames(umi))

umi <- calculateQCMetrics(umi)

names(rowData(umi))

#Quality metrics per gene
rowData(umi)$mean_counts[1:5] #average count per gene
rowData(umi)$n_cells_by_counts[1:5] #number of cells with non zero counts per gene
rowData(umi)$pct_dropout_by_counts[1:5] #% of dropouts per gene

#Quality metrics per cell
colData(umi)$total_counts[1:5] #total count for each cell
colData(umi)$total_features_by_counts[1:5] #number of genes with non zero counts per cell
colData(umi)$pct_counts_Spike[1:5] #% of spike-in counts per cell

#Filter by total counts
hist(umi$total_counts,breaks = 100)
abline(v = 25000, col = "red")
filter_by_total_counts <- (umi$total_counts > 25000)
table(filter_by_total_counts)

#Filter by number of expressed features
hist(umi$total_features_by_counts,breaks = 100)
abline(v = 7000, col = "red")
filter_by_expr_features <- (umi$total_features_by_counts > 7000)
table(filter_by_expr_features)

umi$use <- (filter_by_expr_features & filter_by_total_counts)

#Check filtered cells in PCA
logcounts(umi) <- log2(counts(umi)+1)
scater::plotPCA(umi, colour_by = "use", shape_by="replicate")
scater::plotPCA(umi, colour_by = "use", shape_by="individual")
scater::plotPCA(umi, shape_by = "use", colour_by="batch")

#Filter lowly expressed genes (at least two cells should contain more than 1 transcript)
filter_genes <- apply(counts(umi[ , colData(umi)$use]), 1, function(x) length(x[x > 1]) >= 2)
table(filter_genes)
rowData(umi)$use <- filter_genes

#Apply filtering
umi[rowData(umi)$use, colData(umi)$use]

#Remove Spike-ins
endog_genes <- !rowData(umi)$is_feature_control
umi <- umi[endog_genes, ]

#Search for possible problems in the data
logcounts(umi) <- log2(counts(umi)+1)
plotPCA(umi, size_by = "total_counts", colour_by="individual") #library size explains highest variability in the data
plotExplanatoryVariables(umi,  variables = c("total_counts","individual", "batch")) #strong batch effect

#Normalize 
cpm(umi) <- calculateCPM(umi, use_size_factors=FALSE)
logcounts(umi) <- log2(cpm(umi)+1)
plotPCA(umi, size_by = "total_counts", colour_by="individual") #library size explains highest variability in the data
plotExplanatoryVariables(umi,  variables = c("total_counts","individual", "batch")) #strong batch effect

#Explotre possible outliers
scater::plotPCA(umi, colour_by="individual", run_args=c(ntop = 50))
scater::plotPCA(umi, colour_by="batch",run_args=c(ntop = 50))














