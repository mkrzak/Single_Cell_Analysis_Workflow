#Tung2017 dataset:
#Three lines of human induced pluripotent stem cell (iPSC)
#Link to study: https://www.nature.com/articles/srep39921
#Data type: Read counts
#Spike-ins: yes

main.path="/home/monika/Desktop/Single_Cell_Analysis_Workflow"
setwd(main.path)

#########################     SCATER     ####################################
library(SingleCellExperiment)
library(scater)

#CREATE SINGLE CELL EXPERIENT OBJECT 
#############################################################################
counts <- read.table("data/Tung/reads.txt", sep = "\t")
anno <- read.table("data/Tung/annotation.txt", sep = "\t", header = TRUE)

head(counts[ , 1:3])
head(anno)

#Create SingleCellExperiment object
reads <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData = anno)

#Filter not expressed genes
keep_feature <- rowSums(counts(reads) > 0) > 0
reads <- reads[keep_feature, ]

isSpike(reads, "Spike") <- grepl("^ERCC-", rownames(reads))

reads <- calculateQCMetrics(reads)

#Quality metrics per gene
rowData(reads)$mean_counts[1:5] #average count per gene (higher than in UMIs)
rowData(reads)$n_cells_by_counts[1:5] #number of cells with non zero counts per gene
rowData(reads)$pct_dropout_by_counts[1:5] #% of dropouts per gene

#Quality metrics per cell
colData(reads)$total_counts[1:5] #total count for each cell (higher than in UMIs)
colData(reads)$total_features_by_counts[1:5] #number of genes with non zero counts per cell
colData(reads)$pct_counts_Spike[1:5] #% of spike-in counts per cell

#Filter by total counts
hist(reads$total_counts,breaks = 100)
abline(v = 1.3e6, col = "red")
filter_by_total_counts <- (reads$total_counts > 1.3e6)
table(filter_by_total_counts)

#Filter by number of expressed features
hist(reads$total_features_by_counts ,breaks = 100)
abline(v = 7000, col = "red")
filter_by_expr_features <- (reads$total_features_by_counts > 7000)
table(filter_by_expr_features)

reads$use <- (filter_by_expr_features & filter_by_total_counts)

#Check filtered cells in PCA
logcounts(reads) <- log2(counts(reads)+1)
scater::plotPCA(reads, colour_by = "use", shape_by="replicate")
scater::plotPCA(reads, colour_by = "use", shape_by="individual")
scater::plotPCA(reads, shape_by = "use", colour_by="batch")

#Filter lowly expressed genes (at least five reads in at least two cells)
filter_genes <- apply(counts(reads[, colData(reads)$use]), 1, function(x) length(x[x > 5]) >= 2)
rowData(reads)$use <- filter_genes
table(filter_genes)

#Apply filtering
reads[rowData(reads)$use, colData(reads)$use]

#Normalize 
cpm(reads) <- calculateCPM(reads, use_size_factors=FALSE)
logcounts(reads) <- log2(cpm(reads)+1)

#Explotre possible outliers
scater::plotPCA(reads, colour_by="individual", run_args=c(ntop = 50))
scater::plotPCA(reads, colour_by="batch",run_args=c(ntop = 50))























