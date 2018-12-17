#Klein2015 dataset:
#Link to study: https://www.cell.com/cell/fulltext/S0092-8674(15)00500-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867415005000%3Fshowall%3Dtrue
#Cells: Mouse Embro Stem cells sampled after days: 0, 2, 4, 7
#Data type: UMI counts
#Spike-ins: no

main.path="/home/monika/Desktop/Single_Cell_Analysis_Workflow"
setwd(main.path)

#########################     SCATER     ####################################
library(SingleCellExperiment)
library(scater)
library(mclust)

#LOAD DATA IN SINGLE CELL EXPERIENT OBJECT AND ACCESS BASIC INFORMATION
#############################################################################
load("data/Klein2015_UMIs.RData")

sce
counts(sce)[1:4,1:4]
#Access colData information
table(colData(sce)$day)

#QUALITY CONTROL and DATA PREPROCESSING
#############################################################################

#Calculate Quality Metrics
sce <- calculateQCMetrics(sce)

#Check quality metrics
rowData(sce)$mean_counts[1:5] #average count per gene
rowData(sce)$n_cells_by_counts[1:5] #number of cells with non zero counts per gene
rowData(sce)$pct_dropout_by_counts[1:5] #% of dropouts per gene

colData(sce)$total_counts[1:5] #total count for each cell
colData(sce)$total_features_by_counts[1:5] #number of genes with non zero counts per cell

#Plot histograms based on quality metrics
hist(sce$total_counts,breaks=50,  xlab="Library size", ylab="Number of cells")
hist(sce$total_features_by_counts, xlab="Number of expressed genes", breaks=50,  ylab="Number of cells")

#Filter Outlier cells based on quality metrics
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", log=TRUE)

sce <- sce[,!(libsize.drop | feature.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), Remaining=ncol(sce))

#Filter genes by the average count per gene
rowData(sce)$ave.counts <- calcAverage(sce, exprs_values = "counts", use_size_factors=FALSE)
to.keep <- rowData(sce)$ave.counts > 0
sce <- sce[to.keep,]
summary(to.keep)
dim(sce)

#DATA EXPLORATION WITH DIMENSIONALITY REDUCTION
#############################################################################
set.seed(1000)

logcounts(sce) <- log2(counts(sce)+1)

sce <- runPCA(sce)
sce <- runTSNE(sce)

#Explore PCA on non normalized counts
plotPCA(sce, colour_by="day", size_by="total_counts", run_args=c(exprs_values = "logcounts"))
plotExplanatoryVariables(sce,  variables = c("total_counts")) #strong batch effect

#Explore PCA after normalizing the counts
cpm(sce) <- calculateCPM(sce, use_size_factors=FALSE)
logcounts(sce) <- log2(cpm(sce)+1)
plotPCA(sce, colour_by="day", size_by="total_counts", run_args=c(exprs_values = "logcounts"))
plotExplanatoryVariables(sce,  variables = c("total_counts")) #strong batch effect


#ADDITIONAL PLOTS
#############################################################################
plot(sce$day, sce$total_counts, xlab="Group", ylab="Total sum of counts per cell")
save(sce, file="data/Klein2015_UMIs_QC.RData")

#########################     Monocle    ####################################
library(monocle)

#Create CellDataSet
#############################################################################
rowData(sce)$gene_short_name <- rownames(sce)

pd <- new("AnnotatedDataFrame", data = data.frame(colData(sce)))
fd <- new("AnnotatedDataFrame", data = data.frame(rowData(sce)))
cds <- newCellDataSet(cellData = counts(sce), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())  
cds <- estimateSizeFactors(cds)
cds


#Cluster cells 
#############################################################################

cds_subset <- reduceDimension(cds_subset, max_components = 2, reduction_method = 'tSNE', verbose = T)
cds_subset <- clusterCells(cds_subset, num_clusters = NULL)

p1 <- plot_cell_clusters(cds_subset, color_by = 'day')
p2 <- plot_cell_clusters(cds_subset, color_by = 'Cluster')
multiplot(p1, p2, cols=2)

adjustedRandIndex(cds$day, cds$Cluster)

#Select highly variable genes and check if it improves clustering
#############################################################################
cds <- estimateDispersions(cds)
disp <- dispersionTable(cds)
head(disp)
hvg <- subset(disp, mean_expression >= 0.5 & dispersion_empirical >= 0.05)
cds <- setOrderingFilter(cds, ordering_genes = hvg$gene_id)
plot_ordering_genes(cds)

cds_subset <- cds[hvg$gene_id,]
cds_subset <- reduceDimension(cds_subset, max_components = 2, reduction_method = 'tSNE', verbose = T)
cds_subset <- clusterCells(cds_subset, num_clusters = NULL)

p1 <- plot_cell_clusters(cds_subset, color_by = 'day')
p2 <- plot_cell_clusters(cds_subset, color_by = 'Cluster')
multiplot(p1, p2, cols=2)

adjustedRandIndex(cds_subset$day, cds_subset$Cluster)







