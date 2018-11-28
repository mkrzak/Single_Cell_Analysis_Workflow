setwd("/media/monika/683375f9-7173-4a68-ba60-f33b0c47a470/Analysis/Single_Cell_Analysis_Workflow")

# library(SingleCellExperiment)
# library(scater)
# library(scran)

library(simpleSingleCell)

#LOAD DATA IN SINGLE CELL EXPERIENT OBJECT AND ACCESS BASIC INFORMATION
#############################################################################
load("data/Zeisel2015_UMIs.RData")
sce
counts(sce)[1:4,1:4]
#Access colData information
table(colData(sce)$tissue)
table(colData(sce)$level1class)

#QUALITY CONTROL
#############################################################################
#Calculate Quality Metrics with scater
sce <- calculateQCMetrics(sce) 

#Check quality metrics
# table(rowData(sce)$is_feature_control_Spike) #is gene a spike-in
# rownames(sce[grep("ERCC", rownames(sce)),])
# rowData(sce)$mean_counts[1:5] #average count per gene
# rowData(sce)$n_cells_by_counts[1:5] #number of cells with non zero counts per gene
# rowData(sce)$pct_dropout_by_counts[1:5] #% of dropouts per gene
# rowData(sce)$total_counts[1:5] #library size
# colData(sce)$total_features_by_counts[1:5] #number of genes with non zero counts per cell
# colData(sce)$pct_counts_Spike[1:5] #% of spike-in counts per cell

#Plot histograms based on quality metrics
hist(sce$total_counts,breaks=20, ylab="Number of cells",  xlab="Library size")
hist(sce$total_features_by_counts, xlab="Number of expressed genes", breaks=20,  ylab="Number of cells")
hist(sce$pct_counts_Spike, xlab="ERCC proportion (%)", ylab="Number of cells", breaks=20)

#Filter Outlier cells based on quality metrics
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", log=TRUE)
spike.drop <- isOutlier(sce$pct_counts_Spike, nmads=3, type="higher")
sce <- sce[,!(libsize.drop | feature.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
    BySpike=sum(spike.drop), Remaining=ncol(sce))
dim(sce)

#Filter genes by the average count per gene
rowData(sce)$ave.counts <- calcAverage(sce, exprs_values = "counts", use_size_factors=FALSE)
hist(log10(rowData(sce)$ave.counts), breaks=100, main="", col="grey",
    xlab=expression(Log[10]~"average count"))
to.keep <- rowData(sce)$ave.counts > 0
sce <- sce[to.keep,]
summary(to.keep)
dim(sce)

#Normalize data
# based on pulled size factors
sce <- computeSpikeFactors(sce, type="Spike", general.use=FALSE)

set.seed(1000)
colData(sce)$clusters <- quickCluster(sce, min.mean=0.1, method="igraph")
sce <- computeSumFactors(sce, cluster=colData(sce)$clusters, min.mean=0.1)
head(sizeFactors(sce))
sce <- normalize(sce, exprs_values = "counts", return_log = TRUE)
sce
logcounts(sce)[1:4,1:4]

#OR using CPM
# logcounts(sce) <- calculateCPM(sce)
# logcounts(sce)[1:4,1:4]

#For read counts
#calculateFPKM()
#calculateTPM()

#DATA EXPLORATION WITH DIMENSIONALITY REDUCTION
#############################################################################

plotPCA(sce, colour_by="level1class")
plotTSNE(sce, colour_by="level1class")

set.seed(1000)
#sce <- runPCA(sce)
#sce <- runTSNE(sce)
#sce <- runMDS(sce)

set.seed(1000)
tsne1 <- plotTSNE(sce, colour_by="Neurod6") 
tsne2 <- plotTSNE(sce, colour_by="Mog") 
multiplot(tsne1, tsne2, cols=2)

pca1 <- plotReducedDim(sce, use_dimred="PCA", colour_by="Neurod6") 
pca2 <- plotReducedDim(sce, use_dimred="PCA", colour_by="Mog")
multiplot(pca1, pca2, cols=2)

plotReducedDim(sce, use_dimred="PCA", ncomponents=3, colour_by ="level1class")


#CLUSTERING CELLS INTO PUTATIVE SUBPOPULATIOS
#############################################################################
#Build graph
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
#Find communities in the graph - clusters
cluster.out <- igraph::cluster_walktrap(snn.gr)
table(cluster.out$membership)
sce$cluster <- factor(cluster.out$membership)
plotTSNE(sce, colour_by="cluster") 

#DETECT MARKER GENES BETWEEN SUBPOPULATIONS
#############################################################################
markers <- findMarkers(sce, sce$cluster, direction="up")[["1"]]
head(markers[,1:6], 3)
top.markers <- rownames(markers)[markers$Top <= 10]

plotHeatmap(sce, features=top.markers, columns=order(cluster.out$membership),
    colour_columns_by="cluster", cluster_cols=FALSE, 
    center=TRUE, symmetric=TRUE, zlim=c(-5, 5))

saveRDS(file="brain_data.rds", sce)


sessionInfo()






