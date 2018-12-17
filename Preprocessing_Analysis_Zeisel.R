#Zeisel2015 dataset:
#Cells from Mouse Cortex and Hippocampus tissues 
#Link to study: http://science.sciencemag.org/content/347/6226/1138
#Data type: UMI counts
#Spike-ins: yes

main.path="/home/monika/Desktop/Single_Cell_Analysis_Workflow"
setwd(main.path)

#########################     SCATER     ####################################
library(SingleCellExperiment)
library(scater)
library(mclust)

#LOAD DATA IN SINGLE CELL EXPERIENT OBJECT AND ACCESS BASIC INFORMATION
#############################################################################
load("data/Zeisel2015_UMIs.RData")

sce
counts(sce)[1:4,1:4]
#Access colData information
names(colData(sce))
table(colData(sce)$tissue)
table(colData(sce)$level1class)

#QUALITY CONTROL and DATA PREPROCESSING
#############################################################################

#Calculate Quality Metrics
sce <- calculateQCMetrics(sce)

#Check quality metrics
table(rowData(sce)$is_feature_control_Spike) #is gene a spike-in
rownames(sce[grep("ERCC", rownames(sce)),])

#Quality metrics per gene
rowData(sce)$mean_counts[1:5] #average count per gene
rowData(sce)$n_cells_by_counts[1:5] #number of cells with non zero counts per gene
rowData(sce)$pct_dropout_by_counts[1:5] #% of dropouts per gene

#Quality metrics per cell
colData(sce)$total_counts[1:5] #total count for each cell
colData(sce)$total_features_by_counts[1:5] #number of genes with non zero counts per cell
colData(sce)$pct_counts_Spike[1:5] #% of spike-in counts per cell

#Plot histograms based on quality metrics
hist(rowData(sce)$pct_dropout_by_counts, xlab="Percentage of dropouts", ylab="Number of genes", breaks=50)

hist(sce$total_counts,breaks=50,  xlab="Library size", ylab="Number of cells")
hist(sce$total_features_by_counts, xlab="Number of expressed genes", breaks=50,  ylab="Number of cells")
hist(sce$pct_counts_Spike, xlab="ERCC proportion (%)", ylab="Number of cells", breaks=50)

#Automaticly filter Outlier cells based on quality metrics
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", log=TRUE)
spike.drop <- isOutlier(sce$pct_counts_Spike, nmads=3, type="higher")
sce <- sce[,!(libsize.drop | feature.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
    BySpike=sum(spike.drop), Remaining=ncol(sce))

#Filter out spike-in Genes
sce=sce[-which(grepl("^ERCC-", rownames(sce))),]
sce

#Filter genes by the average count per gene
rowData(sce)$ave.counts <- calcAverage(sce, exprs_values = "counts", use_size_factors=FALSE)
to.keep <- rowData(sce)$ave.counts > 0
sce <- sce[to.keep,]
summary(to.keep)
dim(sce)

#Normalize counts to remove library size bias
cpm(sce) <- calculateCPM(sce, use_size_factors=FALSE)
logcounts(sce) <- log2(cpm(sce)+1)

sce
cpm(sce)[1:4,1:4]

#DATA EXPLORATION WITH DIMENSIONALITY REDUCTION
#############################################################################
set.seed(1000)
sce <- runPCA(sce)
sce <- runTSNE(sce)

plotPCA(sce, colour_by="level1class", run_args=c(exprs_values = "logcounts"))
plotTSNE(sce, colour_by="level1class", run_args=c(exprs_values = "logcounts"))

set.seed(1000)
tsne0 <- plotTSNE(sce, colour_by="level1class")
tsne1 <- plotTSNE(sce, colour_by="Neurod6")
tsne2 <- plotTSNE(sce, colour_by="Mog") 
multiplot(tsne0, tsne1, tsne2, cols=3)

#ADDITIONAL PLOTS
#############################################################################
plotHighestExprs(sce, exprs_values = "counts")
plotExpression(sce, c("Malat1", rownames(sce)[1]), x = "level1class", colour_by = "level1class", exprs_values = "logcounts")
plotExpression(sce, c("Neurod6", "Mog", "Gad1"), x = "level1class", colour_by = "level1class", exprs_values = "logcounts")

save(sce, file="data/Zeisel2015_UMIs_QC.RData")

#########################     Monocle    ####################################
library(monocle)

#Create CellDataSet
#############################################################################
rowData(sce)$gene_short_name <- rownames(sce)

pd <- new("AnnotatedDataFrame", data = data.frame(colData(sce))) #cell info
fd <- new("AnnotatedDataFrame", data = data.frame(rowData(sce))) #gene info
cds <- newCellDataSet(cellData = counts(sce), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())  
cds <- estimateSizeFactors(cds)
cds

#Cluster cells to find cell populations
#############################################################################
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'tSNE', verbose = T) #reduce dimension with tSNE
cds <- clusterCells(cds, num_clusters = NULL) #Cluster cells

plot_cell_clusters(cds, color_by = 'Cluster') #check clusterization

#Compare with annotated cell types
p1 <- plot_cell_clusters(cds, color_by = 'level1class')
p2 <- plot_cell_clusters(cds, color_by = 'Cluster')
multiplot(p1, p2, cols=2)

adjustedRandIndex(cds$level1class, cds$Cluster)

#Select highly variable genes and check if you improved clustering
#############################################################################
cds <- estimateDispersions(cds)
disp <- dispersionTable(cds)
head(disp)
hvg <- subset(disp, mean_expression >= 0.5 & dispersion_empirical >= 0.1)
cds <- setOrderingFilter(cds, ordering_genes = hvg$gene_id)
plot_ordering_genes(cds)

cds_subset <- cds[hvg$gene_id,]
cds_subset <- reduceDimension(cds_subset, max_components = 2, reduction_method = 'tSNE', verbose = T)
cds_subset <- clusterCells(cds_subset, num_clusters = NULL)

#Compare with annotated cell types
p1 <- plot_cell_clusters(cds_subset, color_by = 'level1class')
p2 <- plot_cell_clusters(cds_subset, color_by = 'Cluster')
multiplot(p1, p2, cols=2)

adjustedRandIndex(cds_subset$level1class, cds_subset$Cluster)



