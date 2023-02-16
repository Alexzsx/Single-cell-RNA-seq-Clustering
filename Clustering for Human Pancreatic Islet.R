install.packages("devtools")
devtools::install_github("VCCRI/CIDR")

library(cidr)

## Read in data 

pancreaticIsletTags <- read.csv("CIDR_data/pancreaticIsletTags.csv")
rownames(pancreaticIsletTags) <- pancreaticIsletTags[,1]
pancreaticIsletTags <- pancreaticIsletTags[,-1]
info <- read.csv("CIDR_data/SraRunTable2.txt",sep="\t")
cellType <- info$assigned_cell_type_s[match(colnames(pancreaticIsletTags),info$Sample_Name_s)]
cellType <- factor(cellType)
types <- levels(cellType)

## Assign each cell type a color
scols <- c("red","green","blue","grey","purple","brown")
cols <- rep(NA,length(cellType))
true_label <- rep(NA,length(cellType))
for (i in 1:length(cols)){
  cols[i] <- scols[which(types==cellType[i])]
  true_label[i] <- which(types==cellType[i])
}



# Standard principal component analysis using prcomp
priorTPM <- 1
pan10 <- pancreaticIsletTags[rowSums(pancreaticIsletTags)>10,]
pan10_lcpm <- log2(t(t(pan10)/colSums(pan10))*1000000+priorTPM)
PC_lcpm <- prcomp(t(pan10_lcpm))
plot(PC_lcpm$x[,c(1,2)], col=cols, pch=10, 
     xlab="PC1", ylab="PC2", 
     main="Principal Component Analysis (prcomp)")
legend("topright", legend=types, col=scols, pch=10)

#### CIDR ######
################

start_time <- Sys.time()
scPan <- scDataConstructor(as.matrix(pancreaticIsletTags))
scPan <- determineDropoutCandidates(scPan)
scPan <- wThreshold(scPan)
scPan <- scDissim(scPan)
scPan <- scPCA(scPan)
scPan <- nPC(scPan)
scPan <- scCluster(scPan)

end_time <- Sys.time()
t_CIDR <- end_time - start_time
print(t_CIDR)
levels(scPan@clusters)
## Use Adjusted Rand Index to measure the accuracy of CIDR clustering
ARI_CIDR <- adjustedRandIndex(scPan@clusters,cols)
ARI_CIDR

pdf("els-cas-templates/els-cas-templates/figs/CIDR_2.pdf", width = 6, height = 6) 
par(xpd=TRUE)
plot(scPan@PC[,c(1,2)], col=c(scols)[true_label], 
     main="CIDR", xlab="PC1", ylab="PC2",pch=20)
legend("bottomleft", legend=c(1:4), col=scols, pch=20, bty = "n")
dev.off()

## 3D
BiocManager::install("plot3D") 
library(plot3D)
scatter3D(scPan@PC[,1], scPan@PC[,2], scPan@PC[,3],
          xlab="PC1", ylab="PC2", zlab="PC3",
          colvar=NULL, col=cols, pch=scPan@clusters, 
          phi=5, theta=55, pch=20)
legend("bottomright", legend=types, col=scols, pch=15)

## 3D - interactive
BiocManager::install("rgl") 
library(rgl)
plot3d(scPan@PC[,1:3], col=cols, 
       xlab="PC1", ylab="PC2", zlab="PC3")




#### SC3 #######
################
library(SingleCellExperiment)
library(SC3)
start_time <- Sys.time()

# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(pancreaticIsletTags),
    logcounts = log2(as.matrix(pancreaticIsletTags) + 1)
  ), 
  colData = cellType
)

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
#sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
# define spike-ins
isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)

# run SC3 with ks clusters
# sce <- sc3(sce, ks = 10, biology = FALSE, svm_num_cells = 5000)
# sce_svm <- sc3_run_svm(sce, ks = 10)
sce <- sc3(sce, ks = 6, biology = FALSE, rand_seed = sample(1:10000,1))

# cell results
col_data <- colData(sce)
# col_data <- colData(sce_svm)

end_time <- Sys.time()
t_SC3 <- end_time - start_time
print(t_SC3)


ARI_SC3 <- adjustedRandIndex(col_data$sc3_6_clusters,cols)
ARI_SC3

sc3_plot_expression(sce, 6, show_pdata = NULL)




#### Seurat ####
################
BiocManager::install("Seurat") 

library(Seurat)
library(dplyr)
library(Matrix)
# Load the PBMC dataset

#pbmc.data <- read.csv("data_sc3_seurat/A_100.csv", row.names = 1);

start_time <- Sys.time()

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pan10))
sparse.size <- object.size(x = pan10)

mincell=0 # Keep all genes expressed in >= mincell cells (parameter)
mingene=0 # Keep all cells with at least mingene detected genes (parameter)

pbmc <- CreateSeuratObject(pan10)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc, npcs=10, ndims.print = 1:5, nfeatures.print = 5)
j=1.75; # a tunable parameter
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:10, nn.eps = 0.5)
results <- FindClusters(object = pbmc, resolution = j)

end_time <- Sys.time()
t_Seurat = end_time - start_time
print(t_Seurat)

ARI_Seurat <- adjustedRandIndex(results@active.ident,cols)
ARI_Seurat



#### RaceID ####
################
BiocManager::install("RaceID") 

library(RaceID)

start_time <- Sys.time()
sc <- SCseq(pancreaticIsletTags)
sc <- filterdata(sc, mintotal=3000, minexpr=5, minnumber=1)
# k-means clustering
sc <- compdist(sc,metric="pearson")
sc <- clustexp(sc, clustnr=20,bootnr=50, cln=0,rseed=17000)
# compute t-SNE map
sc <- comptsne(sc,rseed=15555,perplexity=10)
# detect outliers and redefine clusters
sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,outdistquant=.75)
end_time <- Sys.time()
t_RaceID = end_time - start_time
print(t_RaceID)
## Use Adjusted Rand Index to measure the accuracy of CIDR clustering
ARI_RaceID <- adjustedRandIndex(sc@cluster$kpart,cols)
ARI_RaceID

## diagnostic plots
# gap statistics
plotgap(sc)
# silhouette of k-means clusters
plotsilhouette(sc)
# Jaccard's similarity of k-means clusters
plotjaccard(sc)
# barchart of outlier probabilities
plotoutlierprobs(sc)
# regression of background model
plotbackground(sc)
# dependence of outlier number on probability threshold (probthr)
plotsensitivity(sc)
# heatmap of k-means cluster
clustheatmap(sc,final=FALSE,hmethod="single")
# heatmap of final cluster
clustheatmap(sc,final=TRUE,hmethod="single")

sc@cluster$kpart
sc@tsne
sc <- comptsne(sc)





pdf("els-cas-templates/els-cas-templates/figs/RaceID_2.pdf", width = 6, height = 6) 
plot(sc@tsne,col=c(scols)[true_label],xlab="t-SNE component 1", 
     ylab="t-SNE component 2",pch=20,main="RaceID3")
legend("bottomleft", legend=types, col=scols, pch=20, bty = "n")
dev.off()

#### SIMLR ####
###############
BiocManager::install("SIMLR") 

library(SIMLR)

start_time <- Sys.time()

res_example1 = SIMLR(X=pancreaticIsletTags,c=6, normalize = TRUE)
pdf("els-cas-templates/els-cas-templates/figs/SIMLR_2.pdf", width = 6, height = 6) 
par(xpd=TRUE)
plot(res_example1$ydata,col=c(scols)[true_label],xlab="SIMLR component 1", ylab="SIMLR component 2",pch=20,main="SIMILR")
legend(-30, 55, legend=types, col=scols, pch=20, bty = "n")
dev.off()
cluster = res_example1$y

ARI_SIMLR <- adjustedRandIndex(cluster$cluster,cols)
ARI_SIMLR

end_time <- Sys.time()
t_SIMLR = end_time - start_time
print(t_SIMLR)

#### monocle ####
###############
BiocManager::install("monocle") 
library(monocle)


start_time <- Sys.time()
pan10 <- pancreaticIsletTags[rowSums(pancreaticIsletTags)>10,]
dataset <- as.matrix(pan10)

data <- newCellDataSet(dataset)
data <- estimateSizeFactors(data)
data <- reduceDimension(data, max_components = 3, num_dim = 10, reduction_method = 'tSNE', verbose = T, perplexity =5)
data <- clusterCells(data, num_clusters = 6)

end_time <- Sys.time()
t_monocle = end_time - start_time
print(t_monocle)


ARI_monocle <- adjustedRandIndex(data$Cluster,cols)
ARI_monocle

pdf("els-cas-templates/els-cas-templates/figs/Monocle2_2.pdf", width = 3.5, height = 3.5) 
plot_cell_clusters(data, x = 1, y = 2, color_by = "Cluster",
                   markers = NULL, cell_size = 1.5, show_cell_names = true_label) + 
  ggtitle("Monocle2") + theme(legend.position = "right", title=element_text(size=8), 
                              legend.text = element_text(size = 5), axis.text.x = element_text(size = 6),
                              axis.text.y = element_text(size = 6))
dev.off()



library(ggplot2)


Clusterring_methods<-c("RaceID3", "Monocle2", "SIMLR", "Seurat", "SC3", "CIDR")
ARI<-c(ARI_RaceID, ARI_monocle, ARI_SIMLR, ARI_Seurat, ARI_SC3, ARI_CIDR)
Time <- c(t_RaceID, t_monocle, t_SIMLR, t_Seurat, t_SC3, t_CIDR)
Cdata <-data.frame(Clusterring_methods, ARI, Time)



p <-ggplot(Cdata, aes(Clusterring_methods, ARI))
p +geom_bar(stat = "identity",aes(fill = Clusterring_methods)) + 
  geom_text(aes(label=as.character(round(ARI,3))),vjust=-0.5) +
  xlab("Clusterring methods") + ylab("ARI") + theme_bw() + guides(fill=guide_legend(title="Methods"))+
  ggtitle("Human pancreatic islet (GSE73727)") + 
  theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1) + ylim(0, 1)
ggsave("els-cas-templates/els-cas-templates/figs/Figure3.pdf", width = 6, height = 6, dpi = 600)


p <-ggplot(Cdata, aes(Clusterring_methods, Time))
p +geom_bar(stat = "identity",aes(fill = Clusterring_methods)) + 
  geom_text(aes(label=as.character(round(Time,2))),vjust=-0.5) +
  xlab("Clusterring methods") + ylab("Running Time (s)") + theme_bw() + guides(fill=guide_legend(title="Methods"))+
  ggtitle("Human pancreatic islet (GSE73727)") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
ggsave("els-cas-templates/els-cas-templates/figs/Figure4.pdf", width = 6, height = 6, dpi = 600)


name = "CIDR_data/GSE73727.csv"
write.csv(Cdata, name)
