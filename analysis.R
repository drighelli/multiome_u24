# Downloading dataset
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6543/GSM6543828/suppl/GSM6543828_HBr_15_A.tar.gz
# 
library(TENxMultiomeTools)
# dir.create("outs")
# sce <- read10xMultiome("HBr_15_A/")

library(Matrix)
mat <- readMM("HBr_15_A/outs/filtered_feature_bc_matrix/matrix.mtx")
dim(mat)
bc <- read.table("HBr_15_A/outs/filtered_feature_bc_matrix/barcodes.tsv", sep="\t")
feat <- read.table("HBr_15_A/outs/filtered_feature_bc_matrix/features.tsv", sep="\t")
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays=list("counts"=mat), rowData=feat)
colnames(sce ) <- bc$V1
rownames(sce) <-  feat$V1
sce <- splitAltExps(sce, feat$V3, ref="Gene Expression")

# library("EnsDb.Hsapiens.v86")
# ensdb <- EnsDb.Hsapiens.v86
# seqlevelsStyle(ensdb) <- "UCSC"
# library(Signac)
# library(Seurat)
# annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
# genome(annotation) <- "hg38"
sce <- scuttle::addPerCellQC(sce)
library(scuttle)
colData(sce)
df <- perCellQCFilters(sce, sub.fields=c("altexps_Peaks_detected", 
                            "altexps_Peaks_sum", "altexps_Peaks_percent"))
summary(df$discard)
sce <- sce[,!df$discard]
sce <- logNormCounts(sce)

sce <- computeDoublets(sce, method="griffiths")
colData(sce)
sce <- sce[,colData(sce)$dbl.calls=="singlet"]
library(scater)
sce <- runPCA(sce)
sce <- runTSNE(sce)
plotTSNE(sce)
library(SingleR)

# # ref=celldex::HumanPrimaryCellAtlasData()
# refmd <- read.csv("mereu_cancer_ref/TICAtlas_metadata.csv")
# ref <- read.csv("mereu_cancer_ref/TIC")
# sceref <- SingleCellExperiment(assays=list(counts=ref), colData=refmd)
# colnames(sceref) <- rownames(refmd)
# # rownames(sceref) <- TBD
# scuttle::aggregateAcrossCells(sceref, ids="")
# colData(ref)
# rownames(sce) <- rowData(sce)$V2
# dflabs <- SingleR(test=sce, ref=ref, labels=ref$label.main)
# 
# table(dflabs$pruned.labels)
# colData(sce)$labels <- dflabs$pruned.labels
# colors <- palette.colors(n=17,palette = "Okabe-Ito") #colorblinded palette
# library(Polychrome)
# pal <- createPalette(17, seedcolors=c("#ff0000", "#00ff00", "#0000ff"))
# pal<- Polychrome::glasbey.colors(n=18)
# pal<-pal[-1]
# names(pal) <- NULL
# sce2 <- sce
# nms<-names(table(sce2$labels)[which(table(sce2$labels)>50)])
# sce2 <- sce[,sce$labels%in%nms]
# plotTSNE(sce2, colour_by="labels",) +
#     scale_color_manual(values=pal)
#     # scale_color_viridis_d(option="magma")


# wget https://singlecell.broadinstitute.org/single_cell/data/public/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers?filename=Whole_miniatlas_meta.csv
atl <- read.csv("breast_cancer_miniatlas/Whole_miniatlas_meta.csv")
atl <- atl[-1,]
table(atl$celltype_major)
library(Matrix)
# wget https://singlecell.broadinstitute.org/single_cell/data/public/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers?filename=matrix.mtx.gz
mat <- readMM("breast_cancer_miniatlas/matrix.mtx.gz")
# wget https://singlecell.broadinstitute.org/single_cell/data/public/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers?filename=barcodes.tsv.gz
bc <- read.table("breast_cancer_miniatlas/barcodes.tsv.gz", sep="\t")
# wget https://singlecell.broadinstitute.org/single_cell/data/public/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers?filename=features.tsv.gz
feat <- read.table("breast_cancer_miniatlas/features.tsv.gz", sep="\t")

scebcref <- SingleCellExperiment(assays=list("counts"=mat), rowData=feat, colData=atl)
rownames(scebcref) <- feat$V1
colnames(scebcref) <- atl$NAME
scebcref$ct_pat <- paste0(scebcref$celltype_major, "_", scebcref$Patient)

# scerefag <- scuttle::aggregateAcrossCells(scebcref, ids=scebcref$ct_pat)

library(muscat)
library(BiocParallel)
scebcrefagg <- aggregateData(scebcref, assay="counts", fun="sum", by=c("celltype_major", "Patient"), BPPARAM=MulticoreParam(workers=8))

agg_mat <- lapply(seq_along(assayNames(scebcrefagg)), function(i) assay(scebcrefagg, i)) |> do.call(cbind, args=_)
scebcreffin <- SingleCellExperiment(assays = list(counts = agg_mat), 
                                    colData = data.frame(celltype=rep(assayNames(scebcrefagg), each=ncol(scebcrefagg))))
colnames(scebcreffin) <- paste0(colnames(scebcreffin), "_", rep(assayNames(scebcrefagg), each=ncol(scebcrefagg)))
scebcreffin <- scebcreffin[, colSums(assay(scebcreffin))>0]
# cbind(colnames(scebcreffin), scebcreffin$celltype)

library(edgeR)
filter <- filterByExpr(assay(scebcreffin), min.count = 1, group = scebcreffin$celltype)
scebcreffin <- scebcreffin[filter,]
sizeFactors(scebcreffin) <- calcNormFactors(assay(scebcreffin))
scebcreffin <- logNormCounts(scebcreffin)
scebcreffin <- runPCA(scebcreffin)
plotPCA(scebcreffin, colour_by = "celltype")

library(BiocParallel)
library(scran)
rownames(sce) <- rowData(sce)$V2
filtered <- sce[intersect(rownames(sce), rownames(scebcreffin)),]
cl <- quickCluster(filtered, assay.type = "logcounts")
table(cl)
filtered <- computeSumFactors(filtered, clusters = cl, BPPARAM = MulticoreParam(workers=8))
filtered <- logNormCounts(filtered)
filtered <- runPCA(filtered)
filtered <- runTSNE(filtered)

nn.clusters <- clusterCells(filtered, use.dimred="PCA")
filtered$clusters <- nn.clusters
plotTSNE(filtered, colour_by = "clusters")

dflabs <- SingleR(test=filtered, ref=scebcreffin, labels=scebcreffin$celltype, clusters = filtered$clusters,
                  BPPARAM = MulticoreParam(workers=8))
dflabs

# clusters 2, 5, 6, 10 --> CAFs
# clusters 1, 7 --> cancer epithelial

filtered$celltype <- filtered$clusters
levels(filtered$celltype) <- dflabs$pruned.labels

plotTSNE(filtered, colour_by = "clusters", text_by = "celltype")

gr <- as.data.frame(matrix(unlist(strsplit(rowData(altExp(filtered))$V1, split = ":|-")), byrow = TRUE,ncol=3))
colnames(gr) <- c("chr", "start", "end")
gr <- makeGRangesFromDataFrame(gr)
names(gr) <- rowData(altExp(filtered))$V1
# library(ChIPseeker)
# BiocManager::install("org.Hs.eg.db")

rranno <- annotatePeak(gr, tssRegion=c(-3000, 3000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")@anno
rranno[which(rranno$SYMBOL == "COL1A1"),]

tapply(counts(altExp(filtered)[names(rranno[which(rranno$SYMBOL == "COL1A1"),]),]), filtered$clusters, sum)
filteredpb <- aggregateAcrossCells(altExp(filtered), ids=filtered$clusters)
counts(filteredpb)[names(rranno[which(rranno$SYMBOL == "COL1A1"),]),]


rranno[which(rranno$SYMBOL == "ETS1"),][5]
counts(filteredpb)[names(rranno[which(rranno$SYMBOL == "ETS1"),][5]),]
# library("Signac")
# sig <- CreateChromatinAssay(counts = counts(altExp(filtered)), sep=c(":", "-"))
# library("Seurat")
# # seu <- as.Seurat(filtered)
# seu <- CreateSeuratObject(counts=sig, assay="ATAC")
# seu <- RunTFIDF(seu)
# seu <- FindTopFeatures(seu, min.cutoff = "q0")
# seu <- RunSVD(seu)
# seu <- FindNeighbors(seu, reduction = "lsi", dims = 2:30)
# seu <- FindClusters(seu, algorithm = 3)
# 
# seu@meta.data$seurat_clusters <- filtered$clusters
#     
# seu@active.ident <- filtered$clusters
# seu <- RunUMAP(seu, reduction = "lsi", dims = 2:30)    
# DimPlot(seu, label=TRUE)
# ga <- GeneActivity(seu)
# granges(seu)
# library("BSgenome.Hsapiens.UCSC.hg38")
# seu <- RegionStats(seu, genome=BSgenome.Hsapiens.UCSC.hg38)
# seu <- LinkPeaks(seu, peak.assay = "ATAC", genes.use = c("ETS1", "COL1A1"))

library("trackViewer")

library(scran) 
atac <- altExp(filtered)
atac <- logNormCounts(atac)
mi <- scoreMarkers(atac, filtered$clusters)
mi[[6]]
library(scater)
chosen <- mi[[6]]
ordered <- chosen[order(chosen$rank.logFC.cohen),]
topRanked <- ordered[ordered$rank.logFC.cohen<=20,]
atac$clusters <- filtered$clusters
plotGroupedHeatmap(atac, features=rownames(topRanked), group="clusters", center=TRUE, show_rownames = FALSE, annotation_col=data.frame(celltype=filtered$celltype))

#atac_cafs <- atac[,filtered$clusters %in% c(6, 10)]
atac_cafs <- atac[,filtered$celltype == "CAFs"]
mi2 <- scoreMarkers(atac_cafs, atac_cafs$clusters)
mi2
chosen <- unique(rbind(mi2[["6"]], mi2[["10"]], mi2[["5"]], mi2[["2"]]))
ordered <- chosen[order(chosen$rank.AUC),]
topRanked <- ordered[ordered$rank.AUC<=5,]
plotGroupedHeatmap(atac_cafs, features=unique(rownames(topRanked)), group="clusters", center=TRUE, show_rownames = FALSE, zlim = c(-.1, .1))
plotGroupedHeatmap(atac_cafs, features=unique(rownames(topRanked)), group="clusters", center=FALSE, show_rownames = FALSE, zlim = c(0, .2))
c("COL1A1", "ETS1") %in% rranno[rownames(topRanked),]$SYMBOL
rranno[rownames(ordered)[1],]

# atac_epi <- atac[,filtered$celltype == "Cancer Epithelial"]
# mi3 <- scoreMarkers(atac_epi, atac_epi$clusters)
# mi3
# chosen <- rbind(mi3[["1"]], mi3[["7"]])
# ordered <- chosen[order(chosen$rank.logFC.cohen),]
# topRanked <- ordered[ordered$rank.logFC.cohen<=20,]
# plotGroupedHeatmap(atac_epi, features=rownames(topRanked), group="clusters", center=TRUE, show_rownames = FALSE)
# c("COL1A1", "ETS1") %in% rranno[rownames(topRanked),]$SYMBOL
# rranno[rownames(ordered)[1],]
# plotHeatmap(atac_epi, features=unique(rownames(topRanked)), colour_columns_by = "clusters")


reducedDim(atac, "TSNE") <- reducedDim(filtered, "TSNE")
plotTSNE(atac, colour_by = names(rranno[which(rranno$SYMBOL == "COL1A1"),]))
scater::plotHighestExprs(atac)
scater::plotDots(atac, names(rranno[which(rranno$SYMBOL == "COL1A1"),]), group = "clusters") + ggtitle("COL1A1")
scater::plotDots(atac, names(rranno[which(rranno$SYMBOL == "ETS1"),]), group = "clusters") + ggtitle("ETS1")

library(Gviz)
library(rtracklayer)
library(trackViewer)
rranno <- keepStandardChromosomes(rranno)
dt <- DataTrack(range=rranno[which(rranno$SYMBOL == "ETS1"),],
                genome="hg38", type="hist", name="ETS1")
plotTracks(dt)

counts(filteredpb[names(rranno[which(rranno$SYMBOL == "ETS1"),]),])
save.image()
