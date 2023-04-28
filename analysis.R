# Downloading dataset
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6543nnn/GSM6543828/suppl/GSM6543828_HBr_15_A.tar.gz
# 
library(TENxMultiomeTools)
# dir.create("outs")
# sce <- read10xMultiome("HBr_15_A/")

library(Matrix)
mat <- readMM("HBr_15_A/outs/filtered_feature_bc_matrix/matrix.mtx")
dim(mat)
bc <- read.table("HBr_15_A/outs/filtered_feature_bc_matrix/barcodes.tsv", sep="\t")
feat <- read.table("HBr_15_A/outs/filtered_feature_bc_matrix/features.tsv", sep="\t")

sce <- SingleCellExperiment(assays=list("counts"=mat), rowData=feat)
colnames(sce ) <- bc$V1
rownames(sce) <-  feat$V1
sce <- splitAltExps(sce, feat$V3, ref="Gene Expression")

library("EnsDb.Hsapiens.v86")
ensdb <- EnsDb.Hsapiens.v86
seqlevelsStyle(ensdb) <- "UCSC"
library(Signac)
library(Seurat)
annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotation) <- "hg38"
sce <- scuttle::addPerCellQC(sce)
library(scuttle)
colData(sce)
df <- perCellQCFilters(sce, sub.fields=c("altexps_Peaks_detected", 
                            "altexps_Peaks_sum", "altexps_Peaks_percent"))
summary(df$discard)
sce <- sce[,!df$discard]
sce <- logNormCounts(sce)
sce <- computeDoublets(sce, method="griffiths")
sce <- sce[,colData(sce)$dbl.calls=="singlet"]
library(scater)
sce <- runPCA(sce)
sce <- runTSNE(sce)
plotTSNE(sce)
library(SingleR)

# ref=celldex::HumanPrimaryCellAtlasData()
refmd <- read.csv("mereu_cancer_ref/TICAtlas_metadata.csv")
ref <- read.csv("mereu_cancer_ref/TIC")
sceref <- SingleCellExperiment(assays=list(counts=ref), colData=refmd)
colnames(sceref) <- rownames(refmd)
# rownames(sceref) <- TBD
scuttle::aggregateAcrossCells(sceref, ids="")
colData(ref)
rownames(sce) <- rowData(sce)$V2
dflabs <- SingleR(test=sce, ref=ref, labels=ref$label.main)

table(dflabs$pruned.labels)
colData(sce)$labels <- dflabs$pruned.labels
colors <- palette.colors(n=17,palette = "Okabe-Ito") #colorblinded palette
library(Polychrome)
pal <- createPalette(17, seedcolors=c("#ff0000", "#00ff00", "#0000ff"))
pal<- Polychrome::glasbey.colors(n=18)
pal<-pal[-1]
names(pal) <- NULL
sce2 <- sce
nms<-names(table(sce2$labels)[which(table(sce2$labels)>50)])
sce2 <- sce[,sce$labels%in%nms]
plotTSNE(sce2, colour_by="labels",) +
    scale_color_manual(values=pal)
    # scale_color_viridis_d(option="magma")


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
scebcrefagg <- aggregateData(scebcref, assay="counts", fun="sum", by=c("celltype_major", "Patient"))

# scebcrefagg <- logNorm
# dflabs <- SingleR(test=sce, ref=ref, labels=ref$label.main)

