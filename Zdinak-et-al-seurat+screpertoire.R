#Analysis of scRNAseq data for Zdinak et al
#Final version 11/17/2022 Alok V. Joglekar

#Loading Seurat and all the associated packages
library(Seurat)
library(cowplot)
library(patchwork)
library(tidyr)
library(scRepertoire)
library(dplyr)
library(circlize)
library(scales)
library(data.table)
library(magrittr)
library(monocle3)
library(SeuratData) 
library(SeuratWrappers)
library(DESeq2)
library(EnhancedVolcano)
library(hexbin)
library(ggplot2)
library(viridis)
library(SingleCellExperiment)
library(Polychrome)
library(ggbeeswarm)
library(ggthemes)
library(enrichR)
library(stringr)
library(tidyverse)
library(ggridges)

setwd('INSERT THE WD PATH HERE')
#Copied hto_mtx, umi_mtx, and filtered_contig_annotation files from AJ-32/36/39 experiments
#Renamed all the files with their respective expt names

#####AJ-36 - experiment 1 - 6 week mice
AJ36_umis <- readRDS("AJ36_umi_mtx.RDS")
AJ36_htos <- readRDS("AJ36_hto_mtx.RDS")
joint.bcs <- intersect(colnames(AJ36_umis), colnames(AJ36_htos))
AJ36_umis <- AJ36_umis[, joint.bcs]
AJ36_htos <- as.matrix(AJ36_htos[, joint.bcs])
AJ36_hash <- CreateSeuratObject(counts=AJ36_umis, project = "6W")
AJ36_csv <- read.csv("AJ36_filtered_contig_annotations.csv", stringsAsFactors = F)
#normalizing the original Seurat object
AJ36_hash <- NormalizeData(AJ36_hash)
AJ36_hash <- ScaleData(AJ36_hash, features = VariableFeatures(AJ36_hash))
AJ36_hash[["HTO"]] <- CreateAssayObject(counts = AJ36_htos)
AJ36_hash_renamed_v2 <- RenameCells(AJ36_hash, add.cell.id = "all_all")
head(AJ36_hash_renamed_v2[[]])
#combining TCRs
AJ36_vdj <- combineTCR(AJ36_csv, samples = "all", ID = "all", cells = "T-AB")
rm(AJ36_vdj_gex_v5)
#AJ36_sce <- as.SingleCellExperiment(AJ36_hash_renamed_v2)
#rm(AJ36_vdj)
#AJ36_vdj <- combineTCR(AJ36_csv, samples = "all", ID = "all", cells = "T-AB")

#AJ36_vdj_gex_v6 <- combineExpression(AJ36_vdj, AJ36_sce, cloneCall="aa", cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
AJ36_vdj_gex_v4 <- combineExpression(AJ36_vdj, AJ36_hash_renamed_v2, cloneCall="aa", proportion=FALSE, cloneType = c(Single = 1, Small = 9, Medium = 100, Large = 1000, Hyperexpanded=10000 ))
#rm(AJ36_vdj_gex_v4)
#LOOKS GOOD
saveRDS(AJ36_vdj_gex_v4, file="AJ36_vdj_gex_v4_042421_6W.rds")
### THIS WORKS -- now all clonotypes are annotated as 'single, small, medium'
#head(AJ36_vdj_gex_v4$cloneType, 20)


#####AJ-31 - experiment 2 - 8 week mice
AJ31_umis <- readRDS("AJ31_umi_mtx.RDS")
AJ31_htos <- readRDS("AJ31_hto_mtx.RDS")
joint.bcs <- intersect(colnames(AJ31_umis), colnames(AJ31_htos))
AJ31_umis <- AJ31_umis[, joint.bcs]
AJ31_htos <- as.matrix(AJ31_htos[, joint.bcs])
AJ31_hash <- CreateSeuratObject(counts=AJ31_umis, project = "8W")
AJ31_csv <- read.csv("AJ31_filtered_contig_annotations.csv", stringsAsFactors = F)
#normalizing the original Seurat object
AJ31_hash <- NormalizeData(AJ31_hash)
AJ31_hash <- ScaleData(AJ31_hash, features = VariableFeatures(AJ31_hash))
AJ31_hash[["HTO"]] <- CreateAssayObject(counts = AJ31_htos)
AJ31_hash_renamed_v2 <- RenameCells(AJ31_hash, add.cell.id = "all_all")
head(AJ31_hash_renamed_v2[[]])
#combining TCRs
AJ31_vdj <- combineTCR(AJ31_csv, samples = "all", ID = "all", cells = "T-AB")
AJ31_vdj_gex_v4 <- combineExpression(AJ31_vdj, AJ31_hash_renamed_v2, cloneCall="aa", proportion=FALSE, cloneType = c(Single = 1, Small = 9, Medium = 100, Large = 1000, Hyperexpanded=10000 ))
#LOOKS GOOD
head(AJ31_vdj_gex_v4$cloneType, 20)
saveRDS(AJ31_vdj_gex_v4, file="AJ31_vdj_gex_v4_042421_8W.rds")

#####AJ-39 - experiment 3 - 10 week mice
AJ39_umis <- readRDS("AJ39_umi_mtx.RDS")
AJ39_htos <- readRDS("AJ39_hto_mtx.RDS")
joint.bcs <- intersect(colnames(AJ39_umis), colnames(AJ39_htos))
AJ39_umis <- AJ39_umis[, joint.bcs]
AJ39_htos <- as.matrix(AJ39_htos[, joint.bcs])
AJ39_hash <- CreateSeuratObject(counts=AJ39_umis, project = "10W")
AJ39_csv <- read.csv("AJ39_filtered_contig_annotations.csv", stringsAsFactors = F)
#normalizing the original Seurat object
AJ39_hash <- NormalizeData(AJ39_hash)
AJ39_hash <- ScaleData(AJ39_hash, features = VariableFeatures(AJ39_hash))
AJ39_hash[["HTO"]] <- CreateAssayObject(counts = AJ39_htos)
AJ39_hash_renamed_v2 <- RenameCells(AJ39_hash, add.cell.id = "all_all")
head(AJ39_hash_renamed_v2[[]])
#combining TCRs
AJ39_vdj <- combineTCR(AJ39_csv, samples = "all", ID = "all", cells = "T-AB")
AJ39_vdj_gex_v4 <- combineExpression(AJ39_vdj, AJ39_hash_renamed_v2, cloneCall="aa", proportion=FALSE, cloneType = c(Single = 1, Small = 9, Medium = 100, Large = 1000, Hyperexpanded=10000 ))
#LOOKS GOOD
head(AJ39_vdj_gex_v4$cloneType, 20)
saveRDS(AJ39_vdj_gex_v4, file="AJ39_vdj_gex_v4_042421_10W.rds")

## HTO-demultiplexing BEFORE merging
#AJ-39
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
AJ39_vdj_gex_v4 <- NormalizeData(AJ39_vdj_gex_v4, assay = "HTO", normalization.method = "CLR")
AJ39_vdj_gex_v4 <- HTODemux(AJ39_vdj_gex_v4, assay = "HTO", positive.quantile = 0.99)
table(AJ39_vdj_gex_v4$HTO_classification.global)
Idents(AJ39_vdj_gex_v4) <- "HTO_maxID"
Idents(AJ39_vdj_gex_v4) <- "HTO_classification.global"
AJ39_vdj_gex_v4.subset <- subset(AJ39_vdj_gex_v4, idents = "Negative", invert = TRUE)
AJ39_vdj_gex_singlet <- subset(AJ39_vdj_gex_v4, idents = "Singlet")
#extracted singlets

#AJ-31
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
AJ31_vdj_gex_v4 <- NormalizeData(AJ31_vdj_gex_v4, assay = "HTO", normalization.method = "CLR")
AJ31_vdj_gex_v4 <- HTODemux(AJ31_vdj_gex_v4, assay = "HTO", positive.quantile = 0.99)
table(AJ31_vdj_gex_v4$HTO_classification.global)
Idents(AJ31_vdj_gex_v4) <- "HTO_maxID"
Idents(AJ31_vdj_gex_v4) <- "HTO_classification.global"
AJ31_vdj_gex_v4.subset <- subset(AJ31_vdj_gex_v4, idents = "Negative", invert = TRUE)
AJ31_vdj_gex_singlet <- subset(AJ31_vdj_gex_v4, idents = "Singlet")
#extracted singlets

#AJ-36
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
AJ36_vdj_gex_v4 <- NormalizeData(AJ36_vdj_gex_v4, assay = "HTO", normalization.method = "CLR")
AJ36_vdj_gex_v4 <- HTODemux(AJ36_vdj_gex_v4, assay = "HTO", positive.quantile = 0.99)
table(AJ36_vdj_gex_v4$HTO_classification.global)
Idents(AJ36_vdj_gex_v4) <- "HTO_maxID"
Idents(AJ36_vdj_gex_v4) <- "HTO_classification.global"
AJ36_vdj_gex_v4.subset <- subset(AJ36_vdj_gex_v4, idents = "Negative", invert = TRUE)
AJ36_vdj_gex_singlet <- subset(AJ36_vdj_gex_v4, idents = "Singlet")
#extracted singlets

#Merging the three Seurat objects
AJ42_singlet <- merge(AJ36_vdj_gex_singlet, y=c(AJ31_vdj_gex_singlet, AJ39_vdj_gex_singlet), add.cell.ids= c("6W", "8W", "10W"), project = "AJ_42", merge.data = TRUE)
AJ42_singlet
head(colnames(AJ42_singlet))
tail(colnames(AJ42_singlet))
table(AJ42_singlet$orig.ident)
saveRDS(AJ42_singlet, file='AJ42_singlet_010322')
AJ42_singlet <- readRDS("AJ42_singlet_010322")
#start from this point after further analysis

### Adding more QC steps -- removing cells with >5% mt transcripts
AJ42_singlet[["percent.mt"]] <- PercentageFeatureSet(AJ42_singlet, pattern = "^mt-")
VlnPlot(AJ42_singlet, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
#looks good
QCplot1 <- FeatureScatter(AJ42_singlet, feature1="nCount_RNA", feature2="percent.mt")
QCplot2 <- FeatureScatter(AJ42_singlet, feature1="nCount_RNA", feature2="nFeature_RNA")
QCplot1 + QCplot2
AJ42_singlet_QC <- subset (AJ42_singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
saveRDS(AJ42_singlet_QC, file='AJ42_singlet_QC_010322')

## Re-reading the AJ_42_single_QC object as AJ42_singlet
rm(AJ42_singlet, AJ42_singlet_QC)
rm(QCplot1, QCplot2)
AJ42_singlet <- readRDS("AJ42_singlet_QC_010322")

#Normalizing data within each project
AJ42_list_project <- SplitObject(AJ42_singlet, split.by = "orig.ident")
AJ42_list_project <- lapply(X = AJ42_list_project, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
#Find anchors
AJ42_immune.anchors <- FindIntegrationAnchors(object.list = AJ42_list_project, dims = 1:20)
AJ42_immune.combined <- IntegrateData(anchorset = AJ42_immune.anchors, dims = 1:20)
DefaultAssay(AJ42_immune.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
AJ42_immune.combined <- ScaleData(AJ42_immune.combined, verbose = FALSE)
AJ42_immune.combined <- RunPCA(AJ42_immune.combined, npcs = 30, verbose = FALSE)
DimHeatmap(AJ42_immune.combined, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(AJ42_immune.combined)
# t-SNE and Clustering
AJ42_immune.combined <- RunUMAP(AJ42_immune.combined, reduction = "pca", dims = 1:20)
AJ42_immune.combined <- FindNeighbors(AJ42_immune.combined, reduction = "pca", dims = 1:20)
AJ42_immune.combined <- FindClusters(AJ42_immune.combined, resolution = 0.5)
# Visualization
DimPlot(AJ42_immune.combined, reduction = "umap", group.by = "HTO_maxID")
DimPlot(AJ42_immune.combined, reduction = "umap", label = TRUE)
AJ42_immune.combined$orig.ident <- factor(x = AJ42_immune.combined$orig.ident, levels = c("6W", "8W", "10W"))
DimPlot(AJ42_immune.combined, reduction = "umap", split.by = "orig.ident")
AJ42_immune.combined@meta.data$new_ident <- paste(AJ42_immune.combined@meta.data$orig.ident, AJ42_immune.combined@meta.data$HTO_maxID, sep="_", collapse = NULL)
#Here, new identities are a concatenation of the Week+Mouse
#sorting the object with new_idents
AJ42_immune.combined$new_ident <- factor(x = AJ42_immune.combined$new_ident, levels = c("6W_mouse-1", "6W_mouse-2", "6W_mouse-3", "6W_mouse-4", "8W_mouse-1", "8W_mouse-2", "8W_mouse-3", "8W_mouse-4", "10W_mouse-1", "10W_mouse-2", "10W_mouse-3"))
DimPlot(AJ42_immune.combined, reduction = "umap", ncol=4, split.by = "new_ident")

############## adding the screpertoire analysis
# https://ncborcherding.github.io/vignettes/vignette.html
#revisualizing Seurat
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
slot (AJ42_immune.combined, "meta.data")$cloneType <- factor(slot(AJ42_immune.combined, "meta.data")$cloneType, 
                                                             levels = c("Hyperexpanded (1000 < X <= 10000)", "Large (100 < X <= 999)", 
                                                                        "Medium (9 < X <= 100)", "Small (1 < X <= 9)", 
                                                                        "Single (0 < X <= 1)", NA))

DimPlot(AJ42_immune.combined, split.by = "orig.ident", group.by = "cloneType") + scale_color_manual(values=colorblind_vector(4), na.value="grey")
DimPlot(AJ42_immune.combined, split.by = "orig.ident")
DimPlot(AJ42_immune.combined, group.by = "cloneType") + scale_color_manual(values=colorblind_vector(4), na.value="grey")

occupiedscRepertoire(AJ42_immune.combined, x.axis = "orig.ident") + theme(title = NULL, axis.text.x = element_text(angle = 45, hjust = 1))
AJ42_immune.combined_byweek <- expression2List(AJ42_immune.combined, group = "orig.ident")
clonalProportion(AJ42_immune.combined_byweek, split=c(10, 50, 100),cloneCall = "aa")+ theme(title = NULL, axis.text.x = element_text(angle = 45, hjust = 1))
AJ42_immune.combined_bymouse <- expression2List(AJ42_immune.combined, group = "new_ident")
clonalProportion(AJ42_immune.combined_bymouse, split=c(10, 30, 100, 300, 1000, 3000),cloneCall = "aa")+ theme(title = NULL, axis.text.x = element_text(angle = 45, hjust = 1))

###################### Subsetting on CD4s
DefaultAssay(AJ42_immune.combined) <- "RNA"
pA <- FeatureScatter(AJ42_immune.combined, feature1 = "Cd8b1", feature2 = "Cd4")
pA
cd4 <- CellSelector(pA,AJ42_immune.combined)
AJ42_cd4 <- subset(cd4, idents = "SelectedCells")
FeaturePlot(AJ42_immune.combined, features = c("Cd3e", "Cd4", "Cd8b1", "Foxp3"))
FeaturePlot(AJ42_cd4, features = c("Cd3e", "Cd4", "Cd8b1", "Foxp3"))
#looks good
saveRDS(AJ42_cd4, file="AJ42-cd4-gated-010322")

#Saved
#Heatmap by cluster
DefaultAssay(AJ42_cd4) <- "RNA"
AJ42_cd4 <- ScaleData(object = AJ42_cd4, features = rownames(AJ42_cd4))
#FInding all the markers for each cluster:
AJ42_cd4_markers <- FindAllMarkers(AJ42_cd4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Heatmap by top 
AJ42_cd4_markers_top10 <- AJ42_cd4_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AJ42_cd4, features = AJ42_cd4_markers_top10$gene) + NoLegend()
write.csv(AJ42_cd4_markers %>% group_by(cluster), "AJ42_cd4_markers.csv")

#Annotate clusters
Idents(AJ42_cd4) <- "seurat_clusters"

#scRepertoire by cluster
DimPlot(AJ42_cd4, group.by = "cloneType") + scale_color_manual(values=palette5(4), na.value="grey")
#saved
DimPlot(AJ42_cd4, group.by = "cloneType", split.by = "new_ident", ncol=4) + scale_color_manual(values=palette5(4), na.value="grey")
occupiedscRepertoire(AJ42_cd4, x.axis = "orig.ident") + theme(title = NULL, axis.text.x = element_text(angle = 45, hjust = 1))
occupiedscRepertoire(AJ42_cd4, x.axis = "seurat_clusters") + theme(title = NULL, axis.text.x = element_text(angle = 45, hjust = 1))
table <- occupiedscRepertoire(AJ42_cd4, x.axis = "seurat_clusters", exportTable = TRUE) 
write.csv(table, "table.csv")
p3 <- DimPlot(AJ42_cd4, group.by = "seurat_clusters", label=TRUE) 
p2 <- FeaturePlot(AJ42_immune.combined, features="Cd4")
p1 <- DimPlot(AJ42_immune.combined, group.by = "seurat_clusters", label=TRUE) 
p4 <- DimPlot(AJ42_cd4, group.by = "cloneType") + scale_color_manual(values=palette5(4), na.value="grey")

(p1 | p2) / (p3 | p4)
DimPlot(AJ42_immune.combined, group.by = "seurat_clusters", split.by = "new_ident", ncol=4) 

p6<- DimPlot(AJ42_cd4, group.by = "cloneType") + scale_color_manual(values=colorblind_vector(4), na.value="grey")
p3 | p6

#c0, 3,4,5 have expansion, 6 has some
AJ42_cd4_clonotype_freq <- AJ42_cd4@meta.data %>% as.data.table
write.csv(AJ42_cd4_clonotype_freq, "AJ42_cd4_clonotype_freq.csv")
#Saved

########## Subsetting based on expansion
Idents(AJ42_cd4) <- "cloneType"
#AJ42_cd4_expanded <- subset(AJ42_cd4, idents=c("Medium (9 < X <= 100)", "Small (1 < X <= 9)", 
#                                               "Single (0 < X <= 1)"))

#### DEG based on expansion
#DEGs
DefaultAssay(AJ42_cd4) <- "RNA"
Idents(AJ42_cd4) <- "cloneType"
table(Idents(AJ42_cd4))

#single vs small+medium
cd4_DEGs <- FindMarkers(AJ42_cd4, ident.1 = c( 'Medium (9 < X <= 100)', 'Small (1 < X <= 9)'), ident.2 = "Single (0 < X <= 1)", test.use = "DESeq2", max.cells.per.ident = 500)
#rm(cd8_treatment_deg_deseq)
EnhancedVolcano(cd4_DEGs, lab = rownames(cd4_DEGs), x = 'avg_log2FC', y = 'p_val', title = 'Expanded vs Unexpanded CD4+ T cells',  pCutoff = 10e-5, FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0, drawConnectors = TRUE)
write.csv(cd4_DEGs, "AJ42_cd4_Exp-vs-Unexp_DEG_DESeq2.csv")

#single vs small+medium
cd4_DEGs_SvsM <- FindMarkers(AJ42_cd4, ident.1 = 'Medium (9 < X <= 100)', ident.2 = "Small (1 < X <= 9)", test.use = "DESeq2", max.cells.per.ident = 500)
#rm(cd8_treatment_deg_deseq)
EnhancedVolcano(cd4_DEGs, lab = rownames(cd4_DEGs), x = 'avg_log2FC', y = 'p_val', title = 'Medium vs Small Expanded CD4+ T cells',  pCutoff = 10e-5, FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0, drawConnectors = TRUE)
write.csv(cd4_DEGs, "AJ42_cd4_M-vs-S_DEG_DESeq2.csv")
## Use test = DESeq2, LR, wilcox, poisson, bimod 


## Plotting individual clonotypes on psuedotime
AJ42_cd4_TCR_InsChg <- highlightClonotypes(AJ42_cd4, cloneCall= "aa", sequence = c("CAASNNNNAPRF_CASSRDGGRAAEQFF", "CAARNMGYKLTF_CASSQGQGQDTQYF", "CAMRGYGNEKITF_CTCSGGQQDTQYF", "CASLSNNRLTL_CASSSGSQDTQYF"))
InsChg <- DimPlot(AJ42_cd4_TCR_InsChg, group.by = "highlight",  order=TRUE, label=FALSE, pt.size=2) + NoLegend() + scale_color_brewer(palette="Dark2", na.value="grey80") + labs(title = "Ins-ChgA-HIP")

AJ42_cd4_TCR_InsIAPP <- highlightClonotypes(AJ42_cd4, cloneCall= "aa", sequence = c("CAMRGYNQGKLIF_CASSRDSSYEQYF", "CAASNNNNAPRF_CASSMDGGRAETLYF"))
InsIAPP <- DimPlot(AJ42_cd4_TCR_InsIAPP, group.by = "highlight",  order=TRUE, label=FALSE, pt.size=2) + NoLegend() + scale_color_brewer(palette="Dark2", na.value="grey80") + labs(title = "Ins-IAPP-HIP")

AJ42_cd4_TCR_Ins <- highlightClonotypes(AJ42_cd4, cloneCall= "aa", sequence = c( "CAASKVNSGGSNYKLTF_CTCSAGTGSERLFF", "CAASKNYAQGLTF_CASSLVPGYYAEQFF"))
Ins <- DimPlot(AJ42_cd4_TCR_Ins, group.by = "highlight",  order=TRUE, label=FALSE, pt.size=2) + NoLegend() + scale_color_brewer(palette="Dark2", na.value="grey80") + labs(title = "InsB9:23")


