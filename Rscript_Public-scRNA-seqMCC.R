#Public scRNA-seq Merkel Cell Carcinoma x Epigenetics Analysis

#MUST INSTALL PACKAGES FIRST, LIBRARY IMPORTS
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)
library(harmony)
library(clustree)
library(fgsea)
library(msigdbr)
library(Seurat)
library(DESeq2)
library(tibble)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(tidyverse)
library(stringr)
library(EnhancedVolcano)
library(extrafont)

#IMPORT FONTS
font_import(pattern = "Times", prompt = FALSE)
loadfonts(device = "pdf")

#SET WORKING DIRECTORY/FOLDER FILE PATH
#setwd("/Users/sophiahou/Desktop/SRI/Project/Data/GSE226438_RAW")

#LOAD SAMPLES (READ 10X FILES AND COMBINE INTO ONE GENE EXPRESSION MATRIX)
#Pick the correct filepath -- where downloaded GEO data is stored on computer
T98A.data <- Read10X(data.dir = "/Users/sophiahou/Desktop/SRI/Project/Data/GSE226438_RAW/T98A/")
T98B.data <- Read10X(data.dir = "/Users/sophiahou/Desktop/SRI/Project/Data/GSE226438_RAW/T98B/")
T99.data <- Read10X(data.dir = "/Users/sophiahou/Desktop/SRI/Project/Data/GSE226438_RAW/T99/")
T100.data <- Read10X(data.dir = "/Users/sophiahou/Desktop/SRI/Project/Data/GSE226438_RAW/T100/")
T101.data <- Read10X(data.dir = "/Users/sophiahou/Desktop/SRI/Project/Data/GSE226438_RAW/T101/")
T102.data <- Read10X(data.dir = "/Users/sophiahou/Desktop/SRI/Project/Data/GSE226438_RAW/T102/")
T104.data <- Read10X(data.dir = "/Users/sophiahou/Desktop/SRI/Project/Data/GSE226438_RAW/T104/")
T112A.data <- Read10X(data.dir = "/Users/sophiahou/Desktop/SRI/Project/Data/GSE226438_RAW//T112A/")
T112B.data <- Read10X(data.dir = "/Users/sophiahou/Desktop/SRI/Project/Data/GSE226438_RAW/T112B/")
#Skip T96 and T111 (treated), filtered to only treatment-naive

#CREATE SEURAT OBJECTS FOR EACH SAMPLE
T98A <- CreateSeuratObject(
  counts = T98A.data,
  project = "T98A",
  min.cells = 3, #Removes genes detected in fewer than 3 cells (noise).
  min.features = 200) #Removes cells expressing fewer than 200 genes (empty droplets or dead cells).

T98B <- CreateSeuratObject(
  counts = T98B.data,
  project = "T98B",
  min.cells = 3,
  min.features = 200)

T99 <- CreateSeuratObject(
  counts = T99.data,
  project = "T99",
  min.cells = 3,
  min.features = 200)

T100 <- CreateSeuratObject(
  counts = T100.data,
  project = "T100",
  min.cells = 3, 
  min.features = 200)

T101 <- CreateSeuratObject(
  counts = T101.data,
  project = "T101",
  min.cells = 3,
  min.features = 200)

T102 <- CreateSeuratObject(
  counts = T102.data,
  project = "T102",
  min.cells = 3,
  min.features = 200)

T104 <- CreateSeuratObject(
  counts = T104.data,
  project = "T104",
  min.cells = 3,
  min.features = 200)

T112A <- CreateSeuratObject(
  counts = T112A.data,
  project = "T112A",
  min.cells = 3,
  min.features = 200)

T112B <- CreateSeuratObject(
  counts = T112B.data,
  project = "T112B",
  min.cells = 3,
  min.features = 200)

#ADD METADATA LABELS TO SAMPLE OBJECTS: virus (VP/VN), treatment (naive), site (lymph node, gland, or skin), sex (female/male), status (primary/metastasis)
#T98A
T98A$sample <- "T98A"
T98A$virus <- "VN"
T98A$treatment <- "Naive"
T98A$site <- "Skin"
T98A$sex <- "Male"
T98A$status <- "Primary"

#T98B
T98B$sample <- "T98B"
T98B$virus <- "VN"
T98B$treatment <- "Naive"
T98B$site <- "Lymph_Node"
T98B$sex <- "Male"
T98B$status <- "Metastasis"

#T99
T99$sample <- "T99"
T99$virus <- "VP"
T99$treatment <- "Naive"
T99$site <- "Lymph_Node"
T99$sex <- "Male"
T99$status <- "Metastasis"

#T100
T100$sample <- "T100"
T100$virus <- "VP"
T100$treatment <- "Naive"
T100$site <- "Skin"
T100$sex <- "Female"
T100$status <- "Primary"

#T101
T101$sample <- "T101"
T101$virus <- "VP"
T101$treatment <- "Naive"
T101$site <- "Skin"
T101$sex <- "Male"
T101$status <- "Primary"

#T102
T102$sample <- "T102"
T102$virus <- "VN"
T102$treatment <- "Naive"
T102$site <- "Lymph_Node"
T102$sex <- "Male"
T102$status <- "Metastasis"

#T104
T104$sample <- "T104"
T104$virus <- "VN"
T104$treatment <- "Naive"
T104$site <- "Parotid_Gland"
T104$sex <- "Female"
T104$status <- "Metastasis"

#T112A
T112A$sample <- "T112A"
T112A$virus <- "VP"
T112A$treatment <- "Naive"
T112A$site <- "Skin"
T112A$sex <- "Female"
T112A$status <- "Primary"

#T112B
T112B$sample <- "T112B"
T112B$virus <- "VP"
T112B$treatment <- "Naive"
T112B$site <- "Lymph_Node"
T112B$sex <- "Female"
T112B$status <- "Metastasis"

#TUMOR T98A: QC
dim(T98A)

#MITOCHONDRIAL PERCENTAGE: avoid high percent.mt
T98A[["percent.mt"]] <- PercentageFeatureSet(
  T98A,
  pattern = "^MT-")

#VIOLIN PLOTS: visualize QC
#nFeature_RNA: number of genes detected in each cell
#nCount_RNA: total RNA molecules (UMIs) detected
#percent.mt: percent mitochondrial DNA
VlnPlot(
  T98A,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts")

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T98A,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T98A,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#MAD CALCULATION
#MEDIANS
feature.med <- median(T98A$nFeature_RNA)
count.med <- median(T98A$nCount_RNA)
mt.med <- median(T98A$percent.mt)
#MADS
feature.mad <- mad(T98A$nFeature_RNA)
count.mad <- mad(T98A$nCount_RNA)
mt.mad <- mad(T98A$percent.mt)
#THRESHOLDS
feature.low <-feature.med - 3*feature.mad
feature.high <-feature.med + 3*feature.mad
count.low <- count.med - 3*count.mad
count.high <- count.med + 3*count.mad
mt.high <- mt.med + 3*mt.mad

#FILTER
T98A <- subset(
  T98A,
  subset = nFeature_RNA > feature.low &
    nFeature_RNA < feature.high &
    nCount_RNA > count.low &
    nCount_RNA < count.high &
    percent.mt < mt.high)

#VIOLIN PLOTS
VlnPlot(
  T98A,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3)

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T98A,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T98A,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#How many cells remain after QC
dim(T98A) 

#TUMOR 98B
dim(T98B) 

#MITOCHONDRIAL PERCENTAGE
T98B[["percent.mt"]] <- PercentageFeatureSet(
  T98B,
  pattern = "^MT-")

#VIOLIN PLOTS
VlnPlot(
  T98B,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts")

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T98B,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T98B,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#MAD CALCULATION
#MEDIANS
feature.med <- median(T98B$nFeature_RNA)
count.med <- median(T98B$nCount_RNA)
mt.med <- median(T98B$percent.mt)
#MADS
feature.mad <- mad(T98B$nFeature_RNA)
count.mad <- mad(T98B$nCount_RNA)
mt.mad <- mad(T98B$percent.mt)
#THRESHOLDS
feature.low <-feature.med - 3*feature.mad
feature.high <-feature.med + 3*feature.mad
count.low <- count.med - 3*count.mad
count.high <- count.med + 3*count.mad
mt.high <- mt.med + 3*mt.mad

#FILTER
T98B <- subset(
  T98B,
  subset = nFeature_RNA > feature.low &
    nFeature_RNA < feature.high &
    nCount_RNA > count.low &
    nCount_RNA < count.high &
    percent.mt < mt.high)

#VIOLIN PLOTS
VlnPlot(
  T98B,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3)

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T98B,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T98B,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#How many cells remain after QC
dim(T98B) 

#TUMOR T99: QC
dim(T99) 

#MITOCHONDRIAL PERCENTAGE:\
T99[["percent.mt"]] <- PercentageFeatureSet(
  T99,
  pattern = "^MT-")

#VIOLIN PLOTS
VlnPlot(
  T99,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts")

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T99,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T99,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#MAD CALCULATION
#MEDIANS
feature.med <- median(T99$nFeature_RNA)
count.med <- median(T99$nCount_RNA)
mt.med <- median(T99$percent.mt)
#MADS
feature.mad <- mad(T99$nFeature_RNA)
count.mad <- mad(T99$nCount_RNA)
mt.mad <- mad(T99$percent.mt)
#THRESHOLDS
feature.low <-feature.med - 3*feature.mad
feature.high <-feature.med + 3*feature.mad
count.low <- count.med - 3*count.mad
count.high <- count.med + 3*count.mad
mt.high <- mt.med + 3*mt.mad

#FILTER
T99 <- subset(
  T99,
  subset = nFeature_RNA > feature.low &
    nFeature_RNA < feature.high &
    nCount_RNA > count.low &
    nCount_RNA < count.high &
    percent.mt < mt.high)

#VIOLIN PLOTS
VlnPlot(
  T99,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3)

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T99,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T99,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#How many cells remain after QC
dim(T99) 

#TUMOR T100: QC
dim(T100) 

#MITOCHONDRIAL PERCENTAGE:\
T100[["percent.mt"]] <- PercentageFeatureSet(
  T100,
  pattern = "^MT-")

#VIOLIN PLOTS
VlnPlot(
  T100,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts")

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T100,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T100,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#MAD CALCULATION
#MEDIANS
feature.med <- median(T100$nFeature_RNA)
count.med <- median(T100$nCount_RNA)
mt.med <- median(T100$percent.mt)
#MADS
feature.mad <- mad(T100$nFeature_RNA)
count.mad <- mad(T100$nCount_RNA)
mt.mad <- mad(T100$percent.mt)
#THRESHOLDS
feature.low <-feature.med - 3*feature.mad
feature.high <-feature.med + 3*feature.mad
count.low <- count.med - 3*count.mad
count.high <- count.med + 3*count.mad
mt.high <- mt.med + 3*mt.mad

#FILTER
T100 <- subset(
  T100,
  subset = nFeature_RNA > feature.low &
    nFeature_RNA < feature.high &
    nCount_RNA > count.low &
    nCount_RNA < count.high &
    percent.mt < mt.high)

#VIOLIN PLOTS
VlnPlot(
  T100,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3)

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T100,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T100,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#How many cells remain after QC
dim(T100) 

#TUMOR T101: QC
dim(T101) 

#MITOCHONDRIAL PERCENTAGE:\
T101[["percent.mt"]] <- PercentageFeatureSet(
  T101,
  pattern = "^MT-")

#VIOLIN PLOTS
VlnPlot(
  T101,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts")

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T101,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T101,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#MAD CALCULATION
#MEDIANS
feature.med <- median(T101$nFeature_RNA)
count.med <- median(T101$nCount_RNA)
mt.med <- median(T101$percent.mt)
#MADS
feature.mad <- mad(T101$nFeature_RNA)
count.mad <- mad(T101$nCount_RNA)
mt.mad <- mad(T101$percent.mt)
#THRESHOLDS
feature.low <-feature.med - 3*feature.mad
feature.high <-feature.med + 3*feature.mad
count.low <- count.med - 3*count.mad
count.high <- count.med + 3*count.mad
mt.high <- mt.med + 3*mt.mad

#FILTER
T101 <- subset(
  T101,
  subset = nFeature_RNA > feature.low &
    nFeature_RNA < feature.high &
    nCount_RNA > count.low &
    nCount_RNA < count.high &
    percent.mt < mt.high)

#VIOLIN PLOTS
VlnPlot(
  T101,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3)

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T101,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T101,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#How many cells remain after QC
dim(T101) 

#TUMOR T102
#MITOCHONDRIAL PERCENTAGE:\
T102[["percent.mt"]] <- PercentageFeatureSet(
  T102,
  pattern = "^MT-")

#VIOLIN PLOTS
VlnPlot(
  T102,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts")

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T102,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T102,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#MAD CALCULATION
#MEDIANS
feature.med <- median(T102$nFeature_RNA)
count.med <- median(T102$nCount_RNA)
mt.med <- median(T102$percent.mt)
#MADS
feature.mad <- mad(T102$nFeature_RNA)
count.mad <- mad(T102$nCount_RNA)
mt.mad <- mad(T102$percent.mt)
#THRESHOLDS
feature.low <-feature.med - 3*feature.mad
feature.high <-feature.med + 3*feature.mad
count.low <- count.med - 3*count.mad
count.high <- count.med + 3*count.mad
mt.high <- mt.med + 3*mt.mad

#FILTER
T102 <- subset(
  T102,
  subset = nFeature_RNA > feature.low &
    nFeature_RNA < feature.high &
    nCount_RNA > count.low &
    nCount_RNA < count.high &
    percent.mt < mt.high)

#VIOLIN PLOTS
VlnPlot(
  T102,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3)

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T102,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T102,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#TUMOR T104: QC
#MITOCHONDRIAL PERCENTAGE:
T104[["percent.mt"]] <- PercentageFeatureSet(
  T104,
  pattern = "^MT-")

#VIOLIN PLOTS
VlnPlot(
  T104,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts")

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T104,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T104,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#MAD CALCULATION
#MEDIANS
feature.med <- median(T104$nFeature_RNA)
count.med <- median(T104$nCount_RNA)
mt.med <- median(T104$percent.mt)
#MADS
feature.mad <- mad(T104$nFeature_RNA)
count.mad <- mad(T104$nCount_RNA)
mt.mad <- mad(T104$percent.mt)
#THRESHOLDS
feature.low <-feature.med - 3*feature.mad
feature.high <-feature.med + 3*feature.mad
count.low <- count.med - 3*count.mad
count.high <- count.med + 3*count.mad
mt.high <- mt.med + 3*mt.mad

#FILTER
T104 <- subset(
  T104,
  subset = nFeature_RNA > feature.low &
    nFeature_RNA < feature.high &
    nCount_RNA > count.low &
    nCount_RNA < count.high &
    percent.mt < mt.high)

#VIOLIN PLOTS
VlnPlot(
  T104,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3)

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T104,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T104,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#How many cells remain after QC
dim(T104) 

#TUMOR T112A: QC
#skipped here because below is filtering for figures again.
dim(T112A) 

#MITOCHONDRIAL PERCENTAGE:\
T112A[["percent.mt"]] <- PercentageFeatureSet(
  T112A,
  pattern = "^MT-")

#VIOLIN PLOTS
VlnPlot(
  T112A,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts")

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T112A,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T112A,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#MAD CALCULATION
#MEDIANS
feature.med <- median(T112A$nFeature_RNA)
count.med <- median(T112A$nCount_RNA)
mt.med <- median(T112A$percent.mt)
#MADS
feature.mad <- mad(T112A$nFeature_RNA)
count.mad <- mad(T112A$nCount_RNA)
mt.mad <- mad(T112A$percent.mt)
#THRESHOLDS
feature.low <-feature.med - 3*feature.mad
feature.high <-feature.med + 3*feature.mad
count.low <- count.med - 3*count.mad
count.high <- count.med + 3*count.mad
mt.high <- mt.med + 3*mt.mad

#FILTER
T112A <- subset(
  T112A,
  subset = nFeature_RNA > feature.low &
    nFeature_RNA < feature.high &
    nCount_RNA > count.low &
    nCount_RNA < count.high &
    percent.mt < mt.high)

#VIOLIN PLOTS
VlnPlot(
  T112A,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3)

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T112A,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T112A,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#How many cells remain after QC
dim(T112A) 


#Representative visualization for QC-filtering, T112A
#Initial
print(dim(T112A))
Idents(T112A) <- "T112A"
T112A[["percent.mt"]] <- PercentageFeatureSet(T112A, pattern = "^MT-")

#Pre-QC plots
vln_pre <- VlnPlot(
  T112A,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  combine = FALSE,
  pt.size = 0.1)
vln_pre[[1]] <- vln_pre[[1]] + 
  ggtitle("Pre-QC: Gene Count") + 
  ylab("Number of Genes") + 
  xlab(NULL) + 
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none", 
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    plot.title = element_text(size = 11))
vln_pre[[2]] <- vln_pre[[2]] + 
  ggtitle("Pre-QC: UMI Count") + 
  ylab("Number of UMIs") + 
  xlab(NULL) + 
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none", 
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    plot.title = element_text(size = 11))
vln_pre[[3]] <- vln_pre[[3]] + 
  ggtitle("Pre-QC: mtRNA %") + 
  ylab("Percentage of mtRNA") + 
  xlab(NULL) + 
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none", 
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    plot.title = element_text(size = 11))
wrap_plots(vln_pre, ncol = 3)
plot1_pre <- FeatureScatter(T112A, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("Pre-QC: UMI vs. Gene Count") +
  xlab("UMI Count") +
  ylab("Gene Count") +
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none", 
    plot.title = element_text(size = 11))
plot2_pre <- FeatureScatter(T112A, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("Pre-QC: UMI vs. mtRNA %") +
  xlab("UMI Count") +
  ylab("% mtRNA") +
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none", 
    plot.title = element_text(size = 11))
plot1_pre + plot2_pre

#MAD
feature.med <- median(T112A$nFeature_RNA)
count.med <- median(T112A$nCount_RNA)
mt.med <- median(T112A$percent.mt)

feature.mad <- mad(T112A$nFeature_RNA)
count.mad <- mad(T112A$nCount_RNA)
mt.mad <- mad(T112A$percent.mt)

feature.low <- feature.med - 3 * feature.mad
feature.high <- feature.med + 3 * feature.mad
count.low <- count.med - 3 * count.mad
count.high <- count.med + 3 * count.mad
mt.high <- mt.med + 3 * mt.mad

T112A <- subset(
  T112A,
  subset = nFeature_RNA > feature.low &
    nFeature_RNA < feature.high &
    nCount_RNA > count.low &
    nCount_RNA < count.high &
    percent.mt < mt.high)

#POST-QC PLOTS
vln_post <- VlnPlot(
  T112A,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  combine = FALSE,
  pt.size = 0.1)
vln_post[[1]] <- vln_post[[1]] + 
  ggtitle("Post-QC: Gene Count") + 
  ylab("Number of Genes") + 
  xlab(NULL) + 
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none", 
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    plot.title = element_text(size = 11))
vln_post[[2]] <- vln_post[[2]] + 
  ggtitle("Post-QC: UMI Count") + 
  ylab("Number of UMIs") + 
  xlab(NULL) + 
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none", 
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    plot.title = element_text(size = 11))
vln_post[[3]] <- vln_post[[3]] + 
  ggtitle("Post-QC: mtRNA %") + 
  ylab("Percentage of mtRNA") + 
  xlab(NULL) + 
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none", 
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    plot.title = element_text(size = 11))
wrap_plots(vln_post, ncol = 3)
plot1_post <- FeatureScatter(T112A, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("Post-QC: UMI vs. Gene Count") +
  xlab("UMI Count") +
  ylab("Gene Count") +
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none", 
    plot.title = element_text(size = 11))
plot2_post <- FeatureScatter(T112A, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("Post-QC: UMI vs. mtRNA %") +
  xlab("UMI Count") +
  ylab("% mtRNA") +
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none", 
    plot.title = element_text(size = 11))
plot1_post + plot2_post



#TUMOR T112B: QC
#MITOCHONDRIAL PERCENTAGE:
T112B[["percent.mt"]] <- PercentageFeatureSet(
  T112B,
  pattern = "^MT-")

#VIOLIN PLOTS
VlnPlot(
  T112B,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts")

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T112B,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T112B,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#MAD CALCULATION
#MEDIANS
feature.med <- median(T112B$nFeature_RNA)
count.med <- median(T112B$nCount_RNA)
mt.med <- median(T112B$percent.mt)
#MADS
feature.mad <- mad(T112B$nFeature_RNA)
count.mad <- mad(T112B$nCount_RNA)
mt.mad <- mad(T112B$percent.mt)
#THRESHOLDS
feature.low <-feature.med - 3*feature.mad
feature.high <-feature.med + 3*feature.mad
count.low <- count.med - 3*count.mad
count.high <- count.med + 3*count.mad
mt.high <- mt.med + 3*mt.mad

#FILTER
T112B <- subset(
  T112B,
  subset = nFeature_RNA > feature.low &
    nFeature_RNA < feature.high &
    nCount_RNA > count.low &
    nCount_RNA < count.high &
    percent.mt < mt.high)

#VIOLIN PLOTS
VlnPlot(
  T112B,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3)

#SCATTER PLOTS
plot1 <- FeatureScatter(
  T112B,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(
  T112B,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt")
plot1 + plot2

#MERGE ALL TUMORS
mcc <- merge(
  x = T98A,
  y = list(T98B,T99,T100,T101,T102,T104,T112A,T112B))

#information
ncol(mcc)
length(Cells(mcc))
table(mcc$sample, mcc$virus)

#additional info
dim(mcc)
table(mcc$sample)
table(mcc$virus)
table(mcc$treatment)
table(mcc$site)
View(mcc@meta.data)

#NORMALIZE
mcc <- NormalizeData(
  mcc,
  normalization.method = "LogNormalize",
  scale.factor = 10000)

#FIND VARIABLE GENES
mcc <- FindVariableFeatures(
  mcc,
  selection.method = "vst",
  nfeatures = 3000)

#SCALE
mcc <- ScaleData(mcc)

#PRINCIPAL COMPONENT ANALYSIS (PCA)
#Positive genes: genes with the highest positive loading
#Negative genes: genes with the highest negative loading
#Genes on opposite sides are basically what separate cells along that axis.
mcc <- RunPCA(
  mcc,
  npcs = 50)

#ELBOW PLOT
ElbowPlot(mcc)
VizDimLoadings(
  mcc,
  dims = 1:4,
  reduction = "pca")
DimHeatmap(
  mcc,
  dims = 1:6,
  cells = 500,
  balanced = TRUE)

#PCA COLORED BY 
p1 <- DimPlot(
  mcc,
  reduction = "pca",
  group.by = "sample",
  pt.size = 0.5)
p2 <- DimPlot(
  mcc,
  reduction = "pca",
  group.by = "virus",
  pt.size = 0.5)
p3 <- DimPlot(
  mcc,
  reduction = "pca",
  group.by = "site",
  pt.size = 0.5)
p4 <- DimPlot(
  mcc,
  reduction = "pca",
  group.by = "sex",
  pt.size = 0.5)
p5 <- DimPlot(
  mcc,
  reduction = "pca",
  group.by = "status",
  pt.size = 0.5)
p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 2)

#RUN HARMONY
mcc <- RunHarmony(
  object = mcc,
  group.by.vars = "sample")

#Pre-Harmony Plot
p1 <- DimPlot(
  mcc,
  reduction = "pca",
  group.by = "sample",
  pt.size = 0.3) + 
  ggtitle("Pre-Harmony by Sample") + 
  theme(text = element_text(family = "Times New Roman"))

#Post-Harmony Plot
p2 <- DimPlot(
  mcc,
  reduction = "harmony",
  group.by = "sample",
  pt.size = 0.3) + 
  ggtitle("Post-Harmony by Sample") + 
  theme(text = element_text(family = "Times New Roman"))
p1
p2

#ELBOW PLOT POST-HARMONY
ElbowPlot(
  mcc,
  reduction = "harmony")

#NEIGHBORS
mcc <- FindNeighbors(
  mcc,
  reduction = "harmony",
  dims = 1:10)

#CLUSTERING
mcc <- FindClusters(
  mcc,
  resolution = c(0.2,0.4,0.6,0.8)) #testing multiple

clustree::clustree(mcc)

mcc <- FindClusters(
  mcc,
  resolution = 0.4)

#SET RES
Idents(mcc) <- "RNA_snn_res.0.4"

#UMAP, HARMONY
mcc <- RunUMAP(
  mcc,
  reduction = "harmony",
  dims = 1:10)

#VISUALIZE CLUSTER W RES 0.4
DimPlot(
  mcc,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.4)

#CHECK WHAT DRIVES CLUSTERS
DimPlot(
  mcc,
  reduction = "umap",
  group.by = "sample",
  pt.size = 0.2)
DimPlot(
  mcc,
  reduction = "umap",
  group.by = "virus",
  pt.size = 0.5)
DimPlot(
  mcc,
  reduction = "umap",
  group.by = "status",
  pt.size = 0.5)

#CONFIRM HARMONY
DimPlot(
  mcc,
  reduction = "umap",
  split.by = "sample",
  ncol = 3)

#FIND CLUSTER MARKERS
mcc <- JoinLayers(mcc)

cluster.markers <- FindAllMarkers(
  object = mcc,
  #defaults
  assay = "RNA", 
  only.pos = TRUE,
  min.pct = 0.25, 
  logfc.threshold = 0.25,
  test.use = "wilcox")
head(cluster.markers)

#MARKER TABLE FOR ANNOTATION (10)
top10.markers <- cluster.markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10)
View(top10.markers)

#DOT PLOT W MARKER GENES
DotPlot(
  mcc,
  features = list(
    MCC_Neuroendocrine = c("KRT20", "CHGA","CHGB", "INSM1","SOX2",  "POU4F3"),
    T_cells = c("CD3D","CD3E","TRAC", "IL7R"),
    NK_Cytotoxic = c("NKG7", "GNLY", "GZMB", "PRF1", "KLRD1"), 
    Activated_T = c("TNFRSF4", "TNFRSF18", "CD69", "RGS1"), 
    B_cells = c("MS4A1",  "CD79A", "CD74", "HLA-DRA"),
    Myeloid = c("LYZ", "LST1", "FCER1G", "TYROBP"),
    Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN","LUM"),
    Endothelial = c("PECAM1","VWF","KDR","EMCN"),
    Cycling = c("MKI67", "TOP2A", "BIRC5", "MCM6"), 
    Stress = c("HSPA6", "HSPA1A", "DNAJB1", "FOS", "FOSB"))) + RotatedAxis()
#Extra
DotPlot(
  mcc,
  features = c("KRT20", "CHGA","CHGB","INSM1", 
               "SOX2","POU4F3","CD3D","CD3E","TRAC","MS4A1","CD79A","NKG7",
               "GNLY","LYZ","COL1A1","COL3A1","KRT15", "MAGEA4"),
  group.by = "RNA_snn_res.0.4") +RotatedAxis()

#FEATURE PLOT W MARKER GENES
#MCC/NEUROENDOCRINE
FeaturePlot(
  mcc,
  features = c("KRT20", "CHGA", "CHGB", "INSM1", "SOX2", "POU4F3"),
  reduction = "umap",
  ncol = 3)
#T CELLS
FeaturePlot(
  mcc,
  features = c("CD3D", "CD3E", "TRAC", "IL7R"),
  reduction = "umap",
  ncol = 2)
#NK / CYTOTOXIC
FeaturePlot(
  mcc,
  features = c("NKG7", "GNLY", "GZMB", "PRF1", "KLRD1"),
  reduction = "umap",
  ncol = 3)
#ACTIVATED T CELLS
FeaturePlot(
  mcc,
  features = c("TNFRSF4", "TNFRSF18", "CD69", "RGS1"),
  reduction = "umap",
  ncol = 2)
#B CELLS
FeaturePlot(
  mcc,
  features = c("MS4A1", "CD79A", "CD74", "HLA-DRA"),
  reduction = "umap",
  ncol = 2)
#MYELOID CELLS
FeaturePlot(
  mcc,
  features = c("LYZ", "LST1", "FCER1G", "TYROBP"),
  reduction = "umap",
  ncol = 2)
#FIBROBLASTS
FeaturePlot(
  mcc,
  features = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  reduction = "umap",
  ncol = 3)
#ENDOTHELIAL CELLS
FeaturePlot(
  mcc,
  features = c("PECAM1", "VWF", "KDR", "EMCN"),
  reduction = "umap",
  ncol = 2)
#CYCLING/PROLIFERATION
FeaturePlot(
  mcc,
  features = c("MKI67", "TOP2A", "BIRC5", "MCM6"),
  reduction = "umap",
  ncol = 2)
#STRESS
FeaturePlot(
  mcc,
  features = c("HSPA6", "HSPA1A", "DNAJB1", "FOS", "FOSB"),
  reduction = "umap",
  ncol = 3)

#CLUSTER VISUALIZATION
DimPlot(
  mcc,
  reduction = "umap",
  group.by = "RNA_snn_res.0.4",
  label = TRUE,
  repel = TRUE)

#CLUSTER 6 FURTHER ANALYSIS
FeaturePlot(
  mcc,
  features = c(
    "KRT20","CHGA","INSM1","SOX2","HSPA6","FOS"), #neuroendocrine + stress
  ncol = 3,
  order = TRUE)
#Immune markers
FeaturePlot(
  mcc,
  features = c(
    "CD3D","TRAC","NKG7","LYZ")) #no
FeaturePlot(
  mcc,
  features = c(
    "KRT20","CHGA","INSM1","HSPA6",
    order = TRUE,
    cells.highlight = WhichCells(
      mcc,
      idents = "6")))

#CLUSTER 13 FURTHER ANALYSIS
cluster13.markers <- FindMarkers(
  mcc,
  ident.1 = "13",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox")
head(cluster13.markers, 30)

FeaturePlot(
  mcc,
  features = c("KRT20","CHGA","CHGB","INSM1","SOX2","POU4F3","MKI67","TOP2A","HSPA6","FOS"),
  reduction = "umap",
  ncol = 5)
VlnPlot(
  mcc,
  features = c("KRT20","INSM1","CD3D","CD8A","CCL5","NKG7"),
  idents = "13",
  pt.size = 0)
VlnPlot(
  mcc,
  features = c("KRT20","INSM1","CD3D","CD8A","CCL5","NKG7"),
  idents = c("2","13"),
  pt.size = 0)
FeaturePlot(
  mcc,
  features = c("CD8A","CD3D","CCL5","NKG7","KRT20","INSM1"),
  reduction = "umap",
  ncol = 3,
  order = TRUE)

#CELL TYPE ANNOTATION (1)
#INITIAL CELL-TYPE ANNOTATION
mcc$cell_type <- "Unassigned"
#MCC neuroendocrine
mcc$cell_type[mcc$RNA_snn_res.0.4 %in% c("0", "1")] <- "Neuroendocrine MCC"
#Cycling clusters
mcc$cell_type[mcc$RNA_snn_res.0.4 %in% c("3", "5", "9")] <- "Cycling MCC"
#T cells
mcc$cell_type[mcc$RNA_snn_res.0.4 == "2"] <- "T cell"
#B cells
mcc$cell_type[mcc$RNA_snn_res.0.4 %in% c("4", "12")] <- "B cell"
#Activated B cells
mcc$cell_type[mcc$RNA_snn_res.0.4 == "11"] <- "Activated B cell"
#NK cytotoxic
mcc$cell_type[mcc$RNA_snn_res.0.4 == "7"] <- "NK Cytotoxic"
#Activated T cells
mcc$cell_type[mcc$RNA_snn_res.0.4 == "8"] <- "Activated T cell"
#Fibroblasts
mcc$cell_type[mcc$RNA_snn_res.0.4 == "10"] <- "Fibroblast"
#Myeloid
mcc$cell_type[mcc$RNA_snn_res.0.4 == "14"] <- "Myeloid"
#Cytotoxic/activated immune
mcc$cell_type[mcc$RNA_snn_res.0.4 == "15"] <- "Activated Cytotoxic T/NK"
#Stressed tumor cells
mcc$cell_type[mcc$RNA_snn_res.0.4 == "6"] <- "Stressed MCC"
#Inflamed MCC 
mcc$cell_type[mcc$RNA_snn_res.0.4 == "13"] <- "Inflamed MCC"

#CHECK ANNOTATIONS
table(mcc$cell_type)

table(
  mcc$RNA_snn_res.0.4,
  mcc$cell_type) 

DimPlot(
  mcc,
  reduction = "umap",
  group.by = "cell_type",
  repel = TRUE,
  pt.size = 0.4) + labs(
    title = "TME Cell Types") + 
  theme(text = element_text(family = "Times New Roman"))

#CHECK ANNOTATION AGAINST VIRUS STATUS
DimPlot(
  mcc,
  reduction = "umap",
  group.by = "sample",
  pt.size = 0.3,
  alpha = 0.5) +
  labs(
    title = "TME by Sample") +
  theme(
    text = element_text(family = "Times New Roman"))

DimPlot(
  mcc,
  reduction = "umap",
  group.by = "virus",
  pt.size = 0.3,
  alpha = 0.5) +
  labs(
    title = "TME by Subtype") +
  scale_color_discrete(
    labels = c("VP" = "Virus-positive", "VN" = "Virus-negative")) +
  theme(
    text = element_text(family = "Times New Roman"))

#CELL TYPE MARKER GENE HEATMAP
Idents(mcc) <- "cell_type"
celltype.markers <- FindAllMarkers(
  object = mcc,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox")
top5 <- celltype.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5,
    with_ties = FALSE
  ) %>%
  ungroup()
genes_to_plot <- unique(top5$gene)
mcc <- ScaleData(
  mcc,
  features = genes_to_plot)
mat <- GetAssayData(
  mcc,
  layer = "scale.data"
)[genes_to_plot, ]
loadfonts(quiet = TRUE)
group_vec <- factor(
  mcc$cell_type,
  levels = sort(unique(mcc$cell_type)))
cell_order <- order(group_vec)
mat_ordered <- mat[, cell_order]
group_ordered <- group_vec[cell_order]
celltype.colors <- setNames(
  scales::hue_pal()(length(levels(group_ordered))),
  levels(group_ordered))
bottom_anno <- HeatmapAnnotation(
  Cell_Type = group_ordered,
  which = "column",
  col = list(
    Cell_Type = celltype.colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    title = "Cell Type",
    title_gp = gpar(
      fontfamily = "Times New Roman",
      fontface = "plain",
      fontsize = 11),
    labels_gp = gpar(
      fontfamily = "Times New Roman",
      fontface = "plain",
      fontsize = 10)))
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "#F7F7F7", "#B2182B"))

#Heatmap
Heatmap(
  mat_ordered,
  name = "Scaled Expression",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  column_split = group_ordered,
  column_title = NULL,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(
    fontfamily = "Times New Roman",
    fontsize = 10,
    fontface = "plain"),
  bottom_annotation = bottom_anno,
  column_gap = unit(1.5, "mm"),
  use_raster = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(
      fontfamily = "Times New Roman",
      fontface = "plain",
      fontsize = 11),
    
    labels_gp = gpar(
      fontfamily = "Times New Roman",
      fontface = "plain",
      fontsize = 10)))

#CREATE TUMOR CELL ONLY OBJECT (INCLUDE ONLY MCC, EXCLUDE NON-MCC LIKE IMMUNE)
tumor_cells <- subset(
  mcc,
  subset = cell_type %in% c(
    "Neuroendocrine MCC",
    "Cycling MCC",
    "Stressed MCC",
    "Inflamed MCC"))

table(tumor_cells$cell_type)
table(tumor_cells$virus)
table(tumor_cells$status)
table(tumor_cells$sample)

DimPlot(
  tumor_cells,
  reduction = "umap",
  repel = TRUE)

table(
  mcc$RNA_snn_res.0.4,
  mcc$sample)
table(
  mcc$cell_type,
  mcc$sample)
DimPlot(
  mcc,
  reduction = "umap",
  group.by = "cell_type",
  split.by = "sample",
  ncol = 3)

#REPROCESS TUMOR ONLY OBJECT
tumor_cells <- NormalizeData(tumor_cells)
tumor_cells <- FindVariableFeatures(
  tumor_cells,
  selection.method = "vst",
  nfeatures = 3000)
tumor_cells <- ScaleData(tumor_cells)
tumor_cells <- RunPCA(
  tumor_cells,
  npcs = 30)
tumor_cells <- RunHarmony(
  tumor_cells,
  group.by.vars = "sample")
tumor_cells <- FindNeighbors(
  tumor_cells,
  reduction = "harmony",
  dims = 1:10)
tumor_cells <- FindClusters(
  tumor_cells,
  resolution = 0.4)
tumor_cells <- RunUMAP(
  tumor_cells,
  reduction = "harmony",
  dims = 1:10)

#TUMOR-ONLY CLUSTER ANNOTATION
tumor_markers <- FindAllMarkers(
  tumor_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25)

top10.tumor.markers <- tumor_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10)
View(top10.tumor.markers)

#CLUSTERS
#Check clusters/markers
markers <- FindAllMarkers(
  tumor_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25)

markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10)

DimPlot(
  tumor_cells,
  reduction = "umap",
  group.by = "seurat_clusters",
  repel = TRUE)

tumor_cells <- subset(
  tumor_cells,
  subset = seurat_clusters %notin% c("8","11")) #Removed because non-MCC according to markers.


#check remaining clusters
table(tumor_cells$seurat_clusters)

DimPlot(
  tumor_cells,
  reduction = "umap",
  group.by = "seurat_clusters",
  repel = TRUE)


#ANNOTATION OF NEW TUMOR SUBSET (2)
tumor_cells$Tumor_State <- "Unassigned"
tumor_cells$Tumor_State[
  tumor_cells$seurat_clusters == "0"] <- "Transcriptional maintenance MCC"
tumor_cells$Tumor_State[
  tumor_cells$seurat_clusters == "1"] <- "Differentiated MCC"
tumor_cells$Tumor_State[
  tumor_cells$seurat_clusters == "2"] <- "Stress response MCC"
tumor_cells$Tumor_State[
  tumor_cells$seurat_clusters == "3"] <- "Neuroendocrine signaling MCC"
tumor_cells$Tumor_State[
  tumor_cells$seurat_clusters == "4"] <- "Metabolic MCC"
tumor_cells$Tumor_State[
  tumor_cells$seurat_clusters == "5"] <- "Cycling S-phase MCC"
tumor_cells$Tumor_State[
  tumor_cells$seurat_clusters == "6"] <- "Secretory neuroendocrine MCC"
tumor_cells$Tumor_State[
  tumor_cells$seurat_clusters == "7"] <- "Cycling G2/M MCC"
tumor_cells$Tumor_State[
  tumor_cells$seurat_clusters == "9"] <- "Transitional MCC"
tumor_cells$Tumor_State[
  tumor_cells$seurat_clusters == "10"] <- "Metabolic stress MCC"

table(
  tumor_cells$Tumor_State,
  tumor_cells$virus)

DimPlot(
  tumor_cells,
  group.by = "Tumor_State",
  repel = TRUE) + 
  labs(title = "Tumor Cell State") +
  theme(
    text = element_text(family = "Times New Roman"))

#Check MCC markers are retained
FeaturePlot(
  tumor_cells,
  features = c("KRT20","CHGA","CHGB","INSM1","SOX2","POU4F3","MKI67","TOP2A"),
  ncol = 4,
  order = TRUE)

#CHECK VP vs. VN
DimPlot(
  tumor_cells,
  group.by = "virus") + 
  labs(title = "Tumor Cells by Subtype") +
  scale_color_discrete(
    labels = c("VN" = "Virus-negative", "VP" = "Virus-positive") ) +
  theme(
    text = element_text(family = "Times New Roman"))

DimPlot(
  tumor_cells,
  group.by = "sample",
  repel = TRUE) + 
  labs(title = "Tumor Cells by Sample") +
  theme(
    text = element_text(family = "Times New Roman"))


#CELL STATE HEATMAP
genes_to_plot <- unique(top5$gene)
tumor_cells <- ScaleData(
  tumor_cells,
  features = genes_to_plot)
mat <- GetAssayData(
  tumor_cells,
  layer = "scale.data"
)[genes_to_plot, ]
group_vec <- factor(
  tumor_cells$Tumor_State,
  levels = sort(unique(tumor_cells$Tumor_State)))
cell_order <- order(group_vec)
mat_ordered <- mat[, cell_order]
group_ordered <- group_vec[cell_order]
unique_groups <- levels(group_ordered)
group_colors <- setNames(
  scales::hue_pal()(length(unique_groups)),
  unique_groups)
bottom_anno <- HeatmapAnnotation(
  Tumor_State = group_ordered,
  col = list(
    Tumor_State = group_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    title = "Tumor State",
    title_gp = gpar(
      fontfamily = "Times New Roman",
      fontsize = 11,
      fontface = "plain"),
    labels_gp = gpar(
      fontfamily = "Times New Roman",
      fontsize = 10,
      fontface = "plain")))
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "#F7F7F7", "#B2182B"))
ComplexHeatmap::Heatmap(
  mat_ordered,
  name = "Scaled Expression",
  col = col_fun,
  # Split by tumor state
  column_split = group_ordered,
  column_title = NULL,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(
    fontfamily = "Times New Roman",
    fontsize = 11,
    fontface = "plain"),
  bottom_annotation = bottom_anno,
  heatmap_legend_param = list(
    title_gp = gpar(
      fontfamily = "Times New Roman",
      fontsize = 11,
      fontface = "plain"),
    labels_gp = gpar(
      fontfamily = "Times New Roman",
      fontsize = 10,
      fontface = "plain")),
  column_gap = unit(1.5, "mm"),
  use_raster = FALSE)

#RUN PSEUDOBULK
pseudo_bulk <- AggregateExpression(
  tumor_cells,
  assays = "RNA",
  group.by = "sample",
  return.seurat = TRUE)

pseudo_bulk <- NormalizeData(pseudo_bulk)
pseudo_bulk <- FindVariableFeatures(pseudo_bulk)
pseudo_bulk <- ScaleData(pseudo_bulk)
pseudo_bulk <- RunPCA(pseudo_bulk, npcs = 8)

expr_mat <- GetAssayData(pseudo_bulk, layer = "data")
sample_cor <- cor(as.matrix(expr_mat), method = "pearson")
pheatmap::pheatmap(sample_cor)

#DIFFERENTIAL GENE EXPRESSION (DGE) ANALYSIS 
counts_mat <- GetAssayData(pseudo_bulk, assay = "RNA", layer = "counts")
sample_meta <- pseudo_bulk@meta.data

sample_info <- tumor_cells@meta.data %>%
  select(sample, virus) %>%
  distinct()
sample_meta <- pseudo_bulk@meta.data %>%
  left_join(sample_info, by = "sample")
rownames(sample_meta) <- colnames(pseudo_bulk)

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_mat)),
  colData = sample_meta,
  design = ~ virus)

dds <- DESeq(dds)
res <- results(dds, contrast = c("virus", "VP", "VN"))
res <- res[order(res$padj), ]

sum(res$padj < 0.05, na.rm = TRUE)
head(as.data.frame(res), 20)

#avg_log2FC > 0: higher in VP
#avg_log2FC < 0: higher in VN
#p_val_adj: statistical significance after correction
#pct.1/pct.2: percentage of cells expressing the gene

#DESeq2 volcano pot
EnhancedVolcano(
  as.data.frame(res),
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 3,
  selectLab = c(
    "EZHIP", "EZH2", "BMI1", "SUZ12",
    "SOX2", "ATOH1", "INSM1",
    "CHGA", "CHGB", "KRT20"),
  title = "VP vs VN pseudobulk differential expression")
#top deg genes
head(as.data.frame(res), 30)

#CORE POLYCOMB GENES DGE
#Core polycomb genes DE? No.
prc_genes <- c("EZH2","SUZ12","EED","BMI1","RING1","RNF2")
prc_DE <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  filter(gene %in% prc_genes)
prc_DE

#GENE SET ENRICHMENT ANALYSIS (GSEA)
#Ranked gene list
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[!is.na(res_df$stat), ]
ranks <- res_df$stat
names(ranks) <- res_df$gene
ranks <- sort(ranks, decreasing = TRUE)

#PRC1 GSEA
prc1_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CGP") %>%
  filter(
    grepl(
      "BMI1|RING1|RNF2|CBX7|CBX8|PCGF|PRC1",
      gs_name,
      ignore.case = TRUE)) %>%
  split(x = .$gene_symbol,
        f = .$gs_name)

fgsea_prc1 <- fgsea(
  pathways = prc1_sets,
  stats = ranks,
  minSize = 10,
  maxSize = 500)

fgsea_prc1 <- fgsea_prc1[
  order(fgsea_prc1$padj),]

fgsea_prc1


#PRC2 GSEA
prc2_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CGP") %>%
  filter(
    grepl(
      "EZH2|SUZ12|PRC2",
      gs_name,
      ignore.case = TRUE)) %>%
  split(x = .$gene_symbol,
        f = .$gs_name)
fgsea_prc2 <- fgsea(
  pathways = prc2_sets,
  stats = ranks,
  minSize = 10,
  maxSize = 500)
fgsea_prc2 <- fgsea_prc2[
  order(fgsea_prc2$padj),]
fgsea_prc2

#HALLMARK GSEA (not used)
hallmark_sets <- msigdbr(
  species = "Homo sapiens",
  category = "H") %>%
  split(
    x = .$gene_symbol,
    f = .$gs_name)

hallmark_fgsea <- fgsea(
  pathways = hallmark_sets,
  stats = ranks,
  minSize = 15,
  maxSize = 500)
hallmark_fgsea <- hallmark_fgsea[
  order(hallmark_fgsea$padj),]
hallmark_fgsea


#PATHWAY ENRICHMENT (SIGNIFICANT ONLY) HORIZONTAL BAR PLOT
fgseaRes <- bind_rows(
  fgsea_prc1 %>% mutate(Complex = "PRC1"),
  fgsea_prc2 %>% mutate(Complex = "PRC2"))
polycomb_bar <- fgseaRes %>%
  filter(padj < 0.05) %>%
  filter(grepl("PRC|POLYCOMB|BMI1|EZH2|SUZ12|H3K27", pathway, ignore.case = TRUE)) %>%
  mutate(
    Complex = case_when(
      grepl("PRC1|BMI1", pathway, ignore.case = TRUE) ~ "PRC1",
      TRUE ~ "PRC2"),
    pathway = reorder(pathway, NES))

#HORIZONTAL BAR PLOT
ggplot(polycomb_bar, aes(x = NES, y = pathway, fill = Complex)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("PRC1" = "#F4A234", "PRC2" = "#bd83ec")) +
  coord_cartesian(xlim = c(min(polycomb_bar$NES, na.rm = TRUE), 2.5)) +
  scale_x_continuous(breaks = seq(-3, 2.5, by = 0.5)) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  labs(
    title = "PRC1- and PRC2-Target Gene Signature Enrichment",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    fill = "Complex") +
  theme(
    text = element_text(family = "Times New Roman"),
    # Align title to the full plot border instead of just the inner axes panel
    plot.title.position = "plot", 
    plot.title = element_text(family = "Times New Roman", face = "bold", hjust = 0.5, size = 13),
    axis.text.y = element_text(family = "Times New Roman", size = 11, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 10, color = "black"),
    legend.text = element_text(family = "Times New Roman"),
    legend.title = element_text(family = "Times New Roman"),
    legend.position = "top")
#+ Fisher tests in P4 code (additional, not included)

#EZHIP EXPRESSION COMPARISONS
Idents(tumor_cells) <- "status"
FindMarkers(tumor_cells, 
            ident.1 = "Primary", 
            ident.2 = "Metastasis",
            features = "EZHIP", 
            test.use = "wilcox")
FeaturePlot(
  tumor_cells, 
  features = "EZHIP", 
  order = TRUE) + 
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(family = "Times New Roman", face = "bold", hjust = 0.5),
    axis.title = element_text(family = "Times New Roman"),
    axis.text = element_text(family = "Times New Roman", color = "black"),
    legend.text = element_text(family = "Times New Roman"),
    legend.title = element_text(family = "Times New Roman")) +
  labs(title = "EZHIP Expression")


#EZHIP expression by cluster
tumor_cells$Tumor_State[tumor_cells$seurat_clusters == "0"] <- "Transcriptional maintenance MCC"
tumor_cells$Tumor_State[tumor_cells$seurat_clusters == "1"] <- "Differentiated MCC"
tumor_cells$Tumor_State[tumor_cells$seurat_clusters == "2"] <- "Stress response MCC"
tumor_cells$Tumor_State[tumor_cells$seurat_clusters == "3"] <- "Neuroendocrine signaling MCC"
tumor_cells$Tumor_State[tumor_cells$seurat_clusters == "4"] <- "Metabolic MCC"
tumor_cells$Tumor_State[tumor_cells$seurat_clusters == "5"] <- "Cycling S-phase MCC"
tumor_cells$Tumor_State[tumor_cells$seurat_clusters == "6"] <- "Secretory neuroendocrine MCC"
tumor_cells$Tumor_State[tumor_cells$seurat_clusters == "7"] <- "Cycling G2/M MCC"
tumor_cells$Tumor_State[tumor_cells$seurat_clusters == "9"] <- "Transitional MCC"
tumor_cells$Tumor_State[tumor_cells$seurat_clusters == "10"] <- "Metabolic stress MCC"

alphabetical_states <- sort(unique(na.omit(tumor_cells$Tumor_State)))
tumor_cells$Tumor_State <- factor(tumor_cells$Tumor_State, levels = alphabetical_states)
Idents(tumor_cells) <- "Tumor_State"
p <- VlnPlot(
  tumor_cells, 
  features = "EZHIP", 
  pt.size = 0)
max_val <- max(FetchData(tumor_cells, vars = "EZHIP")[, 1], na.rm = TRUE)

p + 
  geom_tile(
    data = p$data, 
    aes(x = ident, y = -0.05 * max_val, fill = ident), 
    height = 0.05 * max_val, 
    color = NA,     
    show.legend = FALSE) +
  coord_cartesian(ylim = c(-0.05 * max_val, max_val), clip = "off") +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(family = "Times New Roman", face = "bold", hjust = 0.5, size = 13),
    axis.title.y = element_text(family = "Times New Roman", size = 11, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 10, color = "black"),
    axis.title.x = element_text(family = "Times New Roman", size = 11, color = "black", vjust = -1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key = element_blank(),
    legend.text = element_text(family = "Times New Roman", size = 10),
    legend.title = element_text(family = "Times New Roman", size = 11) ) +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  labs(
    title = "EZHIP by Tumor Cell State",
    x = "Tumor State",
    fill = "Tumor State")


#EZHIP expression by subtype
df <- FetchData(tumor_cells, vars = c("EZHIP", "virus"))
pink_palette <- c("VN" = "#F4ACF9", "VP" = "#F4ACF9")

ggplot(df, aes(x = virus, y = EZHIP, color = virus)) +
  geom_jitter(
    position = position_jitter(width = 0.25), 
    alpha = 0.4, 
    size = 1.0) +
  scale_color_manual(values = pink_palette) +
  scale_x_discrete(labels = c("VN" = "Virus-negative", "VP" = "Virus-positive")) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  labs(
    title = "EZHIP Expression by Subtype",
    x = "Subtype",
    y = "EZHIP Expression") +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(family = "Times New Roman", hjust = 0.5, face = "bold", size = 14),
    legend.position = "none",
    axis.text.x = element_text(family = "Times New Roman", color = "black", size = 11),
    axis.text.y = element_text(family = "Times New Roman", color = "black"),
    axis.title.x = element_text(family = "Times New Roman", margin = margin(t = 8)),
    axis.title.y = element_text(family = "Times New Roman", margin = margin(r = 8)))

#PRC1/PRC2-TARGET GENE SIGNATURES
prc1_targets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CGP") %>%
  filter(
    grepl(
      "BMI1|RING1|RNF2|CBX7|CBX8|PCGF|PRC1",
      gs_name,
      ignore.case = TRUE)) %>%
  pull(gene_symbol) %>%
  unique()

prc2_targets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CGP") %>%
  filter(
    grepl(
      "EZH2|SUZ12|PRC2",
      gs_name,
      ignore.case = TRUE)) %>%
  pull(gene_symbol) %>%
  unique()

prc1_targets <- intersect(
  prc1_targets,
  rownames(tumor_cells))
prc2_targets <- intersect(
  prc2_targets,
  rownames(tumor_cells))

head(prc1_targets,20)
head(prc2_targets,20)

#CALCULATE SCORES (PRC1_UCell / PRC2_UCell)
gene_sets <- list(
  PRC1 = prc1_targets,
  PRC2 = prc2_targets)

tumor_cells <- AddModuleScore_UCell(
  tumor_cells,
  features = gene_sets,
  maxRank = 5000) #


#PRC1/2 BY SUBTYPE
#COMPARE VP/VN
VlnPlot(
  tumor_cells,
  features = c(
    "PRC1_UCell",
    "PRC2_UCell"),
  group.by = "virus",
  pt.size = 0,
  ncol = 2)

wilcox.test(
  PRC1_UCell ~ virus,
  data = FetchData(
    tumor_cells,
    vars = c(
      "virus",
      "PRC1_UCell")))
wilcox.test(
  PRC2_UCell ~ virus,
  data = FetchData(
    tumor_cells,
    vars = c(
      "virus",
      "PRC2_UCell")))

#CORRELATION EZHIP/PRC
ezhip <- FetchData(
  tumor_cells,
  vars = "EZHIP")[,1]
is.numeric(ezhip)

cor.test(
  ezhip,
  tumor_cells$PRC1_UCell,
  method = "spearman")

cor.test(
  ezhip,
  tumor_cells$PRC2_UCell,
  method = "spearman")

FeatureScatter(
  tumor_cells,
  feature1 = "EZHIP",
  feature2 = "PRC1_UCell")

FeatureScatter(
  tumor_cells,
  feature1 = "EZHIP",
  feature2 = "PRC2_UCell")

FeatureScatter(
  tumor_cells,
  feature1 = "PRC1_UCell",
  feature2 = "PRC2_UCell")

#CORRELATION BETWEEN PRC1/PRC2
cor.test(
  tumor_cells$PRC1_UCell,
  tumor_cells$PRC2_UCell,
  method = "spearman")


#DEFINE EZHIP-STATUS (POS/NEG)
tumor_cells$EZHIP_status <- ifelse(
  FetchData(tumor_cells, vars = "EZHIP")[,1] > 0,
  "EZHIP_positive",
  "EZHIP_negative")

df_long <- tumor_cells@meta.data %>%
  mutate(
    EZHIP_status = factor(
      EZHIP_status,
      levels = c("EZHIP_negative", "EZHIP_positive"),
      labels = c("EZHIP-negative", "EZHIP-positive"))) %>%
  pivot_longer(
    cols = c(PRC1_UCell, PRC2_UCell),
    names_to = "Signature",
    values_to = "Score")


#PRC1/2 BY EZHIP STATUS
ggplot(df_long, aes(x = EZHIP_status, y = Score, fill = EZHIP_status)) +
  geom_jitter(
    aes(color = EZHIP_status), 
    width = 0.15, 
    size = 0.2, 
    alpha = 0.15, 
    show.legend = FALSE) +
  geom_violin(
    trim = FALSE, 
    alpha = 0.7, 
    color = "black", 
    linewidth = 0.4,
    draw_quantiles = 0.5 ) +
  facet_wrap(~ Signature, scales = "free_y") +
  scale_fill_manual(
    name = "EZHIP status",
    values = c("EZHIP-negative" = "#E67768", "EZHIP-positive" = "#4DBBD5")) +
  scale_color_manual(
    values = c("EZHIP-negative" = "#C0392B", "EZHIP-positive" = "#16A085")) +
  labs(
    title = "PRC1/PRC2 Target Signature Score by EZHIP Status",
    x = "EZHIP status",
    y = "Signature Score") +
  theme_classic(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.title.x = element_text(face = "plain", margin = margin(t = 15)),
    axis.title.y = element_text(face = "plain"),
    axis.text = element_text(face = "plain"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "plain"))


#CORRELATION EZHIP/PRC
cor.test(
  ezhip,
  tumor_cells$PRC1_UCell,
  method = "spearman")

cor.test(
  ezhip,
  tumor_cells$PRC2_UCell,
  method = "spearman")

#EZHIP BY SAMPLE
#Single cell level
fisher.test(table(tumor_cells$virus, tumor_cells$EZHIP_status))

#tumor level
tumor_level <- tumor_cells@meta.data %>%
  dplyr::group_by(sample, virus) %>%
  dplyr::summarise(
    EZHIP_detected = any(EZHIP_status == "EZHIP_positive"),
    .groups = "drop")
table(tumor_level$virus, tumor_level$EZHIP_detected)
fisher.test(table(tumor_level$virus, tumor_level$EZHIP_detected))

