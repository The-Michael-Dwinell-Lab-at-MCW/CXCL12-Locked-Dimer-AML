- 👋 Donovan's standard workflow for scRNA-seq analysis
- 📫 How to reach me: ddrouillard@mcw.edu

```{r Install Packages}
# The first step is to install packages necessary for analysis. Here are the packages we typically use
install.packages('Seurat')
install.packages("tidyverse")
install.packages("dplyr")
install.packages("magrittr")
install.packages("data.table")
install.packages("Matrix")
install.packages("devtools")
install.packages("RcppArmadillo")
install.packages("Rcpp")
install.packages("scales")
install.packages("pheatmap")
install.packages("gplots")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("tibble")
install.packages("hdf5r")
install.packages("ggpubr")
install.packages("RColorBrewer")
install.packages("sctransform")
install.packages("MASS")
install_github("ctlab/fgsea")
install.packages("xtable")
BiocManager::install("topGO")
install.packages("DBI")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("org.Mm.eg.db")
# This is the mouse genome, which can be used in GO analyses
BiocManager::install("org.Hs.eg.db")
# This is the human genome, which can be used in GO analyses
BiocManager::install("DOSE")
BiocManager::install("MAST")
install.packages("devtools")
```

```{r Load Packages}
# Downloaded packages are not automatically available for use in R
# The library() function loads packages for use in one session of R
# If you restart R, you don't need to install the packages again, but you do need to library() them
library(Seurat)
library(tidyverse)
library(dplyr)
library(magrittr)
library(data.table)
library(Matrix)
library(devtools)
library(RcppArmadillo)
library(Rcpp)
library(scales)
library(pheatmap)
library(gplots)
library(ggplot2)
library(cowplot)
library(tibble)
library(data.table)
library(hdf5r)
library(ggpubr)
library(RColorBrewer)
library(sctransform)
library(xtable)
library(topGO)
library(installr)
library(MASS)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(DOSE)
library(MAST)
library(DBI)
# The following function updates R if necessary
if(!require(installr)) {
  install.packages("installr"); 
  require(installr)
} #load / install+load installr
updateR()
```

```{r Load data}
# Set your working directory (wd) to where your barcodes, matrix, and feature files are stored
setwd("C:/Users/location of file")

# Double check your working directory is where you want it to be
getwd()

# Here is where you name the samples as they will be referenced in R studio. Capitilzation matters!
V1 <- Read10X(data.dir = "Vehicle1")
V2 <- Read10X(data.dir = "Vehicle2")
T1 <- Read10X(data.dir = "Treated1")
T2 <- Read10X(data.dir = "Treated2")

# The samples are currently stored as a matrix object in R (matrix of genes and cells)
# We can check how many cells and genes are available for analysis with the dim() function
dim(V1)
dim(V2)
dim(T1)
dim(T2)

# We need to make a Seurat object from the matrix objects
V1 <- CreateSeuratObject(V1, project = "V1")
V2 <- CreateSeuratObject(V2, project = "V2")
T1 <- CreateSeuratObject(T1, project = "T1")
T2 <- CreateSeuratObject(T2, project = "T2")

#Add sample names to the barcodes so the samples don't get mixed up when they are integrated (combined into one Seurat object)
V1 <- RenameCells(V1, add.cell.id = "V1")
V2 <- RenameCells(V2, add.cell.id = "V2")
T1 <- RenameCells(T1, add.cell.id = "T1")
T2 <- RenameCells(T2, add.cell.id = "T2")

#Add metadata to each object specifying which treatment they received, so we can sort by treatment in the future analysis
V1$Treatment <- "Vehicle"
V2$Treatment <- "Vehicle"
T1$Treatment <- "Treated"
T2$Treatment <- "Treated"
```

```{r Quality Control}
# QC is done by both gene expression and % of mitochondrial genes
# Cells with too few genes are excluded (low read quality), too many genes (doublets), or high percentage of mitochondrial genes (apoptotic cells)

V1$percent.mt <- PercentageFeatureSet(V1, pattern = "^mt-")
V2$percent.mt <- PercentageFeatureSet(V2, pattern = "^mt-")
T1$percent.mt <- PercentageFeatureSet(T1, pattern = "^mt-")
T2$percent.mt <- PercentageFeatureSet(T2, pattern = "^mt-")

# PercentageFeatureSet calculates the frequency of mitochondrial genes. It adds that frequency on a per cell basis to the metadata in a new column called "percent.mt"
# The follow commmands allow us to visualize the number of genes, how many total genes, and percentage mitochondrial DNA

VlnPlot(V1, features =c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =3)
(which(V1$percent.mt <5) %>% length()) / length(V1$percent.mt) * 100
#76.12821% of cells in this sample have less than 5% mitochondrial genes
(which(V1$percent.mt <10) %>% length()) / length(V1$percent.mt) * 100
#87.68936% of cells in this sample have less than 10% mitochondrial genes

VlnPlot(V2, features =c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =3)
(which(V2$percent.mt <5) %>% length()) / length(V2$percent.mt) * 100
#86.39755% of cells in this sample have less than 5% mitochondrial genes
(which(V2$percent.mt <10) %>% length()) / length(V2$percent.mt) * 100
#94.07882% of cells in this sample have less than 10% mitochondrial genes

VlnPlot(T1, features =c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =3)
(which(T1$percent.mt <5) %>% length()) / length(T1$percent.mt) * 100
#84.62991% of cells in this sample have less than 5% mitochondrial genes
(which(T1$percent.mt <10) %>% length()) / length(T1$percent.mt) * 100
#93.92946% of cells in this sample have less than 10% mitochondrial genes

VlnPlot(T2, features =c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =3)
(which(T2$percent.mt <5) %>% length()) / length(T2$percent.mt) * 100
#82.46869% of cells in this sample have less than 5% mitochondrial genes
(which(T2$percent.mt <10) %>% length()) / length(T2$percent.mt) * 100
#91.34615% of cells in this sample have less than 10% mitochondrial genes

#I'll go with a 10% mito gene cutoff, so we get a larger percentage of V1. Up to 20% is the highest I've seen acceptable.
#I'll use a nFeature_RNA cutoff of 4000, to avoid inclusion of doublets but maintaining a large portion of each sample. The lowest cutoff I have seen is 3500

V1 <- subset(V1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
V2 <- subset(V2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
T1 <- subset(T1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
T2 <- subset(T2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

# We can check the number of cells avaiable for analysis after the QC
dim(V1)
dim(V2)
dim(T1)
dim(T2)

```{r scTransform}
#May need to increase size limits of variables R allows to be used for parallel processing (default is 500 MB. 1000MB = 1 GB). This will depend on the RAM in your computer
options(future.globals.maxSize = 4000 * 1024^2)

#scTransform replaces Normalize, FindVariableFeatures, and ScaleData functions.

V1 <- SCTransform (V1, vars.to.regress = "percent.mt", verbose = FALSE)
V2 <- SCTransform (V2, vars.to.regress = "percent.mt", verbose = FALSE)
T1 <- SCTransform (T1, vars.to.regress = "percent.mt", verbose = FALSE)
T2 <- SCTransform (T2, vars.to.regress = "percent.mt", verbose = FALSE)
# Can ignore all warnings that say "iteration limit reached"
# https://github.com/ChristophH/sctransform/issues/25
```


```{r Integrate Data}
#Combining multiple samples requires integration to account for batch effects. Do not use the merge function. Integration can be completed as follows, but also using other packages such as "Harmony"
Samples.list <- list(V1, V2, T1, T2)

Samples.features <- SelectIntegrationFeatures(object.list = Samples.list, nfeatures = 3500)
#You can use less features if your computer is struggling, try not to go below 2000

Samples_Reg.list <- PrepSCTIntegration(object.list = Samples.list, anchor.features = Samples.features, verbose = FALSE)

#Typical integration commands would look as follows (do not run until reading rest of integration chunk)
#Following command may take some time (~30 minutes) to run, but can run on 32GB RAM for this sample set which contains a considerable number of cells
Samples.anchors <- FindIntegrationAnchors(object.list = Samples_Reg.list, normalization.method = "SCT", anchor.features = Samples.features, verbose = FALSE)

#However, the dataset is very large and standard integration is too computationally intensive for a 32GB RAM computer
#Samples <- IntegrateData(anchorset = Samples.anchors, normalization.method = "SCT", verbose = FALSE)


#Instead, we can perform reference based integration that use 1-2 samples as a reference that all samples are compared to, rather than comparing all samples to each other
#This is less computationally intensive, but not as thorough as first method. We will use one vehicle and one treated as a reference
#Example and explanation can be found here: https://satijalab.org/seurat/articles/integration_large_datasets.html
Samples_ref.anchors <- FindIntegrationAnchors(object.list = Samples_Reg.list, reference = c(2, 3), normalization.method = "SCT", anchor.features = Samples.features, verbose = FALSE)

Samples <- IntegrateData(anchorset = Samples_ref.anchors, normalization.method = "SCT", verbose = FALSE)

#Still not running? Run the following three commands to check/help clear up some space
memory.size()
memory.limit()
gc()

#At this point, we can remove individual sample Seurat objects and integration intermediates to free up RAM (ONLY DO THIS AFTER INTEGRATION IS COMPLETE)
rm(V1, V2, T1, T2, Samples.features, Samples.list, Samples_Reg.list, Samples.anchors, Samples_ref.anchors)
rm(cluster0, T_cell_subset, T_cell_subset_markers)
```


``` {r Cell cycle sorting and factor levels}
DefaultAssay(Samples) <- "SCT"

# Load genes for cell cycle regression
CC_genes <- Seurat::cc.genes.updated.2019

# These genes are currently in human gene format (all caps)
# We need to change them to mouse gene format (first letter caps, all other letters lowercase)
CC_genes <- lapply(CC_genes, function(x){
  paste0(
  substring(x, 1, 1),
  tolower(substring(x, 2)))})

#Calculate cell cycle phases
Samples <- CellCycleScoring(Samples, s.features = CC_genes$s.genes, g2m.features = CC_genes$g2m.genes)

#Get cell cycle score based on treatment
table(Samples$Treatment)
```

``` {r Samples UMAP}
DefaultAssay(Samples) <- "integrated"
#Use the integrated data matrix ONLY to run functions involved in PCA/UMAP calculations
#Setting the default assay to integrated is not required for viewing UMAP plots

Samples <- RunPCA(Samples, npcs = 50) #50 PCs is default for scTransform
Samples <- RunUMAP(Samples, dims = 1:50)
Samples <- FindNeighbors(Samples, reduction = "pca", dims = 1:50)

#to determine appropriate resolution for clustering, you can use multiK. This is not necessary
#https://github.com/siyao-liu/MultiK
#install.packages("sigclust")
#install_github("siyao-liu/MultiK")
#library(sigclust)
#library(MultiK)


Samples <- FindClusters(Samples, resolution = 1.2)
#The more cells, the higher the resolution required. Do not overcluster, but do not undercluster. May take some guess/check

#ggsave is a command to save the last figure run in R. Save UMAPs as PNGs, and all other figures as PDFs

DimPlot(Samples, label = T)
ggsave("Images/01_UMAP.png", dpi = 300)

DimPlot(Samples, split.by = "Treatment", ncol = 2, label = T)
ggsave("Images/01_UMAP_Treatment.png", dpi = 300)

DimPlot(Samples, group.by = "Phase")
ggsave("Images/01_UMAP_Phase.png", dpi = 300)

```{r Normalize RNA assay after scTransform}
DefaultAssay(Samples) <- "RNA"
#Need to set the default assay to RNA for anything with real gene expression. Using the integrated assay will result in incorrect gene expression (i.e. negative values on violin plots)
```

```{r Identifying differentially expressed genes (DEGs)}
# If we don't know what genes are important, we can analyze differentially expressed genes (DEGs) in each cluster. Not required
Samples_markers <- FindAllMarkers(Samples, only.pos = T)
# Setting only.pos = T tells Seurat to only look for genes that are upregulated in each cluster
Samples_markers <- FindAllMarkers(Samples)
Samples_markers <- cbind(Row.names =rownames(Samples_markers), Samples_markers)
#The following command will give us the top 10 differentially expressed genes per cluster. This can be useful for identifying which cells are making up the clusters
#You can use more more or less than 10 genes if desired. Just change the "n=#" in the following command (I recommend changing the assigned name too to match) 
top_10_per_cluster <- Samples_markers%>%group_by(cluster)%>%top_n(n=10, wt = avg_log2FC)
Idents(Samples)
#The following will plot a heatmap of the differentially expressed genes. We only use 100 cells per cluster (downsample =#) instead of all the cells to make it run faster. Recommended if a large dataset
DoHeatmap(subset(Samples, downsample = 100), features = top_10_per_cluster$gene, size =2.5, assay = 'SCT', slot = 'scale.data')
ggsave("Images/01_Heatmap_cluster_genes.png", dpi = 300, width = 8, height = 12)

#Find DEGs by treatment
levels(Samples)
Idents(Samples) <- "Treatment"
Idents(Samples) <- 'seurat_clusters'
Treated.markers <- FindMarkers(Samples, ident.1 = "Treated", ident.2 = "Vehicle", logfc.threshold = log(2), min.pct = 0.25, only.pos = TRUE)
Treated.markers %>% group_by("Treatment") %>% top_n(n = 100, wt = avg_log2FC)
table(Treated.markers %>% group_by("Treatment") %>% top_n(n = 50, wt = avg_log2FC))
top100 <- Treated.markers %>% group_by("Treatment") %>% top_n(n = 100, wt = avg_log2FC)
top100pval <- subset(top100, rowSums(top100[5] <0.05) >0)
Treated.markers
```

````{r Clustering}
#Can use a dot plot to examine known cell-specific genes (i.e. CD3 = T cells)
#If you have a lot of genes on one DotPlot, you can rotate the axis labels to save space
DotPlot(Samples, features = c('Ms4a1', 'Cd79a', 'Jchain', 'Cd79b', 'Blnk', 'Slamf7', 'Fcrl5', 'Cd14', 'Cd68', 'Cd163', 'Csf1r', 'Fcgr3', 'Fcgr2b', 'Fcgr1', 'Msr1', 'Mrc1', 'Arg1' , 'Nos2', 'Cd40', 'Cd80', 'Cd83', 'Ccr7', 'Ccr6', 'Flt3', 'Itgax', 'Itgae', 'Cd86', 'Ncam1', 'Klrb1', 'Abcc4', 'Adcyap1', 'Cpa3', 'Ctsg', 'Mpo', 'Tpsb2', 'Cd3d', 'Cd3e', 'Cd8a', 'Cd44', 'Sell', 'Cxcr6', 'Cx3cr1', 'Foxp3', 'Il2ra', 'Gzma', 'Pdcd1', 'Lag3', 'Tigit', 'Ctla4', 'Cxcr3', 'Ccr5', 'Il17a', 'Ifng', 'Ccr10', 'Havcr2', 'Ly6g', 'Ly6c1', 'Xcr1', 'Anpep', 'Cd33', 'Clec12a', 'Bst2', 'Klra17', 'Ptprc', 'Siglech', 'Crem', 'Krt19', 'Vwf', 'Lum')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#Combine this information with the DEGs heatmap you created to assign cell identities
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)
new.cluster.ids <- c("PMN-MDSC/Neutrophil","PMN-MDSC/Neutrophil","PMN-MDSC/Neutrophil","M-MDSC/Macrophage","M-MDSC/Macrophage","M-MDSC/Macrophage", "Mast cell","PMN-MDSC/Neutrophil","Cd4+ T cell","Cd8+ T cell","M-MDSC/Macrophage","Epithelial",
                     "Intermediate MDSC","PMN-MDSC/Neutrophil", "DC","DC","B cell", "M-MDSC/Macrophage", "Cd4+ T cell", "DC", "PMN-MDSC/Neutrophil","DC","M-MDSC/Macrophage","NK cell","M-MDSC/Macrophage",
                     "M-MDSC/Macrophage",  "PMN-MDSC/Neutrophil", "pDC", "Eosinophil", "Basophil", "B cell", "Cd4+ T cell", 
                     "M-MDSC/Macrophage", "Activated Mast cell")
Samples <- RenameIdents(object = Samples, new.cluster.ids)
DimPlot(Samples, label = T, split.by = "Treatment")

# Let's examine the top 10 DEGs per cluster
Samples_markers %>% group_by("Vehicle", "Treatment") %>% top_n(10, wt = avg_log2FC)
```

``` {r Graphs}
#Dot Plot Poster Quality Example
DotPlot(Han_epithelial, features = c('GATA6', 'CCL20')) + 
  theme(text = element_text(face = "bold"),
        axis.text.x=element_text(size=20),
        axis.title = element_text(size=20,face="bold"),
        axis.title.y.right = element_text(size=20,face="bold"),
        axis.text.y=element_text(hjust=1, size=20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        axis.line = element_line(size=2))
ggsave("Images/CCL20_Classical.png", bg="white", dpi = 300)
```

```{r Subsetting populations}
#Example subset will be on CD8 T cells
#Tcell subsetting: https://www.nature.com/articles/s41467-019-12464-3

Cd8T_cell_subset <- subset(Samples, idents = "Cd8+ T cell")
Cd8T_cell_subset <- SCTransform (Cd8T_cell_subset, vars.to.regress = "percent.mt", verbose = FALSE)
Cd8T_cell_subset <- RunPCA(Cd8T_cell_subset, npcs = 50) #50 PCs is default for scTransform
Cd8T_cell_subset <- RunUMAP(Cd8T_cell_subset, dims = 1:50)
Cd8T_cell_subset <- FindNeighbors(Cd8T_cell_subset, reduction = "pca", dims = 1:50)
Cd8T_cell_subset <- FindClusters(Cd8T_cell_subset, resolution = 0.2)
#low resolution because low number of cells
DimPlot(Cd8T_cell_subset, split.by = "Treatment", ncol = 2, label = T)
#check genes of interest just for funsies
DotPlot(Cd8T_cell_subset, features =c('Cxcr3', 'Ccr5', 'Sell', 'Cd44', 'Pdcd1', 'Tigit', 'Ctla4'))
VlnPlot(Cd8T_cell_subset, features =c('Cxcr3', 'Ccr5', 'Sell', 'Cd44', 'Pdcd1', 'Tigit', 'Ctla4'), split.by = "Treatment")
ggsave("Images/01_CD8_subsetting_genes_Vln.png", dpi = 300)
Cd8T_markers <- FindAllMarkers(Cd8T_cell_subset)
top_10_Cd8T <- Cd8T_markers%>%group_by(cluster)%>%top_n(n=10, wt = avg_log2FC)
DoHeatmap(Cd8T_cell_subset, features = top_10_Cd8T$gene, size =5.5, assay = 'SCT', slot = 'scale.data')
ggsave("Images/01_CD8_Top10_heatmap.png", dpi = 300)
FeaturePlot(Cd8T_cell_subset, features = c('Cxcr3', 'Ccr5', 'Sell', 'Cd44', 'Pdcd1', 'Tigit', 'Ctla4'), split.by = "Treatment", label = T)
Cd8T_cell_subset[["Cd8_collapsed"]] <- Idents(object = Cd8T_cell_subset)
Idents(Cd8T_cell_subset) <- "seurat_clusters"
table(Cd8T_cell_subset$seurat_clusters, Cd8T_cell_subset$Cd8_collapsed)
?table
table(Cd8T_cell_subset$seurat_clusters)
DimPlot(Cd8T_cell_subset, split.by = "Treatment", label = T)
ggsave("Images/CD8_subset_by_treatment.png", dpi = 300)
md <- Cd8T_cell_subset@meta.data %>% as.data.table
md[, .N, by = c("Treatment", "seurat_clusters")]
md[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N")


#Can use this to generate KEGG Terms
PMN_MDSC_Neutrophil_Subset <- subset(Samples, idents = "PMN-MDSC/Neutrophil")
Idents(object = PMN_MDSC_Neutrophil_Subset) <- "Treatment"
levels(PMN_MDSC_Neutrophil_Subset)
TreatmentDElist_PMN <-FindMarkers(object = PMN_MDSC_Neutrophil_Subset, test.use = "MAST", only.pos = T, ident.1 = "Vehicle", ident.2 = "Treated")
TreatmentDElist_PMN <- cbind(Row.Names =rownames(TreatmentDElist_PMN), TreatmentDElist_PMN)
PMN_gene_list <- TreatmentDElist_PMN$Row.Names
head(PMN_gene_list)
PMN_gene_list
TreatmentDElist_PMN_gene_convert <- bitr(PMN_gene_list, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Mm.eg.db", drop = F)
TreatmentDElist_PMN_gene_convert$SYMBOL[is.na(TreatmentDElist_PMN_gene_convert$ENTREZID)]
z <- enrichKEGG(gene=TreatmentDElist_PMN_gene_convert$ENTREZID, organism = 'mmu', pvalueCutoff = .05)
z <- setReadable(z, OrgDb = "org.Mm.eg.db", keyType = 'ENTREZID')
barplot(z, showCategory = 20)
#Need to check if this is genes upregulated in treated or vehicle, depends on which group is assigned ident.1 vs ident.2 in FindMarkers function

<!---
pdanowhere/pdanowhere is a ✨ special ✨ repository because its `README.md` (this file) appears on your GitHub profile.
You can click the Preview link to take a look at your changes.
--->
