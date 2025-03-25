
##############################################
# Load Required Libraries
##############################################
library(dplyr)        # Data manipulation
library(Seurat)       # Single-cell RNA-seq analysis
library(patchwork)    # Combining ggplot2 plots
library(ggplot2)      # Plotting
library(cowplot)      # Additional plotting functions
library(RColorBrewer) # Color palettes
library(pheatmap)     # Heatmap visualization
library(readr)        # Data reading
library(tidyr)        # Data tidying
library(viridis)      # Color scales

##############################################
# Define Directories Containing 10X Data
##############################################
# Find folders ending with "hss" and "hcs" in the current directory
hss_folders <- list.files('./', pattern = 'hss$')
hsc_folders <- list.files('./', pattern = 'hcs$')

##############################################
# Create Seurat Objects for HSS Data
##############################################
# Load 10X data from each folder and create Seurat objects.
# Only cells with at least 60 counts and 200 features are kept.
hssList <- lapply(hss_folders, function(folder) { 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder, 
                     min.cells = 60, 
                     min.features = 200)
})

# Manually assign group labels based on file order.
# (Adjust indices if the folder order changes.)
hssList[[1]]@meta.data$group <- "control"
hssList[[2]]@meta.data$group <- "3-del"
hssList[[3]]@meta.data$group <- "3-del"
hssList[[4]]@meta.data$group <- "5-del"
hssList[[5]]@meta.data$group <- "5-del"

##############################################
# Create Seurat Objects for HCS Data
##############################################
hcsList <- lapply(hsc_folders, function(folder) { 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder, 
                     min.cells = 60, 
                     min.features = 200)
})

# Manually assign group labels.
hcsList[[1]]@meta.data$group <- "control"
hcsList[[2]]@meta.data$group <- "3-del"
hcsList[[3]]@meta.data$group <- "3-del"
hcsList[[4]]@meta.data$group <- "5-del"
hcsList[[5]]@meta.data$group <- "5-del"

##############################################
# Merge Seurat Objects for HSS and HCS
##############################################
# Merge individual objects into a single Seurat object.
# Note: For hss, check if the correct folder names (hss_folders) should be used 
# in add.cell.ids instead of hsc_folders.
hss <- merge(hssList[[1]], 
             y = c(hssList[[2]], hssList[[3]], hssList[[4]], hssList[[5]]),
             add.cell.ids = hsc_folders, 
             project = "human")

hcs <- merge(hcsList[[1]], 
             y = c(hcsList[[2]], hcsList[[3]], hcsList[[4]], hcsList[[5]]),
             add.cell.ids = hsc_folders, 
             project = "human")

##############################################
# Calculate Mitochondrial and MK Gene Percentages
##############################################
# For quality control, calculate the percentage of mitochondrial (MT) and MK (marker) genes.
hss[["percent.mt"]] <- PercentageFeatureSet(hss, pattern = "^MT-")
hss[["percent.mk"]] <- PercentageFeatureSet(hss, pattern = "^MK-")

hcs[["percent.mt"]] <- PercentageFeatureSet(hcs, pattern = "^MT-")
hcs[["percent.mk"]] <- PercentageFeatureSet(hcs, pattern = "^MK-")

##############################################
# Filter and Normalize Data
##############################################
# Filter out cells based on gene counts and mitochondrial/MK percentages.
hss <- subset(hss, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &
                     nCount_RNA > 2000 & percent.mt < 30 & percent.mk > 1)
# Normalize using log-normalization.
hss <- NormalizeData(hss, normalization.method = "LogNormalize", scale.factor = 10000)

hcs <- subset(hcs, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &
                     nCount_RNA > 2000 & percent.mt < 30 & percent.mk > 1)
hcs <- NormalizeData(hcs, normalization.method = "LogNormalize", scale.factor = 10000)

##############################################
# Scale Data, Run PCA, Clustering and Dimensionality Reduction (HSS)
##############################################
# Scale data for hss object
hss <- ScaleData(hss)
# Select MK genes for PCA (using genes that match the pattern "MK-")
mk.genes <- rownames(hss)[grep("MK-", rownames(hss))]
hss <- RunPCA(object = hss, features = mk.genes)
# Find nearest neighbors and clusters using the PCA results (dims 1:15)
hss <- FindNeighbors(hss, reduction = "pca", dims = 1:15)
hss <- FindClusters(hss, resolution = 0.2)
# Run t-SNE and UMAP for visualization
hss <- RunTSNE(hss, dims = 1:15, do.fast = TRUE)
DimPlot(hss, reduction = "tsne", label = TRUE, label.size = 6)
hss <- RunUMAP(hss, dims = 1:15)
DimPlot(hss, reduction = "umap", label = TRUE, label.size = 6)

# Identify marker genes for each cluster and export results
hss.markers <- FindAllMarkers(object = hss, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# Select top 50 markers per cluster based on average log fold-change and export
hss.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> hss.top50

##############################################
# Scale Data, Run PCA, Clustering and Dimensionality Reduction (HCS)
##############################################
hcs <- ScaleData(hcs)
mk.genes <- rownames(hcs)[grep("MK-", rownames(hcs))]
hcs <- RunPCA(object = hcs, features = mk.genes)
# t-SNE and clustering for hcs object
hcs <- FindNeighbors(hcs, reduction = "pca", dims = 1:15)
hcs <- FindClusters(hcs, resolution = 0.2)
hcs <- RunTSNE(hcs, dims = 1:15, do.fast = TRUE)
DimPlot(hcs, reduction = "tsne", label = TRUE, label.size = 6)
hcs <- RunUMAP(hcs, dims = 1:15)
DimPlot(hcs, reduction = "umap", label = TRUE, label.size = 6)

# Identify marker genes for hcs clusters and export results
hcs.markers <- FindAllMarkers(object = hcs, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
hcs.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> hcs.top50

##############################################
# Rename Cluster Identities to Cell Types
##############################################
# For the hss object, assign cell type labels to clusters based on prior knowledge.
hss <- RenameIdents(hss, 
                    `0` = "GABAergic Neurons", 
                    `1` = "Glutamatergic Neurons", 
                    `2` = "Glutamatergic Neurons", 
                    `3` = "Astroglia", 
                    `4` = "Radial Glia", 
                    `5` = "OPC", 
                    `6` = "GABAergic Neurons", 
                    `7` = "Intermediate Progenitors")
DimPlot(hss, reduction = "tsne", label = TRUE)
DimPlot(hss, reduction = "umap", label = TRUE, label.size = 6)

# For the hcs object, assign cell type labels to clusters.
hcs <- RenameIdents(hcs, 
                    `0` = "Glutamatergic Neurons", 
                    `1` = "Glutamatergic Neurons", 
                    `2` = "Astroglia", 
                    `3` = "Radial Glia", 
                    `4` = "Choroid Plexus", 
                    `5` = "GABAergic Neurons", 
                    `6` = "GABAergic Neurons", 
                    `7` = "Radial Glia", 
                    `8` = "Intermediate Progenitors")
DimPlot(hcs, reduction = "tsne", label = TRUE)
DimPlot(hcs, reduction = "umap", label = TRUE, label.size = 6)

##############################################
# Custom Color Palette and UMAP Plotting
##############################################
celltype.colors <- c("GABAergic Neurons" = "#F8766D", 
                     "Glutamatergic Neurons" = "#B79F00", 
                     "Astroglia" = "#00BA38", 
                     "Radial Glia" = "#00BFC4",
                     "Intermediate Progenitors" = "#F564E3", 
                     "Choroid Plexus" = "#00B9E3", 
                     "unknown" = "white")

# Plot UMAP with custom cell type colors
DimPlot(hcs, group.by = "ident", cols = celltype.colors)
DimPlot(hss, group.by = "ident", cols = celltype.colors)

##############################################
# Create Bar Plots for Cell Type Proportions
##############################################
# For hss: Create a frequency table and plot proportions by sample group.
cell_counts_hss <- table(hss@meta.data$group, hss@active.ident)
cell_counts_hss <- as.data.frame(cell_counts_hss)
cell_counts_hss$Var1 <- factor(cell_counts_hss$Var1, levels = c("control", "5-del", "3-del"))
ggplot(cell_counts_hss, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.75) +
  xlab("Sample") +
  ylab("Proportion") +
  theme_classic() +
  theme(legend.title = element_blank())

# For hcs: Do the same as above.
cell_counts_hcs <- table(hcs@meta.data$group, hcs@active.ident)
cell_counts_hcs <- as.data.frame(cell_counts_hcs)
cell_counts_hcs$Var1 <- factor(cell_counts_hcs$Var1, levels = c("control", "5-del", "3-del"))
ggplot(cell_counts_hcs, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.75) +
  xlab("Sample") +
  ylab("Proportion") +
  theme_classic() +
  theme(legend.title = element_blank())

##############################################
# Feature and Dot Plots for Selected Genes
##############################################
# Define genes of interest for feature visualization.
features <- c("FOXG1", "NRXN1", "EMX1", "DLX2")
FeaturePlot(hss, features = features)
FeaturePlot(hcs, features = features)

# Create dot plots for a set of merged markers.
merged_marker <- c("TUBB", "FOXG1", "NRXN1", "TBR1", "GRIA1", "GRIA2", "DLX1", "GAD1", "GAD2", "S100B")
DotPlot(hss, cols = c("lightgrey", "blue", "red"), assay = "RNA", 
        features = merged_marker, scale = FALSE, split.by = "group")
DotPlot(hcs, cols = c("lightgrey", "blue", "red"), assay = "RNA", 
        features = merged_marker, scale = FALSE, split.by = "group")

##############################################
# Differential Expression Analysis
##############################################
# For GABAergic neurons in hss: Subset cells and identify differentially expressed genes.
GABA.cells <- subset(hss, idents = "GABAergic Neurons")
# Note: The following two calls to FindMarkers will overwrite each other.
# Consider storing them in different variables if both results are needed.
GABA.DEGs <- FindMarkers(GABA.cells, ident.1 = "2607hss", ident.2 = c("581hss", "641hss"), 
                         group.by = 'orig.ident', verbose = FALSE)
GABA.DEGs <- FindMarkers(GABA.cells, ident.1 = "2607hss", ident.2 = c("972hss", "973hss"), 
                         group.by = 'orig.ident', verbose = FALSE)

# For Glutamatergic neurons in hcs: Subset cells and identify differentially expressed genes.
Glu.cells <- subset(hcs, idents = "Glutamatergic Neurons")
Glu.DEGs <- FindMarkers(Glu.cells, ident.1 = "2607hcs", ident.2 = c("581hcs", "641hcs"), 
                        group.by = 'orig.ident', verbose = FALSE)
Glu.DEGs <- FindMarkers(Glu.cells, ident.1 = "2607hcs", ident.2 = c("972hcs", "973hcs"), 
                        group.by = 'orig.ident', verbose = FALSE)

##############################################
# Heatmap Visualization for Fold-Change Data
##############################################
# Define a  magma color palette from the viridis package.
color_palette <- viridis::magma(50)

# Plot heatmap for down-regulated pathways in GABAergic cells (HSS)
hss_down.foldchange <- read.table("GABA.hss_down.foldchange.txt", head = FALSE, sep = "\t")
colnames(hss_down.foldchange) <- c("pathway", "GABA_5del", "GABA_3del")
rownames(hss_down.foldchange) <- hss_down.foldchange[,1]
hss_down.foldchange <- hss_down.foldchange[,-(1:3)]
pheatmap(hss_down.foldchange, cluster_col = FALSE, cluster_row = FALSE, 
         color = color_palette, border_color = "white", breaks = seq(1, 6, length.out = 51))

# Plot heatmap for up-regulated pathways in GABAergic cells (HSS)
hss_up.foldchange <- read.table("GABA.hss_up.foldchange.txt", head = FALSE, sep = "\t")
colnames(hss_up.foldchange) <- c("pathway", "GABA_5del", "GABA_3del")
rownames(hss_up.foldchange) <- hss_up.foldchange[,1]
hss_up.foldchange <- hss_up.foldchange[,-1]
pheatmap(hss_up.foldchange, cluster_col = FALSE, cluster_row = FALSE, 
         color = color_palette, border_color = "white", breaks = seq(1, 6, length.out = 51))

# Plot heatmap for up-regulated pathways in Glutamatergic cells (HCS)
hcs_up.foldchange <- read.table("Glu.hcs_up.foldchange.txt", head = FALSE, sep = "\t")
colnames(hcs_up.foldchange) <- c("pathway", "Glu_5del", "Glu_3del")
rownames(hcs_up.foldchange) <- hcs_up.foldchange[,1]
hcs_up.foldchange <- hcs_up.foldchange[,-1]
pheatmap(hcs_up.foldchange, cluster_col = FALSE, cluster_row = FALSE, 
         color = color_palette, border_color = "white", breaks = seq(1, 6, length.out = 51))

# Plot heatmap for down-regulated pathways in Glutamatergic cells (HCS)
hcs_down.foldchange <- read.table("Glu.hcs_down.foldchange.txt", head = FALSE, sep = "\t")
colnames(hcs_down.foldchange) <- c("pathway", "Glu_5del", "Glu_3del")
rownames(hcs_down.foldchange) <- hcs_down.foldchange[,1]
hcs_down.foldchange <- hcs_down.foldchange[,-1]
pheatmap(hcs_down.foldchange, cluster_col = FALSE, cluster_row = FALSE, 
         color = color_palette, border_color = "white", breaks = seq(1, 6, length.out = 51))
