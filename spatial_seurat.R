# Load required libraries
library(Seurat)        # For spatial data analysis and visualization
library(dplyr)         # For data manipulation
library(ggplot2)       # For plotting
library(patchwork)     # For arranging multiple plots
library(SCDC)          # For deconvolution analysis
library(CellChat)      # For cell-cell communication analysis

###### Spatial Transcriptomics Data Loading and Initial Processing (3P) ########

# Set working directory and seed for reproducibility
setwd("~/Library/CloudStorage/OneDrive-UniversitätdesSaarlandes/Desktop/singlecellbioinformatics/project_3_dataset")
set.seed(10)

# Load Section 1 spatial data
# - data.dir: Directory containing spatial data files
# - filename: H5 file with feature-barcode matrix
# - slice/image: Metadata identifiers for spatial coordinates
spatial_data1 <- Load10X_Spatial(
  data.dir = "Section_1/",
  filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5",
  slice = "slice1",
  image = Read10X_Image("Section_1/spatial/")  # Load spatial coordinates
)

# Load Section 2 spatial data with same parameters
spatial_data2 <- Load10X_Spatial(
  data.dir = "Section_2/",
  filename = "V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5",
  slice = "slice2",
  image = Read10X_Image("Section_2/spatial/")
)

# Calculate mitochondrial gene percentage for QC (common QC metric)
spatial_data1[["percent.mt"]] <- PercentageFeatureSet(spatial_data1, pattern = "^MT-")
spatial_data2[["percent.mt"]] <- PercentageFeatureSet(spatial_data2, pattern = "^MT-")

# Visualize QC metrics and save plots
VlnPlot(spatial_data1, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
ggsave("QC_Violin_Section1.png", width = 10, height = 6, dpi = 300)
VlnPlot(spatial_data2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
ggsave("QC_Violin_Section2.png", width = 10, height = 6, dpi = 300)

######### Data Preprocessing (4P) ##########

# Filter cells based on QC metrics (remove outliers)
spatial_data1 <- subset(
  spatial_data1,
  subset = nFeature_Spatial > 2000 & nFeature_Spatial < 7500 &
    nCount_Spatial > 10000 & nCount_Spatial < 45000
)

spatial_data2 <- subset(
  spatial_data2,
  subset = nFeature_Spatial > 2000 & nFeature_Spatial < 7500 &
    nCount_Spatial > 7500 & nCount_Spatial < 40000
)

# Re-check QC metrics after filtering
VlnPlot(spatial_data1, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
ggsave("QC_Violin_Section1_Filtered.png", width = 10, height = 6, dpi = 300)
VlnPlot(spatial_data2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
ggsave("QC_Violin_Section2_Filtered.png", width = 10, height = 6, dpi = 300)

# Inspect Seurat object structure
spatial_data1  # Shows object summary (counts, assays, metadata)
spatial_data2

# Access gene expression data (stored in 'Spatial' assay)
spatial_data1@assays$Spatial  # Raw counts
spatial_data2@assays$Spatial

# Access spatial coordinate images
spatial_data1@images  # Contains tissue position information
spatial_data2@images

# Visualize spatial expression of random genes
genes <- rownames(spatial_data1)
random_genes <- sample(genes, 2)
SpatialFeaturePlot(spatial_data1, features = random_genes, slot = "counts", pt.size.factor = 1, alpha = c(0.1, 1))
ggsave("Random_Genes_Section1.png", width = 8, height = 6, dpi = 300)

random_genes2 <- sample(genes, 2)
SpatialFeaturePlot(spatial_data2, features = random_genes2, slot = "counts", pt.size.factor = 1, alpha = c(0.1, 1))
ggsave("Random_Genes_Section2.png", width = 8, height = 6, dpi = 300)

# Normalize data using SCTransform (regularized negative binomial regression)
spatial_data1 <- SCTransform(spatial_data1, assay = "Spatial", verbose = FALSE)
spatial_data2 <- SCTransform(spatial_data2, assay = "Spatial", verbose = FALSE)

####### Dimensionality Reduction, Clustering, and Visualization (3P) #########

# Perform PCA for dimensionality reduction
spatial_data1 <- RunPCA(spatial_data1, verbose = FALSE)
spatial_data2 <- RunPCA(spatial_data2, verbose = FALSE)

# Determine significant PCs using elbow plot
ElbowPlot(spatial_data1, ndims = 30)
ggsave("Elbow_Section1.png", width = 6, height = 4, dpi = 300)
ElbowPlot(spatial_data2, ndims = 30)
ggsave("Elbow_Section2.png", width = 6, height = 4, dpi = 300)

# Run UMAP for 2D visualization
spatial_data1 <- RunUMAP(spatial_data1, dims = 1:20, verbose = FALSE)
spatial_data2 <- RunUMAP(spatial_data2, dims = 1:20, verbose = FALSE)

# Cluster cells using Louvain algorithm
spatial_data1 <- FindNeighbors(spatial_data1, dims = 1:20, verbose = FALSE)
spatial_data1 <- FindClusters(spatial_data1, verbose = FALSE)

spatial_data2 <- FindNeighbors(spatial_data2, dims = 1:20, verbose = FALSE)
spatial_data2 <- FindClusters(spatial_data2, verbose = FALSE)

# Visualize clusters in UMAP space
DimPlot(spatial_data1, reduction = "umap", group.by = "seurat_clusters") + ggtitle("Section1 Clusters")
ggsave("UMAP_Clusters_Section1.png", width = 8, height = 6, dpi = 300)
DimPlot(spatial_data2, reduction = "umap", group.by = "seurat_clusters") + ggtitle("Section2 Clusters")
ggsave("UMAP_Clusters_Section2.png", width = 8, height = 6, dpi = 300)

# Visualize clusters on spatial coordinates
SpatialDimPlot(spatial_data1, group.by = "seurat_clusters") + ggtitle("Section1 Spatial Clusters")
ggsave("Spatial_Clusters_Section1.png", width = 8, height = 6, dpi = 300)
SpatialDimPlot(spatial_data2, group.by = "seurat_clusters") + ggtitle("Section2 Spatial Clusters")
ggsave("Spatial_Clusters_Section2.png", width = 8, height = 6, dpi = 300)

####### Differential Expression Analysis (8P) ########

# Identify cluster markers using Wilcoxon test
clusters1 <- FindMarkers(
  spatial_data1, 
  ident.1 = 0,        # Compare cluster 0 against all others
  ident.2 = NULL, 
  test.use = "wilcox"
)  
head(clusters1)

clusters2 <- FindMarkers(
  spatial_data2, 
  ident.1 = 0, 
  ident.2 = NULL, 
  test.use = "wilcox"
)  
head(clusters2)

# Find positive markers for all clusters
allclusters1 <- FindAllMarkers(spatial_data1, test.use = "wilcox", only.pos = TRUE)  
allclusters2 <- FindAllMarkers(spatial_data2, test.use = "wilcox", only.pos = TRUE)  

# Visualize top marker genes spatially
gene_names <- head(rownames(spatial_data1), 3)
SpatialFeaturePlot(spatial_data1, features = gene_names, ncol = 3)
ggsave("Top_Markers_Section1.png", width = 10, height = 4, dpi = 300)

######## Merging and Integration (7P) ########

# Merge datasets from both sections
merged_data <- merge(spatial_data1, spatial_data2)

# Reprocess merged data
merged_data <- SCTransform(merged_data, assay = "Spatial", verbose = FALSE)
merged_data <- RunPCA(merged_data, verbose = FALSE)
ElbowPlot(merged_data, ndims = 20)
ggsave("Elbow_Merged.png", width = 6, height = 4, dpi = 300)

merged_data <- FindNeighbors(merged_data, dims = 1:15)
merged_data <- FindClusters(merged_data, resolution = 0.5)
merged_data <- RunUMAP(merged_data, dims = 1:15)

# Visualize merged clusters
DimPlot(merged_data, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave("UMAP_Merged.png", width = 8, height = 6, dpi = 300)
SpatialDimPlot(merged_data, group.by = "seurat_clusters", label.size = 5)
ggsave("Spatial_Merged.png", width = 8, height = 6, dpi = 300)

# Check batch effects between original sections
DimPlot(merged_data, reduction = "umap", group.by = "orig.ident")
ggsave("UMAP_BatchCheck.png", width = 8, height = 6, dpi = 300)

# Advanced integration using Seurat's integration workflow
merge.list <- list(spatial_data1, spatial_data2)
merge.list <- lapply(merge.list, SCTransform, assay = "Spatial", method = "poisson")
features <- SelectIntegrationFeatures(object.list = merge.list)
anchors <- FindIntegrationAnchors(object.list = merge.list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)

# Process integrated data
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, dims = 1:30)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

# Save processed objects
saveRDS(spatial_data1, "spatialdata1.rds")
saveRDS(spatial_data2, "spatialdata2.rds")
saveRDS(combined, "combined.rds")

###### Cell-type Identification (8P) ########

# Load reference scRNA-seq data (Allen Brain Atlas)
allen <- readRDS("allen_cortex.rds")

# Transfer cell type labels from reference to spatial data
anchors <- FindTransferAnchors(reference = allen, query = combined, dims = 1:30)
combined <- MapQuery(anchors, combined, reference = allen, refdata = allen$subclass)

# Visualize transferred labels
DimPlot(combined, reduction = "umap", group.by = "predicted.id", label = TRUE) +
  ggtitle("Transferred Cell Types")
ggsave("UMAP_CellTypes.png", width = 8, height = 6, dpi = 300)

SpatialDimPlot(combined, group.by = "predicted.id") +
  ggtitle("Spatial Cell Types")
ggsave("Spatial_CellTypes.png", width = 8, height = 6, dpi = 300)

# Manual annotation using known marker genes
marker_genes <- c("Gfap", "Atp13a4", "Cd68", "Arg2", "Apold1", "A2m")
FeaturePlot(combined, features = marker_genes, cols = c("grey", "red"))
ggsave("Marker_UMAPs.png", width = 10, height = 8, dpi = 300)

######## Deconvolution and Cell-Cell Communication (7P) ########

# Deconvolution with SCDC (requires reference data)
# Note: The following is a conceptual example; actual implementation may vary
# deconv_results <- SCDC_deconv(spatial_data, reference = allen)

# CellChat analysis for signaling networks
cellchat <- createCellChat(
  object = combined,
  meta = combined@meta.data,
  group.by = "predicted.id",
  datatype = "spatial",
  coordinates = GetTissueCoordinates(combined)
)

# Configure database and compute communication probabilities
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- computeCommunProb(cellchat)

# Visualize communication networks
netVisual_circle(cellchat@net$count, title.name = "Interaction Counts")
ggsave("CellChat_Counts.png", width = 8, height = 6, dpi = 300)
netVisual_circle(cellchat@net$weight, title.name = "Interaction Weights")
ggsave("CellChat_Weights.png", width = 8, height = 6, dpi = 300)

# Pathway-specific visualization
netVisual_aggregate(cellchat, signaling = "PSAP", layout = "spatial")
ggsave("PSAP_Signaling.png", width = 8, height = 6, dpi = 300)