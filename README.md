# 🔬 Single-Cell Bioinformatics Project

This repository provides a comprehensive workflow for analyzing **single-cell transcriptomics data** using **R** and various bioinformatics tools. The pipeline includes **data preprocessing, cell-cell communication analysis with CellChat, spatial data visualization with Seurat, and deconvolution analysis with SCDC**.

---

## 🛠 System Setup

### **1. Download and Extract Data**

Download the dataset using `curl` and unzip it:

```bash
curl -O https://icbb-share.s3.eu-central-1.amazonaws.com/single-cell-bioinformatics/scbi_p3.zip
unzip scbi_p3.zip
```

### **2. Clone the Repository and Move Files**

Clone the repository to obtain the necessary R scripts:

```bash
git clone https://github.com/jissksaji/Single-Cell-Bioinformatics.git
```

Move the extracted data into the cloned repository:

```bash
mv project_3_dataset/* Single-Cell-Bioinformatics/
cd Single-Cell-Bioinformatics
```

### **3. Set Up the R Environment**

Ensure your R environment is ready and install the required packages:

```r
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
devtools::install_github("sqjin/CellChat")
```

### **4. Additional Installations for Windows Users**

Windows users may need to install the following packages manually:

- **SCDC**: [github.com/meichendong/SCDC](https://github.com/meichendong/SCDC)
- **glmgampoi**: [github.com/const-ae/glmGamPoi](https://github.com/const-ae/glmGamPoi)
- **biocneighbors**: [bioconductor.org](https://bioconductor.org/packages/release/bioc/html/BiocNeighbors.html)

---

## 📚 Required R Libraries

Before running the analysis, load the necessary libraries:

```r
library(Seurat)        # Spatial data analysis and visualization
library(dplyr)         # Data manipulation
library(ggplot2)       # Plotting
library(patchwork)     # Arranging multiple plots
library(SCDC)          # Deconvolution analysis
library(CellChat)      # Cell-cell communication analysis
```

---

## 📊 Overview of the Analysis Pipeline

1. **Preprocessing**: Load and normalize single-cell RNA-seq data.
2. **Quality Control**: Filter low-quality cells and genes to maintain high data integrity.
3. **Dimensionality Reduction**:
   - Perform **Principal Component Analysis (PCA)** and **Uniform Manifold Approximation and Projection (UMAP)** for clustering.
4. **Cell-Cell Communication Analysis**:
   - Utilize **CellChat** to infer intercellular communication networks and signaling pathways.
5. **Spatial Transcriptomics Visualization**:
   - Integrate with **Seurat** to map spatial data onto tissue slides.
6. **Deconvolution Analysis**:
   - Apply **SCDC** to estimate cell type proportions within bulk RNA-seq data.
7. **Differential Expression Analysis (DEA)**:
   - Identify key marker genes across clusters and spatial regions.
8. **Data Integration and Batch Correction**:
   - Merge multiple datasets while addressing batch effects using Seurat’s **FindIntegrationAnchors** method.

---

## 🌐 References

- [Seurat Documentation](https://satijalab.org/seurat/)
- [CellChat GitHub](https://github.com/sqjin/CellChat)
- [SCDC GitHub](https://github.com/meichendong/SCDC)
- [Bioconductor](https://bioconductor.org/)
- [10X Genomics Spatial Transcriptomics](https://www.10xgenomics.com/solutions/spatial-gene-expression)

---

This project provides a **comprehensive pipeline** for **analyzing single-cell transcriptomics data** with a focus on **spatial data visualization, intercellular communication, and deconvolution techniques**. The full report on this project is also available in this repository. 🚀

