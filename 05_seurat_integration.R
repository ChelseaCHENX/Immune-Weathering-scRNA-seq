# ------------------------------------------------------------
# Purpose:
#   Perform per-lineage integration and clustering of merged PBMC scRNA-seq data.
#   The merged Seurat object contains multiple datasets (dataset_name). For each
#   major cell lineage (and for all cells), we:
#     1) subset to that lineage (based on Azimuth predicted.celltype.l1)
#     2) remove datasets with too few cells (<= 1 cell) to avoid unstable integration
#     3) split the RNA assay by dataset and run standard preprocessing
#     4) integrate across datasets using Harmony on PCA embeddings
#     5) run graph-based clustering and UMAP using the Harmony reduction
#     6) save an integrated Seurat object per lineage for downstream analyses
#
# Inputs:
#   data/02_preproc/seurat_filtered/merge.qs
#     - merged Seurat object with metadata:
#         predicted.celltype.l1 : major lineage labels (e.g., B, CD4 T, Mono, etc.)
#         dataset_name          : dataset / batch identifier
#
# Outputs:
#   data/03_cleandata/itg_{celltype}.qs
#     - integrated and clustered Seurat object for each entry in `celltypes`
#       (including "all")
#
# Notes:
#   - Integration is run within each lineage to reduce confounding from major
#     lineage differences and to improve alignment across batches.
#   - HarmonyIntegration is applied with dims = 10 to compute the "harmony"
#     reduction; clustering and UMAP use dims = 1:30 from the harmony space.
#   - options(future.globals.maxSize) is increased to support large Seurat objects.
#   - In practice, this loop was executed as separate jobs (e.g., via SLURM),
#     one celltype per run, to manage memory/runtime.
# ------------------------------------------------------------

# Header
setwd('/n/data1/dfci/surgonc/mittendorf/lab/chelsea/iw/github')
source('codes/functions.R')
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(glue)

options(future.globals.maxSize = 200 * 1024 * 1024 * 1024)

# iterate (in practice, it was by submitting individual scripts via slurm)
celltypes = c('all', 'B','CD4 T','CD8 T','other T','DC','Mono','NK')

dat = qs::qread("data/02_preproc/seurat_filtered/merge.qs")

for (celltype in celltypes){

    if (celltype == 'all'){obj = dat} else{
    obj = subset(x = dat, 
                 subset = predicted.celltype.l1 == celltype)        
    }

    cell_counts <- table(obj$dataset_name)
    valid_datasets <- names(cell_counts[cell_counts > 1])   # Keep only if >1 cell
    obj <- subset(obj, subset = dataset_name %in% valid_datasets)
    
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$dataset_name)
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj)
    
    obj <- IntegrateLayers(
      object = obj, method = HarmonyIntegration,
      orig.reduction = "pca", new.reduction = "harmony",
      reference = 1, dims=10, verbose = T
    )
    
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]]) 
    
    obj <- FindNeighbors(obj, dims = 1:30, reduction = "harmony")
    obj <- FindClusters(obj, resolution = .1)
    obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony")
    
    qs::qsave(obj, file = glue("data/03_cleandata/itg_{celltype}.qs")) 
}