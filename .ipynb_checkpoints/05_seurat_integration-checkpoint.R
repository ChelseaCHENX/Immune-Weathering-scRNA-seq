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