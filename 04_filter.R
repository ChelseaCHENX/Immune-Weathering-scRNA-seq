# ------------------------------------------------------------
# Purpose:
#   Per-sample filtering and merging of Seurat objects after Azimuth annotation.
#   For each dataset/sample:
#     1) Load Seurat object annotated by Azimuth (predicted cell types + scores)
#     2) Add a dataset-specific cell ID prefix to keep cell barcodes unique
#     3) Apply QC and annotation-confidence filters:
#          - nFeature_RNA > 1000 (retain higher-quality cells)
#          - predicted.celltype.l1.score >= az_thres (retain confident Azimuth labels)
#          - predicted.celltype.l1 != "other" (remove ambiguous lineage calls)
#     4) Save the filtered per-sample object
#     5) Incrementally merge all filtered samples into a single merged object
#   Finally, save the merged Seurat object for downstream integration analyses.
#
# Inputs:
#   data/01_meta/meta_pt.csv
#     - contains dataset_name column listing sample IDs to process
#   data/02_preproc/seurat_tmp/{ds}.rds
#     - per-sample Seurat objects after conversion + Azimuth mapping
#     - expected metadata fields:
#         predicted.celltype.l1
#         predicted.celltype.l1.score
#         nFeature_RNA
#
# Outputs:
#   data/02_preproc/seurat_filtered/{ds}.rds
#     - filtered Seurat object per sample
#   data/02_preproc/seurat_filtered/merge.qs
#     - merged Seurat object across all samples after filtering
#
# Key parameters:
#   az_thres = 0.9
#     - minimum Azimuth prediction score for retaining cells
#
# Notes:
#   - RenameCells(add.cell.id = ds) ensures unique cell names across samples,
#     which prevents barcode collisions during merge().
#   - Uses qs::qsave for faster serialization of large Seurat objects.
#   - rm() + gc() called each iteration to reduce memory footprint.
# ------------------------------------------------------------

library(Seurat)
source('functions.R')
library(future)
options(future.globals.maxSize = 180 * 1024 * 1024 * 1024)
setwd('../')

meta = read.table('data/01_meta/meta_pt.csv', sep=',', header=T) %>% arrange(dataset_name)


az_thres = .9

ct = 0
for (ds in meta$dataset_name){

    ct = ct + 1

    obj = readRDS(glue('data/02_preproc/seurat_tmp/{ds}.rds'))
    obj <- RenameCells(obj, add.cell.id = ds)
    obj@meta.data$dataset_name <- ds
    
    obj = subset(x = obj, subset = nFeature_RNA > 1000) # refine min genes to 1000

    obj = subset(x = obj, 
             subset = predicted.celltype.l1.score >= az_thres)
    obj = subset(x=obj, 
             subset = predicted.celltype.l1 != 'other')

    qs::qsave(obj, 
              file = glue("data/02_preproc/seurat_filtered/{ds}.rds")) 
    
    if (ct == 1) {mg = obj} else {mg = merge(mg, obj)}
    print(glue('{ct}, {ds}, {dim(obj)[1]}, {dim(mg)[1]}'))
    
    rm(obj); gc()  # Free memory after saving
    # if (ct > 2) break
}

qs::qsave(mg, 
          file = glue("data/02_preproc/seurat_filtered/merge.qs")) 