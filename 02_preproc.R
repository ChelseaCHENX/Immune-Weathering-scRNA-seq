source('functions.R')
#use_condaenv('arm_py38') # run on my Mac M1
library(sceasy)
library(reticulate)
library(glue)
library(Azimuth)

args <- commandArgs(trailingOnly = TRUE)
sp <- args[1]

sceasy::convertFormat(glue('../data/scanpy_tmp/{sp}.h5ad'), from="anndata", to="seurat",
                       outFile=glue('../data/seurat_tmp/{sp}.rds'))

srt = readRDS(glue('../data/seurat_tmp/{sp}.rds'))
srt <- RunAzimuth(srt, reference = "pbmcref", do.adt=TRUE)


saveRDS(srt, file = glue('../data/02_preproc/seurat_tmp/{sp}.rds'))

