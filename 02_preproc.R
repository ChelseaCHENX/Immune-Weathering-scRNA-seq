# ------------------------------------------------------------
# Purpose:
#   Convert a Scanpy AnnData object (.h5ad) to a Seurat object (.rds),
#   then run Azimuth reference mapping (pbmcref) to annotate PBMC cell types
#   (and optionally ADT features, if present / do.adt=TRUE).
#
# Inputs:
#   Command-line argument:
#     sp  = sample identifier (used to locate input/output files)
#   File:
#     ../data/scanpy_tmp/{sp}.h5ad    (AnnData object produced in Python/Scanpy)
#
# Outputs:
#   Intermediate:
#     ../data/seurat_tmp/{sp}.rds     (Seurat object converted from AnnData)
#   Final:
#     ../data/02_preproc/seurat_tmp/{sp}.rds
#       - Seurat object after Azimuth mapping to "pbmcref"
#       - Adds predicted annotations (e.g., predicted cell types) and
#         embedding / mapping info returned by RunAzimuth
#
# Notes / Dependencies:
#   - functions.R: project helper functions (loaded first)
#   - sceasy + reticulate are used to convert AnnData -> Seurat
#   - Azimuth is used for reference-based annotation via RunAzimuth()
# ------------------------------------------------------------

source('functions.R')
#use_condaenv('arm_py38')
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

