# Core Seurat and related visualization packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(gridExtra)

# Data manipulation and formatting
library(dplyr)
library(tidyr)
library(gtools)
library(glue)

# Interoperability
library(anndata)

# ---- Optional packages (commented out for now) ----
# library(SeuratData)
# library(SeuratWrappers)
# library(SeuratDisk)
# library(Azimuth)

# library(sceasy)
# library(reticulate)
# library(scCustomize)
# InstallData("pbmcsca")

parse_gmt <- function(gmt_path) {
  db <- readLines(gmt_path)
  geneset <- list()

  for (line in db) {
    # split strictly by TAB
    fields <- strsplit(line, "\t")[[1]]

    set_name <- fields[1]
    genes    <- fields[-c(1, 2)]  # drop name + description

    geneset[[set_name]] <- genes
  }

  return(geneset)
}
