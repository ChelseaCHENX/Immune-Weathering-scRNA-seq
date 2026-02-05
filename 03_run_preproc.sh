#!/bin/bash
# ------------------------------------------------------------
# Purpose:
#   Batch preprocessing pipeline for PBMC scRNA-seq samples.
#   For each sample:
#     1) Run Python-based Scanpy preprocessing (QC, filtering, doublet removal)
#     2) Run R-based Seurat preprocessing and Azimuth annotation
#
# Workflow:
#   - Iterates over Seurat .rds files in ../data/seurat_tmp/
#   - Extracts sample ID (sp_name) from filename
#   - Skips samples that already have a processed Scanpy output (.h5ad)
#   - Runs:
#       a) 01_preproc.py   (Python / Scanpy preprocessing)
#       b) 02_preproc.R    (R / Seurat + Azimuth annotation)
#   - Logs stdout/stderr and timestamps for each step
#
# Inputs:
#   ../data/seurat_tmp/{sp_name}.rds   (used only to define sample list)
#
# Outputs:
#   Python:
#     ../data/scanpy/{sp_name}.h5ad
#   R:
#     ../data/02_preproc/seurat_tmp/{sp_name}.rds
#
# Logging:
#   ../logs/101824_{sp_name}.txt
#   Includes timestamps, progress messages, and error traces.
#
# Error handling:
#   - If Python preprocessing fails → stop pipeline
#   - If R preprocessing fails → stop pipeline
#   - Ensures downstream steps only run after upstream success
#
# Environment:
#   - Uses a fixed conda Python environment for reproducibility
#   - Assumes Rscript and required R packages are available in PATH
#
# Usage:
#   bash 03_run_preproc.sh
# ------------------------------------------------------------

python_path='/n/data1/dfci/surgonc/mittendorf/lab/chelsea/packages/miniconda3/envs/scrna_new/bin/python'

cd /n/data1/dfci/surgonc/mittendorf/lab/chelsea/iw/codes

# Loop through each directory in the specified folder
for sp in ../data/seurat_tmp/*; do

  # Check if it's a directory (ignores files)

    basename="${sp##*/}"      
    sp_name="${basename%.rds}"

    if [ ! -f "../data/scanpy/${sp_name}.h5ad" ]; then

        log_file="../logs/101824_$sp_name.txt"
        # Log the start time
        echo "[$(date)] Starting process for $sp_name" > $log_file

        echo "Processing $sp_name" >> $log_file

        # Run the Python script with the current 'sp' as an argument
        $python_path 01_preproc.py --sp "$sp_name" >> $log_file 2>&1

        # Wait for the Python script to complete before proceeding to R script
        if [ $? -eq 0 ]; then
        # Log time after Python finishes
        
        echo "[$(date)] Python script completed for $sp_name" >> $log_file

        echo "Python script completed successfully for $sp_name, now running R script" >> $log_file

        # Run the R script with the same 'sp' argument
        Rscript 02_preproc.R "$sp_name" >> $log_file 2>&1

        # Check if the R script was successful
        if [ $? -eq 0 ]; then
            # Log time after R script finishes
            echo "[$(date)] R script completed for $sp_name" >> $log_file
            echo "R script completed successfully for $sp_name" >> $log_file
        else
            echo "Error: R script failed for $sp_name" >> $log_file
            exit 1  # Exit if the R script fails
        fi
        else
        echo "Error: Python script failed for $sp_name" >> $log_file
        exit 1  # Exit if the Python script fails
        fi
    fi
    # Break after the first iteration is done
    # break

done

# bash 03_run_preproc.sh 