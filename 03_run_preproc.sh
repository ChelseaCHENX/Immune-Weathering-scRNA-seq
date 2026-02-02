#!/bin/bash

# cd "/Users/chelseachen/Partners HealthCare Dropbox/Chelsea Chen/projects/iw/scRNA/scripts" ## have to run in this folder d/t blanks in file name 

python_path='/n/data1/dfci/surgonc/mittendorf/lab/chelsea/packages/miniconda3/envs/scrna_new/bin/python'#'/home/fac135/miniconda3/envs/scrna/bin/python'

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