#!/bin/bash
#File:
#SBATCH -t 0-4:00:00                       # <-- walltime in format of days-hours:minutes:sec
#SBATCH --account=def-gbader                  # <-- computecanada account
#SBATCH --mem-per-cpu=120G                   # <-- memory
#SBATCH -J scRNA_outfile                 # <-- name of job
#SBATCH -c 2                          # <-- number of CPUs
#SBATCH -o %x-%j.out                      # <-- redirect job output (both stdout and stderr)
#SBATCH --mail-user=charlotte.volk@mail.mcgill.ca      # <-- email address to UHN always
#SBATCH --mail-type=ALL                     # <-- both start and end


module load StdEnv/2020  
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14 # R is a dependency
module load geos/3.10.2
module load dplyr/1.0.7
module load Seurat/4.0.5
module load ggplot2/3.3.5
module load gridExtra/2.3
module load ggpubr/0.4.0
module load scran/1.22.1

Rscript ~/projects/def-gbader/cvolk/scripts/scvi.R
#Rscript ~/test.R
