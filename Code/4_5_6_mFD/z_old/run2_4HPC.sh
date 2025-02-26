#!/bin/bash
#SBATCH --job-name=llew0024_4_Gawdis
#SBATCH --mail-user=john.llewelyn@flinders.edu.au
#SBATCH --mail-type=ALL
#SBATCH --output=/home/llew0024/%x-%j.out.txt
#SBATCH --error=/home/llew0024/%x-%j.err.txt
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs
#SBATCH --time=0-12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=600G

# Load R module
export PATH=/home/llew0024/R/4.4.0/bin:$PATH

# Run R script
Rscript /home/llew0024/4HPC_gawdisANDmFD_wBetaCUT.R

