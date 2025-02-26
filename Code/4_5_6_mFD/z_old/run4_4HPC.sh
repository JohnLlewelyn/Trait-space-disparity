#!/bin/bash
#SBATCH --job-name=llew0024_4_Gawdis
#SBATCH --mail-user=john.llewelyn@flinders.edu.au
#SBATCH --mail-type=ALL
#SBATCH --output=/home/llew0024/%x-%j.out.txt
#SBATCH --error=/home/llew0024/%x-%j.err.txt
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=600G
module load R
Rscript /home/llew0024/4HPC_BetaPara.R
