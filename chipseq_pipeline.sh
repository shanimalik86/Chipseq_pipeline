#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1  
#SBATCH -p normal
#SBATCH --mem=8gb
module load snakemake
snakemake -s Snakefile --cluster-config cluster.yaml --cluster 'sbatch --mem={cluster.mem} -c {cluster.cpus}' -j 16
