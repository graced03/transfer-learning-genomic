#!/bin/bash
#
#SBATCH --job-name=gwas-clumps-chr1-euro-restricted
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=5G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xding20@jh.edu

plink2 \
  --bfile /users/xding/UKBB/data/geno/chr22 \
  --covar /users/xding/UKBB/data/LDL/GWAS/input/cov_LDL_Irish.txt \
  --geno 0.05 \
  --hwe 1e-08 \
  --mach-r2-filter 0.8 2.0 \
  --maf 0.01 \
  --mind 0.05 \
  --glm hide-covar --covar-variance-standardize \
  --out /users/xding/UKBB/data/LDL/GWAS/output/euro_restricted/plink_gwas_ldl_irish_chr22.assoc \
  --pheno /users/xding/UKBB/data/LDL/GWAS/input/pheno_LDL_Irish.txt \
  --snps-only \
  --threads 24


