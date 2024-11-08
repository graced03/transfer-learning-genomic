#!/bin/bash
#
#SBATCH --job-name=extract-hm-snps-and-gwas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=10G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xding20@jh.edu


for chr in {1..22}
do
  echo "Processing chromosome ${chr}..."
  plink2 --pfile /dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr${chr} \
    --extract /users/xding/UKBB/data/hapmap3/euro_rsid_chr${chr}.txt \
    --snps-only \
    --max-alleles 2 \
    --make-pgen \
    --geno 0.05 \
    --hwe 1e-08 \
    --mach-r2-filter 0.8 2.0 \
    --maf 0.01 \
    --mind 0.05 \
    --out /users/xding/UKBB/data/geno/chr${chr} \
    --threads 24

  plink2 \
    --pfile /users/xding/UKBB/data/geno/chr${chr} \
    --covar /users/xding/UKBB/data/LDL/GWAS/input/cov_LDL_Irish.txt \
    --glm hide-covar --covar-variance-standardize \
    --out /users/xding/UKBB/data/LDL/GWAS/output/hapmap3_snp/plink_gwas_ldl_irish_chr${chr}.assoc \
    --pheno /users/xding/UKBB/data/LDL/GWAS/input/pheno_LDL_Irish.txt \
    --no-input-missing-phenotype \
    --snps-only \
    --threads 24
done
