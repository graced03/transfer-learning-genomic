#!/bin/bash
#
#SBATCH --job-name=hm-snps-clump
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --mem-per-cpu=15G
#SBATCH --time=05:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xding20@jh.edu


for chr in {1..22}
do
  echo "Processing chromosome ${chr}..."
  # plink2 \
  #   --pfile /users/xding/UKBB/data/geno/chr${chr} \
  #   --covar /users/xding/UKBB/data/LDL/GWAS/input/cov_LDL_British_Irish.txt \
  #   --glm hide-covar --covar-variance-standardize \
  #   --out /users/xding/UKBB/data/LDL/GWAS/output/euro_restricted/british_irish_snp_assoc/plink_gwas_ldl_british_irish_chr${chr}.assoc \
  #   --pheno /users/xding/UKBB/data/LDL/GWAS/input/pheno_LDL_British_Irish.txt \
  #   --no-input-missing-phenotype \
  #   --snps-only \
  #   --threads 24

  plink2 --pfile /users/xding/UKBB/data/geno/chr${chr} \
  --pheno /users/xding/UKBB/data/LDL/GWAS/input/pheno_LDL_British_Irish.txt \
  --covar /users/xding/UKBB/data/LDL/GWAS/input/cov_LDL_British_Irish.txt \
  --no-input-missing-phenotype \
  --clump /users/xding/UKBB/data/LDL/GWAS/output/euro_restricted/british_irish_snp_assoc/plink_gwas_ldl_british_irish_chr${chr}.assoc.LDL.glm.linear \
  --out /users/xding//UKBB/data/LDL/clump/british_irish_snp_clump/ldl_british_irish_chr${chr}_plink_clumped \
  --clump-kb 2000 \
  --clump-p1 0.0005 \
  --clump-r2 0.6 \
  --threads 22
done
