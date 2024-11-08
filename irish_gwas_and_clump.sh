#!/bin/bash
#
#SBATCH --job-name=gwas-trait-explore-1031
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=15G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xding20@jh.edu

# LDL
base_dir="/users/xding/UKBB/data"

for chr in {1..22}
do
  echo "Processing chromosome ${chr}..."
  plink2 \
    --pfile ${base_dir}/geno/chr${chr} \
    --covar ${base_dir}/LDL/GWAS/input/cov_LDL_unrelated_Irish.txt \
    --glm hide-covar --covar-variance-standardize \
    --out ${base_dir}/LDL/GWAS/output/hapmap3_snp/unrelated_irish_snp_assoc/plink_gwas_LDL_unrelated_irish_chr${chr}.assoc \
    --pheno ${base_dir}/LDL/GWAS/input/pheno_LDL_unrelated_Irish.txt \
    --no-input-missing-phenotype \
    --snps-only \
    --threads 24
done

output_dir="${base_dir}/LDL/clump/unrelated_irish_snp_assoc"
mkdir -p "${output_dir}"

# Loop over each chromosome
for chr in {1..22}; do
    echo "Processing chromosome ${chr}..."
    # Run clumping command for each chromosome and trait
    plink2 --pfile ${base_dir}/geno/chr${chr} \
           --pheno ${base_dir}/LDL/GWAS/input/pheno_LDL_unrelated_Irish.txt \
           --covar ${base_dir}/LDL/GWAS/input/cov_LDL_unrelated_Irish.txt \
           --no-input-missing-phenotype \
           --clump ${base_dir}/LDL/GWAS/output/hapmap3_snp/unrelated_irish_snp_assoc/plink_gwas_LDL_unrelated_irish_chr${chr}.assoc.LDL.glm.linear \
           --out ${output_dir}/LDL_irish_chr${chr}_plink_clumped \
           --clump-kb 250 \
           --clump-p1 0.005 \
           --clump-r2 0.1 \
           --threads 24
done
