#!/bin/bash
#
#SBATCH --job-name=select-unrelated-irish-id-diff-cutoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xding20@jh.edu


mkdir /fastscratch/myscratch/xding/UKBB/irish_genotype/hapmap3
for chr in {1..22}
do
  echo "Processing chromosome ${chr}..."
  plink2 --pfile /users/xding/UKBB/data/geno/chr${chr} \
     --keep /users/xding/UKBB/data/ethnicity_group_ids/irish.txt \
     --make-pgen \
     --out /fastscratch/myscratch/xding/UKBB/irish_genotype/hapmap3/chr${chr}
done

for chr in {1..22}; do
  echo "/fastscratch/myscratch/xding/UKBB/irish_genotype/hapmap3/chr$chr" >> /users/xding/UKBB/data/ethnicity_group_ids/irish_hapmap3_pfile_list.txt
done

plink2 \
  --pmerge-list /users/xding/UKBB/data/ethnicity_group_ids/irish_hapmap3_pfile_list.txt \
  --make-pgen \
  --out /fastscratch/myscratch/xding/UKBB/irish_genotype/all

plink2 \
 --pfile /fastscratch/myscratch/xding/UKBB/irish_genotype/hapmap3/all \
 --king-cutoff 0.05 \
 --out /fastscratch/myscratch/xding/UKBB/irish_genotype/hapmap3/related_pair_cutoff05

plink2 \
 --pfile /fastscratch/myscratch/xding/UKBB/irish_genotype/hapmap3/all \
 --king-cutoff 0.05 \
 --make-rel \
 --out /fastscratch/myscratch/xding/UKBB/irish_genotype/hapmap3/unrelated_cutoff05
