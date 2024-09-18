
plink2 \
  --pfile /dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr1 \
  --covar /users/xding/UKBB/data/LDL/GWAS/input/cov_LDL_Irish.txt \
  --geno 0.05 \
  --hwe 1e-08 \
  --mach-r2-filter 0.8 2.0 \
  --maf 0.01 \
  --mind 0.05 \
  --glm hide-covar
  --out /users/xding/UKBB/data/LDL/GWAS/output/plink_gwas_ldl_irish_chr1.assoc \
  --pheno /users/xding/UKBB/data/LDL/GWAS/input/pheno_LDL_Irish.txt \
  --snps-only \
  --threads 24

  
plink2 --pfile /dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr1 \
  --rm-dup force-first \
  --snps-only \
  --max-alleles 2 \
  --make-pgen \
  --pheno /users/xding/UKBB/data/LDL/GWAS/input/pheno_LDL_Irish.txt \
  --out /users/xding/UKBB/data/geno/

plink2 --pfile /dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr1 \
  --extract /users/xding/UKBB/data/hapmap3/euro_rsid.txt \
  --rm-dup force-first \
  --snps-only \
  --max-alleles 2 \
  --make-pgen \
  --pheno /users/xding/UKBB/data/LDL/GWAS/input/pheno_LDL_Irish.txt \
  --out /users/xding/UKBB/data/geno/
  
plink2 --pfile /dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr22 \
  --extract /users/xding/UKBB/data/hapmap3/euro_rsid_chr22.txt \
  --rm-dup force-first \
  --snps-only \
  --max-alleles 2 \
  --make-bed \
  --out /users/xding/UKBB/data/geno/chr22 \
  --threads 24