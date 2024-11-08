# GWAS and basic thresholding with plink and R on JHPCE

## Step 1: Extract the hapmap3 snps and save the genotype data in the local directory

1.a. Save the snps to .txt files
```
R/extract_hapmap3_snps.R
```

1.b. Extract the snps
```
sbatch extract_hm_snps.sh
```

## Step 2: Construct the covariate and phenotype data tables for GWAS

```
R/create_cov_and_phenotype_datatable.Rmd
```

### Step 2.1: Identify and remove related individuals

```
sbatch identify_and_remove_related_irish_individuals.sh
```

### Step 2.2: construct the covariate and phenotype data with unrelated individuals

```
R/create_cov_and_phenotype_datatable.Rmd
```

## Step 3: GWAS and Clumping with plink2

```
sbatch irish_gwas_and_clump.sh
```

## Step 5: modeling

```

R/baseline_regression_un_related_irish.Rmd
```


Totorial
https://choishingwan.github.io/PRS-Tutorial/plink/
