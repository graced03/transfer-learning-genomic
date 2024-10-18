library(pgenlibr)
library(data.table)
library(dplyr)
library(Matrix)

# Define the function with an additional argument for chromosome numbers
process_gwas_results <- function(data_path, population_name, chromosomes = 1:22, P_value_thred = 0.0005) {
  gwas_results <- list()
  snp_count <- 0

  for (chr in chromosomes) {
    # Read the GWAS SNP table for the specified chromosome
    gwas_snp_table <- fread(paste0(data_path, population_name, "_snp_assoc/plink_gwas_ldl_", population_name, "_chr",
                                   as.character(chr), ".assoc.LDL.glm.linear"))

    # Convert P values to numeric and filter SNPs based on threshold
    gwas_snp_table$P <- as.numeric(gwas_snp_table$P)
    gwas_snp_table_sub <- gwas_snp_table %>% filter(P < P_value_thred)

    # Store the filtered SNP IDs and update SNP count
    gwas_results[[paste0("chr", as.character(chr))]] <- gwas_snp_table_sub$ID
    snp_count <- snp_count + length(gwas_snp_table_sub$ID)
  }

  return(list(results = gwas_results, snp_count = snp_count))
}

process_clumped_snps <- function(data_path, population_name, chromosomes = 1:22) {
  chr_snp_list <- list()
  snp_count <- 0

  for (chr in chromosomes) {
    # Read the clumped SNP table for the specified chromosome
    clump_snp_table <- fread(paste0(data_path, population_name, "_snp_clump/ldl_", population_name, "_chr",
                                    as.character(chr), "_plink_clumped.clumps"))

    # Store the SNP IDs and update SNP count
    chr_snp_list[[paste0("chr", as.character(chr))]] <- clump_snp_table$ID
    snp_count <- snp_count + length(clump_snp_table$ID)
  }

  return(list(chr_snp_list = chr_snp_list, snp_count = snp_count))
}

read_genotype_data <- function(phenotype, geno_dir, chr_snp_list, chromosomes=1:22) {
  # Read phenotype file
  # Loop through each chromosome in chr_list
  for (i in chromosomes) {
    chr <- paste0("chr", as.character(i))
    # Build paths for genotype files
    psam_path <- paste0(geno_dir, chr, ".psam")
    pgen_path <- paste0(geno_dir, chr, ".pgen")
    pvar_path <- paste0(geno_dir, chr, ".pvar")

    # Read .psam and .pvar files
    sample_data <- fread(psam_path)
    pvar_table <- fread(pvar_path)
    pvar_table$Index <- 1:nrow(pvar_table)

    # Initialize PGEN/PVAR objects
    pvar <- NewPvar(pvar_path)
    pgen <- NewPgen(pgen_path, pvar=pvar)

    # Select SNPs based on chr_snp_list
    snp_selected_chr <- pvar_table %>% filter(ID %in% chr_snp_list[[chr]])

    # Read the genotype matrix
    geno_mat <- ReadList(pgen, snp_selected_chr$Index, meanimpute=FALSE)
    colnames(geno_mat) <- chr_snp_list[[chr]]

    # Merge with sample data and phenotype
    merged <- cbind(sample_data[, 2], geno_mat)
    phenotype <- merge(phenotype, merged, by = "IID", all.x = TRUE)
  }
  # Return the final merged phenotype data
  return(phenotype)
}

# Example usage with chromosome control:
# Process only chromosomes 1, 5, and 10 for both populations
# irish_results <- process_gwas_results("./data/LDL/GWAS/output/hapmap3_snp/", "irish", chromosomes = c(1, 5, 10))
# british_results <- process_gwas_results("./data/LDL/GWAS/output/hapmap3_snp/", "british", chromosomes = c(1, 5, 10))

# Access the results
# irish_results$snp_count
# british_results$snp_count

# Example usage:
# chr_list <- 1:22  # List of chromosomes
# chr_snp_list <- list(...)  # A list where each chromosome has a vector of SNP IDs
# phenotype_file <- "~/UKBB/data/LDL/phenotype/dat_LDL_British_Irish.rds"
# geno_dir <- "./data/geno/"

# processed_data <- process_genotype_data(phenotype_file, geno_dir, chr_list, chr_snp_list)

