library(data.table)
library(readr)

# data <- fread("./data/LDL/GWAS/output/hapmap3_snp/irish_snp_assoc/plink_gwas_ldl_irish_chr1.assoc.LDL.glm.linear", 
#               sep = "\t", na.strings = "")
# View(data %>% filter(P <= 0.0005))

#################################### GWAS ######################################
irish_gwas_results <- process_gwas_results("./data/LDL/GWAS/output/hapmap3_snp/", "irish")
british_gwas_results <- process_gwas_results("./data/LDL/GWAS/output/hapmap3_snp/", "british")
irish_gwas_results$snp_count
british_gwas_results$snp_count

################################## Clumped #####################################
british_clumped_snps <- process_clumped_snps("./data/LDL/clump/", "british", chromosomes = 1:22)
british_irish_clumped_snps <- process_clumped_snps("./data/LDL/clump/", "british_irish", chromosomes = 1:22)
british_irish_clumped_snps$snp_count
british_clumped_snps$snp_count

################################ Inspection ####################################

irish_chr_snps_all_sig_0005 <- Reduce(c, irish_gwas_results$results)
british_chr_snps_all_sig_0005 <- Reduce(c, british_gwas_results$results)

length(intersect(irish_chr_snps_all_sig_0005, british_chr_snps_all_sig_0005))
# 300

chr_snp_list <- irish_gwas_results$results
phenotype <- readRDS("~/UKBB/data/LDL/phenotype/dat_LDL_British_Irish.rds")

geno_dir <- "./data/geno/"
phenotype_irish <- phenotype %>% filter(ethnicity == "Irish")
extracted_geno_data_irish_gwas <- read_genotype_data(phenotype_irish, geno_dir, irish_gwas_results$results)
dim(extracted_geno_data_irish_gwas)
saveRDS(phenotype, file = "/fastscratch/myscratch/xding/UKBB/irish_gwas/dat_LDL_Irish_pheotype_genotype_irish_gwas_all_snps_p0005.rds")


phenotype_british <- phenotype %>% filter(ethnicity == "British")

extracted_geno_data_british_gwas <- read_genotype_data(phenotype_british, geno_dir, british_gwas_results$results)
dim(extracted_geno_data_british_gwas)
saveRDS(extracted_geno_data_irish_gwas, file = "/fastscratch/myscratch/xding/UKBB/british_gwas/dat_LDL_British_pheotype_genotype_british_gwas_all_snps_p0005.rds")


#################################### GWAS P-value < 0.005 ######################################
irish_gwas_results <- process_gwas_results("./data/LDL/GWAS/output/hapmap3_snp/", "irish", P_value_thred = 0.005)
irish_gwas_results$snp_count

chr_snp_list <- irish_gwas_results$results
phenotype <- readRDS("~/UKBB/data/LDL/phenotype/dat_LDL_British_Irish.rds")

geno_dir <- "./data/geno/"
extracted_geno_data_irish_gwas <- read_genotype_data(phenotype, geno_dir, irish_gwas_results$results)
dim(extracted_geno_data_irish_gwas)
saveRDS(phenotype, file = "/fastscratch/myscratch/xding/UKBB/irish_gwas/dat_LDL_British_Irish_pheotype_genotype_irish_gwas_all_snps_p005.rds")

