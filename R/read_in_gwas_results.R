library(data.table)
library(readr)

# data <- fread("./data/LDL/GWAS/output/hapmap3_snp/irish_snp_assoc/plink_gwas_ldl_irish_chr1.assoc.LDL.glm.linear", 
#               sep = "\t", na.strings = "")
# View(data %>% filter(P <= 0.0005))

gwas_results <- list()
snp_count <- 0
for (chr in 1:22){
  gwas_snp_table <- fread(paste0("./data/LDL/GWAS/output/hapmap3_snp/irish_snp_assoc/plink_gwas_ldl_irish_chr",
                                  as.character(chr),
                                  ".assoc.LDL.glm.linear"))
  gwas_snp_table$P <- as.numeric(gwas_snp_table$P)
  gwas_snp_table_sub <- gwas_snp_table %>% filter(P < 0.0005)
  gwas_results[[paste0("chr",as.character(chr))]] <- gwas_snp_table_sub$ID
  snp_count <- snp_count + length(gwas_snp_table_sub$ID)
}
snp_count

british_gwas_results <- list()
british_snp_count <- 0
for (chr in 1:22){
  british_gwas_snp_table <- fread(paste0("./data/LDL/GWAS/output/hapmap3_snp/british_snp_assoc/plink_gwas_ldl_british_chr",
                                 as.character(chr),
                                 ".assoc.LDL.glm.linear"))
  british_gwas_snp_table$P <- as.numeric(british_gwas_snp_table$P)
  british_gwas_snp_table_sub <- british_gwas_snp_table %>% filter(P < 0.0005)
  british_gwas_results[[paste0("chr",as.character(chr))]] <- british_gwas_snp_table_sub$ID
  british_snp_count <- british_snp_count + length(british_gwas_snp_table_sub$ID)
}
british_snp_count

################################## British #####################################

british_chr_snp_list <- list()
british_snp_count <- 0
for (chr in 1:22){
  clump_snp_table <- fread(paste0("./data/LDL/clump/british_snp_clump/ldl_british_chr",
                                  as.character(chr),
                                  "_plink_clumped.clumps"))
  british_chr_snp_list[[paste0("chr",as.character(chr))]] <- clump_snp_table$ID
  british_snp_count <- british_snp_count + length(clump_snp_table$ID)
}

############################### British Irish ##################################

british_irish_chr_snp_list <- list()
british_irish_snp_count <- 0
for (chr in 1:22){
  clump_snp_table <- fread(paste0("./data/LDL/clump/british_irish_snp_clump/ldl_british_irish_chr",
                                  as.character(chr),
                                  "_plink_clumped.clumps"))
  british_irish_chr_snp_list[[paste0("chr",as.character(chr))]] <- clump_snp_table$ID
  british_irish_snp_count <- british_irish_snp_count + length(clump_snp_table$ID)
}
british_irish_snp_count


################################ Inspection ####################################

irish_chr_snps_all_sig_005 <- Reduce(c, gwas_results)
british_chr_snps_all_sig_0005 <- Reduce(c, british_gwas_results)
length(irish_chr_snps_all_sig_005)
# 5359
length(british_chr_snps_all_sig_0005)
# 26530
length(intersect(irish_chr_snps_all_sig_005, british_chr_snps_all_sig_0005))
# 768

british_chr_snps <- Reduce(c, british_chr_snp_list)
length(intersect(irish_chr_snps_all_sig_005, british_chr_snps))
# 158
length(british_chr_snps)
# 3256

british_irish_chr_snps <- Reduce(c, british_irish_chr_snp_list)
length(intersect(irish_chr_snps_all_sig_005, british_irish_chr_snps))
# 167
length(british_irish_chr_snps)
# 3336

