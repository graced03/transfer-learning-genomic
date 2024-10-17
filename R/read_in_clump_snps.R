library(pgenlibr)
library(data.table)
library(dplyr)
library(Matrix)

# File paths
# psam_path <- "./data/geno/chr22.psam"
# Read in the sample data
# sample_data <- fread(psam_path)
# dim(sample_data)
# 
# 
# pgen_path <- "./data/geno/chr22.pgen"
# pvar_path <- "./data/geno/chr22.pvar"
# pvar_table <- fread(pvar_path)
# pvar <- NewPvar(pvar_path)
# pgen <- NewPgen(pgen_path, pvar=pvar)
# 
# dim(pvar_table) # 13355     6

clump_snp <- fread("./data/LDL/clump/british_snp_clump/ldl_british_chr1_plink_clumped.clumps")

########################### Snps selected by gwas on both ######################
############################## British and Irish ###############################

british_irish_chr_snp_list <- list()
british_irish_snp_count <- 0
for (chr in 1:22){
  clump_snp_table <- fread(paste0("./data/LDL/clump/british_irish_snp_clump/ldl_british_irish_chr",
                                  as.character(chr),
                                  "_plink_clumped.clumps"))
  british_irish_chr_snp_list[[paste0("chr",as.character(chr))]] <- clump_snp_table$ID
  british_irish_snp_count <- british_irish_snp_count + length(clump_snp_table$ID)
}


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

################################## Irish #######################################

irish_chr_snp_list <- list()
irish_snp_count <- 0
for (chr in 1:22){
  clump_snp_table <- fread(paste0("./data/LDL/clump/irish_snp_clump/ldl_irish_chr",
                                  as.character(chr),
                                  "_plink_clumped.clumps"))
  irish_chr_snp_list[[paste0("chr",as.character(chr))]] <- clump_snp_table$ID
  irish_snp_count <- irish_snp_count + length(clump_snp_table$ID)
}

########################### Look up the intersection############################


british_chr_snps <- Reduce(c, british_chr_snp_list)
irish_chr_snps <- Reduce(c, irish_chr_snp_list)
british_irish_chr_snps <- Reduce(c, british_irish_chr_snp_list)

length(intersect(british_chr_snps, british_irish_chr_snps)) # 2904
length(intersect(british_chr_snps, irish_chr_snps)) # 56
length(intersect(british_irish_chr_snps, irish_chr_snps)) # 58

length(union(british_chr_snps, british_irish_chr_snps)) # 3688
length(union(british_chr_snps, irish_chr_snps)) # 3508
length(union(british_irish_chr_snps, irish_chr_snps)) # 3586
length(union(union(british_irish_chr_snps, irish_chr_snps),british_chr_snps)) # 3937


########################### Extract the SNP data ###############################
# Phenotype file
phenotype <- readRDS("~/UKBB/data/LDL/phenotype/dat_LDL_British_Irish.rds")
# Genotype file directory
geno_path <- "./data/geno/"

for (i in 1:22){
  chr <- paste0("chr", as.character(i))
  psam_path <- paste0(geno_path,chr,".psam")
  sample_data <- fread(psam_path)
  
  pgen_path <- paste0(geno_path,chr,".pgen")
  pvar_path <- paste0(geno_path,chr,".pvar")
  pvar <- NewPvar(pvar_path)
  pgen <- NewPgen(pgen_path, pvar=pvar)
  
  pvar_table <- fread(pvar_path)
  pvar_table$Index <- 1:nrow(pvar_table)
  
  snp_selected_chr <- pvar_table %>% filter(ID %in% british_chr_snp_list[[chr]])
  geno_mat <- ReadList(pgen, snp_selected_chr$Index, meanimpute=F)
  colnames(geno_mat) <- british_chr_snp_list[[chr]]
  
  merged <- cbind(sample_data[,2], geno_mat)
  phenotype <- merge(phenotype, merged, by = "IID", all.x = TRUE)
}

summary(phenotype)
missingness <- colSums(is.na(phenotype))
sumZero <- colSums(phenotype == 0,na.rm = TRUE)

sum(phenotype==0, na.rm=TRUE)
dim(phenotype)[1] * dim(phenotype)[2]

dim(phenotype)
View(phenotype)

# saveRDS(phenotype, file = "~/UKBB/data/LDL/clump_selected_snp_w_phenotype/british_irish_snp_clump/dat_LDL_British_Irish_pheotype_genotype_british_irish_gwas.rds")
saveRDS(phenotype, file = "/fastscratch/myscratch/xding/UKBB/dat_LDL_British_Irish_pheotype_genotype_british_gwas.rds")




