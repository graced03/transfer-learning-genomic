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

# clump_snp <- fread("./data/LDL/clump/irish_snp_clump/ldl_irish_chr1_plink_clumped.clumps")

chr_snp_list <- list()
snp_count <- 0
for (chr in 1:22){
  clump_snp_table <- fread(paste0("./data/LDL/clump/irish_snp_clump/ldl_irish_chr",
                                  as.character(chr),
                                  "_plink_clumped.clumps"))
  chr_snp_list[[paste0("chr",as.character(chr))]] <- clump_snp_table$ID
  snp_count <- snp_count + length(clump_snp_table$ID)
}

# Phenotype file
phenotype <- readRDS("~/UKBB/data/LDL/dat_LDL_British_Irish.rds")
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
  
  snp_selected_chr <- pvar_table %>% filter(ID %in% chr_snp_list[[chr]])
  geno_mat <- ReadList(pgen, snp_selected_chr$Index, meanimpute=F)
  colnames(geno_mat) <- chr_snp_list[[chr]]
  
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
saveRDS(phenotype, file = "~/UKBB/data/LDL/dat_LDL_British_Irish_pheotype_genotype_irish_gwas.rds")

phenotype <- readRDS("~/UKBB/data/LDL/dat_LDL_British_Irish_pheotype_genotype_irish_gwas.rds")

phenotype$sex <- as.numeric(phenotype$sex) - 1
typeof(phenotype$sex)
phenotype_british <- phenotype[which(phenotype$ethnicity == "British"),]
phenotype_irish <- phenotype[which(phenotype$ethnicity == "Irish"),]

X_british <- phenotype_british[,-c("#FID","IID","ethnicity","LDL","used.in.pca.calculation")]
X_irish <- phenotype_irish[,-c("#FID","IID","ethnicity","LDL","used.in.pca.calculation")]

X_british_matrix <- as.matrix(X_british)
X_british_matrix_sparse <- Matrix(X_british_matrix, sparse = TRUE)

saveRDS(X_british_matrix_sparse, "~/UKBB/data/LDL/X_british_matrix_sparse_LDL_pheotype_genotype_irish_gwas.rds")

head(X)
dim(X)

X[,1]



