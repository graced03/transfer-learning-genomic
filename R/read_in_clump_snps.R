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

####################################################################################################

phenotype <- readRDS("~/Documents/UKBB/british_irish_gwas_clumped_snps/dat_LDL_British_Irish_pheotype_genotype_british_irish_gwas.rds")

phenotype$sex <- as.numeric(phenotype$sex) - 1
typeof(phenotype$sex)
phenotype <- as.data.frame(phenotype)
phenotype <- phenotype[complete.cases(phenotype),]
dim(phenotype)
phenotype_british <- phenotype[which(phenotype$ethnicity == "British"),]
phenotype_irish <- phenotype[which(phenotype$ethnicity == "Irish"),]

dim(phenotype_british)
dim(phenotype_irish)

X_british <- phenotype_british[, -which(names(phenotype_british) %in% c("#FID","IID","ethnicity","LDL","used.in.pca.calculation"))]

X_irish <- phenotype_irish[, -which(names(phenotype_irish) %in% c("#FID","IID","ethnicity","LDL","used.in.pca.calculation"))]
# dim(X_irish[complete.cases(X_irish), ])
# dim(X_irish)

dim(X_british)
# 338282    320
# 338282   3348

X_british_matrix <- as.matrix(X_british)
X_british_matrix_sparse <- Matrix(X_british_matrix, sparse = TRUE)
# X_british_matrix_sparse <- as(X_british_matrix_sparse, "RsparseMatrix")
y_british <- phenotype_british$LDL
length(y_british)

X_irish_matrix <- as.matrix(X_irish)
X_irish_matrix_sparse <- Matrix(X_irish_matrix, sparse = TRUE)
# X_irish_matrix_sparse <- as(X_irish_matrix_sparse, "RsparseMatrix")
# 10123   320
# 10123  3353
y_irish <- phenotype_irish$LDL

# saveRDS(X_british_matrix_sparse, "~/Documents/UKBB/X_british_matrix_sparse_LDL_pheotype_genotype_irish_gwas.rds")

saveRDS(y_british, "~/Documents/UKBB/british_irish_gwas_clumped_snps/y_british.rds")
saveRDS(y_irish, "~/Documents/UKBB/british_irish_gwas_clumped_snps/y_irish.rds")

saveRDS(X_british_matrix_sparse@i, "~/Documents/UKBB/british_irish_gwas_clumped_snps/X_british_matrix_sparse_j.rds")
saveRDS(X_british_matrix_sparse@p, "~/Documents/UKBB/british_irish_gwas_clumped_snps/X_british_matrix_sparse_p.rds")
saveRDS(X_british_matrix_sparse@x, "~/Documents/UKBB/british_irish_gwas_clumped_snps/X_british_matrix_sparse_x.rds")

saveRDS(X_irish_matrix_sparse@i, "~/Documents/UKBB/british_irish_gwas_clumped_snps/X_irish_matrix_sparse_i.rds")
saveRDS(X_irish_matrix_sparse@p, "~/Documents/UKBB/british_irish_gwas_clumped_snps/X_irish_matrix_sparse_p.rds")
saveRDS(X_irish_matrix_sparse@x, "~/Documents/UKBB/british_irish_gwas_clumped_snps/X_irish_matrix_sparse_x.rds")


# saveRDS(X_irish, "~/Documents/UKBB/X_british_matrix_sparse_LDL_pheotype_genotype_irish_gwas.rds")


intersect(colnames(dat_LDL_British_Irish_pheotype_genotype_irish_gwas)[18:325],
          colnames(dat_LDL_British_Irish_pheotype_genotype_british_irish_gwas)[18:3353]) # 58 intersections


