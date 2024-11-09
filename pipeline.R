library(bigreadr)
library(stringr)
library(dplyr)

##### read cov data
dat = readRDS('/dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Phenotype/Downloads/ukbCombined_032022_mortality_covid.rds')
col_names <- colnames(dat)
geno_ethnic_group_col <- col_names[str_detect(col_names, "22006")]

# first ten principal components
top_ten_pc_col <- col_names[str_detect(col_names, "22009")][1:10]

# 30780 LDL direct - Biomarker
ldl_col <- col_names[str_detect(col_names, "30780")]
apob_col <- col_names[str_detect(col_names, "30640")]
hdl_col <- col_names[str_detect(col_names, "30760")]
tri_col <- col_names[str_detect(col_names, "30870")]
bmi_col <- col_names[str_detect(col_names, "21001")]
height_col <- col_names[str_detect(col_names, "f.50.0.0")]
# select the columns to include for the phenotype and covariate file
all_columns_to_include <- c("f.eid","f.31.0.0", "f.21022.0.0",
                            "f.21000.0.0",geno_ethnic_group_col,
                            top_ten_pc_col,
                            ldl_col[1],
                            apob_col[1],
                            hdl_col[1],
                            tri_col[1],
                            bmi_col[1],
                            height_col[1])

dat_sub <- dat[,all_columns_to_include]
# rename the columns for readability
psam_chr1 <- fread2("/dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr1.psam")
colnames(dat_sub) <- c("IID", "sex", "age", "ethnicity", "geno_ethnicity",
                       paste0("PC",1:10), "LDL", "ApoB", "HDL", "Tri", "BMI", "height")
phenotype_indicator <- readRDS("/dcs05/legacy-dcl01-arking/data/UK_Biobank/active/pheno/ukbPheno_032022.rds")
# phenotype_indicator
dat <- merge(dat_sub, phenotype_indicator[, c('IID', 'used.in.pca.calculation')], by = "IID", all.x = TRUE)
dat_merged_with_fam <- merge(psam_chr1[,1:2], dat, by = "IID", all.x = TRUE)

dat_merged_with_fam_irish <- dat_merged_with_fam %>% filter(ethnicity == "Irish") %>% filter(used.in.pca.calculation == 1)
ID = dat_merged_with_fam_irish[, c(2, 1)]
write.table(ID, file = "/users/ydun/Irish/id.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
irish_cov = dat_merged_with_fam_irish %>% select(`#FID`, IID, sex, age, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>%
  mutate(sex = ifelse(sex == "Male", 1, 0))
pheno = dat_merged_with_fam_irish %>% select(`#FID`, IID, HDL)
names(pheno)[3] = "trait"
pheno = pheno[!is.na(pheno[,3]),]
readr::write_tsv(irish_cov, "/users/ydun/Irish/cov.txt")
readr::write_tsv(pheno, "/users/ydun/Irish/HDL.txt")

##### get training, tuning, validation ID
set.seed(123)
id = read.table("/users/ydun/Irish/id.txt")
n <- nrow(id)
vec_shuffled <- sample(id[,1])
train_size <- round(0.7 * n)
tune_size <- round(0.15 * n)
train <- vec_shuffled[1:train_size]
tune <- vec_shuffled[(train_size + 1):(train_size + tune_size)]
validation <- vec_shuffled[(train_size + tune_size + 1):n]
train_id = data.frame(FID = train, IID = train)
tuning_id = data.frame(FID = tune, IID = tune)
validation_id = data.frame(FID = validation, IID = validation)
readr::write_tsv(train_id, "/users/ydun/Irish/train.id", col_names = FALSE)
readr::write_tsv(tuning_id, "/users/ydun/Irish/tuning.id", col_names = FALSE)
readr::write_tsv(validation_id, "/users/ydun/Irish/validation.id", col_names = FALSE)

##### get genotype data
id = read.table("/users/ydun/Irish/id.txt")
for (chr in 1:22) {
  plink_code = paste(
    "/dcs04/nilanjan/data/ydun/tools/plink2",
    paste0("--pfile /dcl01/arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr", chr),
    "--make-bed ",
    paste0("--extract /dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb/chr", chr, "/snplist_merged.txt"),
    "--keep /users/ydun/Irish/id.txt",
    "--out", paste0("/users/ydun/Irish/genotype/chr", chr, " &"))
  cat(plink_code, "\n")
}
print('for i in {1..22}; do echo "/users/ydun/Irish/genotype/chr${i}" >> /users/ydun/Irish/genotype/merge_list.txt; done')
system("/dcs04/nilanjan/data/ydun/tools/plink --merge-list /users/ydun/Irish/genotype/merge_list.txt --make-bed --out /users/ydun/Irish/genotype/all")

###### quality control
plink_code = paste(
  "/dcs04/nilanjan/data/ydun/tools/plink2",
  paste0("--bfile /users/ydun/Irish/genotype/all"),
  "--make-bed ",
  "--geno 0.05 --mind 0.05 --hwe 1.0e-7 --maf 0.05 --mach-r2-filter 0.8 2 --snps-only ",
  "--out", paste0("/users/ydun/Irish/genotype/all_qc &"))
cat(plink_code, "\n")


##### get adjusted phenotype in training
cov = fread2("/users/ydun/Irish/cov.txt")
pheno = fread2("/users/ydun/Irish/HDL.txt")
dat = merge(pheno, cov, by = c("#FID", "IID"))
model = lm(dat[[3]] ~ as.matrix(dat[,c(-1, -2, -3)]), data = dat %>% filter(IID %in% train_id[,2]))
pheno = data.frame(FID = dat[,1], IID = dat[,1], y = model$residuals)
trait = "HDL"
readr::write_tsv(pheno, paste0("/users/ydun/Irish/", trait, "/HDL_adjusted.txt"))


##### GWAS
trait = "HDL"
system(paste0("mkdir -p /users/ydun/Irish/", trait, "/GWAS/"))
GWAS_code = paste(paste0('/dcl01/chatterj/data/tools/plink2'),
                  paste0('--bfile /users/ydun/Irish/genotype/all_qc'), # genotype data file (PLINK binary files)
                  paste0('--keep /users/ydun/Irish/train.id'),
                  paste0('--pheno /users/ydun/Irish/', trait, '/HDL_adjusted.txt'), # phenotype data file
                  paste0('--glm '),
                  paste0('--out ', "/users/ydun/Irish/", trait, "/GWAS/HDL_GWAS"))
print(paste0(GWAS_code, " &"))


##### clumping and thresholding
system(paste0("awk '{print $3,$11}' ", "/users/ydun/Irish/", trait, "/GWAS/HDL_GWAS.y.glm.linear > ", "/users/ydun/Irish/", trait, "/GWAS/main.SNP.pvalue"))

### clumping and thresholding
pval_list = c(5e-8, 5e-7, 5e-6, 5e-5, 5e-4, 5e-3, 5e-2)
system(paste0("mkdir -p /users/ydun/Irish/", trait, "/C_T/"))
system(paste0("rm /users/ydun/Irish/", trait, "/C_T/range_list"))
system(paste0("echo '", pval_list[1]," 0 ", pval_list[1], "' > ", "/users/ydun/Irish/", trait, "/C_T//range_list "))
for (pval_i in 1:length(pval_list)) {
  system(paste0("echo '", pval_list[pval_i]," 0 ", pval_list[pval_i], "' >> ", "/users/ydun/Irish/", trait, "/C_T//range_list "))
}
r2_list = c(0.1, 0.2)
for (r2 in r2_list) {
  sumstat = paste0("/users/ydun/Irish/", trait, "/GWAS/HDL_GWAS.y.glm.linear")
  clump_cmd = paste(
    "/dcl01/chatterj/data/tools/plink",
    paste0('--bfile /users/ydun/Irish/genotype/all_qc'),
    "--clump", sumstat,
    "--clump-p1", 1,
    "--clump-r2", r2,
    "--clump-kb", 250,
    "--clump-snp-field ID",
    "--clump-field P",
    "--out", paste0("/users/ydun/Irish/", trait, "/C_T/clump_", "_r2_", r2))
  system(clump_cmd)
}
for (r2 in r2_list) {
  system(paste0("awk 'NR!=1{print $3}' ", "/users/ydun/Irish/", trait, "/C_T/clump_", "_r2_", r2, ".clumped >  /users/ydun/Irish/", trait, "/C_T/valid.r2", r2, ".snp"))
  sumstat = paste0("/users/ydun/Irish/", trait, "/GWAS/HDL_GWAS.y.glm.linear")
  prs_cmd = paste(
    "/dcl01/chatterj/data/tools/plink",
    paste0('--bfile /users/ydun/Irish/genotype/all_qc'), # genotype data file (PLINK binary files)
    paste0('--keep /users/ydun/Irish/tuning.id'),
    "--score", sumstat, "3 5 8 sum header",
    "--q-score-range ", paste0("/users/ydun/Irish/", trait, "/C_T/range_list ", "/users/ydun/Irish/", trait, "/GWAS/main.SNP.pvalue"),
    "--extract ", paste0("/users/ydun/Irish/", trait, "/C_T/valid.r2", r2, ".snp"),
    "--out", paste0("/users/ydun/Irish/", trait, "/C_T/prs_r2_", r2))
    system(prs_cmd)
}


### evaluation
library(dplyr)
cov = fread2("/users/ydun/Irish/cov.txt")
pheno = fread2("/users/ydun/Irish/HDL.txt")
dat = merge(pheno, cov, by = c("#FID", "IID"))
id <- bigreadr::fread2("/users/ydun/Irish/tuning.id")
dat <- dat[dat[,2] %in% id[,2],]
fit <- lm(dat[[3]] ~ as.matrix(dat[,c(-1,-2, -3)]), na.action = na.exclude)
yt <- data.frame(id = dat[,2], yt = resid(fit))
yt <- yt[!is.na(yt$yt),]
pheno <- data.frame(FID=yt$id, IID=yt$id, y=yt$yt)
tuning_result = data.frame(matrix(ncol = 3, nrow = 0))
names(tuning_result) = c('r2', 'pval_threshold', 'p-val'); tuning_result_i = 1

for (r2 in r2_list) {
  for (pval in pval_list) {
    prs_main = bigreadr::fread2(paste0("/users/ydun/Irish/", trait, "/C_T/prs_r2_", r2, ".", pval, ".profile")) %>%
      select(IID, SCORESUM)
    prs_table_all = merge(prs_main, pheno)
    my_lm = lm(y ~ SCORESUM, data = prs_table_all)
    tuning_result[tuning_result_i,] = c(r2, pval, summary(my_lm)$r.squared); tuning_result_i = tuning_result_i + 1
  }
}
