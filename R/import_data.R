library(data.table)
library(stringr)
library(dplyr)
# library(genio)
# install.packages("pgenlibr")
library(pgenlibr)


dat <- readRDS("/users/xding/UKBB/data/UKB_phenotype_sub_09122024.rds")

phenotype_indicator <- readRDS("/dcs05/legacy-dcl01-arking/data/UK_Biobank/active/pheno/ukbPheno_032022.rds")

dat <- merge(dat, phenotype_indicator[,c('IID','used.in.pca.calculation')], by = "IID", all.x = TRUE)
# View(dat %>% filter(!is.na(geno_ethnicity)) %>% group_by(ethnicity) %>% summarise(count = n()))
# only British have confirmed genetic ethnicity

psam_chr1 <- fread("/dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr1.psam")
dat_merged_with_fam <- merge(psam_chr1[,1:2], dat, by = "IID", all.x = TRUE)

###############################data cleaning####################################


dat_LDL_British_Irish <- dat_merged_with_fam %>% 
  filter(used.in.pca.calculation == 1) %>% 
  filter(ethnicity %in% c("British","Irish")) %>%
  select(-geno_ethnicity, -ApoB,-HDL,-Tri) %>%
  filter_all(all_vars(!is.na(.)))

saveRDS(dat_LDL_British_Irish, "./data/LDL/dat_LDL_British_Irish.rds")

cov_LDL_British_Irish <- dat_LDL_British_Irish[,c(2,1,3:15)]
pheno_LDL_British_Irish <- dat_LDL_British_Irish[,c(2,1,16)]

write.table(cov_LDL_British_Irish, file = "./data/LDL/GWAS/input/cov_LDL_British_Irish.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(pheno_LDL_British_Irish, file = "./data/LDL/GWAS/input/pheno_LDL_British_Irish.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

dat_LDL_British <- dat_LDL_British_Irish %>% filter(ethnicity == "British") %>% select(-ethnicity)
cov_LDL_British <- dat_LDL_British[,c(2,1,3:14)]
pheno_LDL_British <- dat_LDL_British[,c(2,1,15)]

dat_LDL_Irish <- dat_LDL_Irish %>% filter(ethnicity == "Irish") %>% select(-ethnicity)
cov_LDL_Irish <- dat_LDL_Irish[,c(2,1,3:14)]
pheno_LDL_Irish <- dat_LDL_Irish[,c(2,1,15)]

# dim(dat_LDL_British)
# dim(dat_LDL_Irish)

write.table(cov_LDL_British, file = "./data/LDL/GWAS/input/cov_LDL_British.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(pheno_LDL_British, file = "./data/LDL/GWAS/input/pheno_LDL_British.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(cov_LDL_Irish, file = "./data/LDL/GWAS/input/cov_LDL_Irish.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(pheno_LDL_Irish, file = "./data/LDL/GWAS/input/pheno_LDL_Irish.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


# all(dat_LDL_British_Irish[,1] == dat_LDL_British_Irish[,2]) confirming all IID and FIDs are the same

# check the size of the datasets
# dat_LDL_British_Irish %>% group_by(ethnicity) %>% summarise(count = n())
# # A tibble: 2 x 2
# ethnicity  count
# <ord>      <int>
#   1 British   338389
# 2 Irish      10123

