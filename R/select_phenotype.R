library(data.table)
library(stringr)
library(dplyr)

dat <- readRDS("/dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Phenotype/Downloads/ukbCombined_032022_mortality_covid.rds")

col_names <- colnames(dat)
# col_names[1] "f.eid"
# dat[,"f.eid"]
# the first column is the IID of the patients, which is the unique identifier

# the following command find the column names for sex, age and ethnicity of the patients
# col_names[str_detect(col_names, "f.31")] # sex
# col_names[str_detect(col_names, "21022")] # age
# col_names[str_detect(col_names, "21000")] # ethnicity

sex <- dat[,"f.31.0.0"]
age <- dat[,"f.21022.0.0"]
ethnicity <- dat[,"f.21000.0.0"]

# Ethnic group
# Genetic ethnic group, indicates samples who self-identified as 'White British' according to Field 21000 
# and have very similar genetic ancestry based on a principal components analysis of the genotypes.
geno_ethnic_group_col <- col_names[str_detect(col_names, "22006")]
# geno_ethnic_group <- dat[,c("f.eid",geno_ethnic_group_col,"f.21000.0.0")]

# first ten principal components
top_ten_pc_col <- col_names[str_detect(col_names, "22009")][1:10]
top_ten_pc <- dat[,top_ten_pc_col]

# 23405 LDL Cholesterol - NMR metabolomics
# 30780 LDL direct - Biomarker
ldl_col <- col_names[str_detect(col_names, "30780")]
ldl <- dat[,ldl_col]

# 30640 Apolipoprotein B
# 23439 Apolipoprotein B
apob_col <- col_names[str_detect(col_names, "30640")]
apob <- dat[,apob_col]

# 30760	HDL cholesterol
hdl_col <- col_names[str_detect(col_names, "30760")]
hdl <- dat[,hdl_col[1]]

# 30870	Triglycerides
tri_col <- col_names[str_detect(col_names, "30870")]
tri <- dat[,tri_col]

# select the columns to include for the phenotype and covariate file
all_columns_to_include <- c("f.eid","f.31.0.0", "f.21022.0.0", 
                            "f.21000.0.0",geno_ethnic_group_col,
                            top_ten_pc_col,
                            ldl_col[1],apob_col[1],hdl_col[1],tri_col[1])
dat_sub <- dat[,all_columns_to_include]
# rename the columns for readability
colnames(dat_sub) <- c("IID","sex","age","ethnicity","geno_ethnicity",
                       paste0("PC",1:10),"LDL","ApoB","HDL","Tri")


saveRDS(dat_sub, "/users/xding/UKBB/data/UKB_phenotype_sub_09122024.rds")

################################################################################
############################### Visualization ##################################
################################################################################
View(geno_ethnic_group)

dat_count <- dat %>%
  group_by(f.22006.0.0, f.21000.0.0) %>%
  summarise(count = n()) %>%
  ungroup()
View(dat_count)

View(dat %>% group_by(f.21000.0.0) %>% summarise(count = n()))
View(dat %>% group_by(f.22006.0.0) %>% summarise(count = n()))
