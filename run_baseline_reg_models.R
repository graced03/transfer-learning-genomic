library(glmnet)
library(data.table)
library(dplyr)

#### load in data
extracted_geno_data_irish_gwas <- readRDS(file = "/fastscratch/myscratch/xding/UKBB/irish_gwas/dat_LDL_British_Irish_pheotype_genotype_irish_gwas_all_snps_p005.rds")
dim(extracted_geno_data_irish_gwas)

complete_cases <- complete.cases(extracted_geno_data_irish_gwas)
extracted_geno_data_irish_gwas <- extracted_geno_data_irish_gwas[complete_cases, ]
extracted_geno_data_irish_gwas$sex <- as.numeric(extracted_geno_data_irish_gwas$sex) - 1

print(dim(extracted_geno_data_irish_gwas))

#### split into british and irish data
phenotype_british <- extracted_geno_data_irish_gwas[which(extracted_geno_data_irish_gwas$ethnicity == "British"),]
X_british <- phenotype_british %>%
  select(-IID, -`#FID`, -ethnicity, -LDL, -used.in.pca.calculation)
y_british <- phenotype_british$LDL

n_samples <- nrow(X_british_complete)
british_test_indices <- sample(1:n_samples, size = round(0.1 * n_samples))

# Split the data into training and testing sets
X_train <- X_british_complete[-british_test_indices, ]
X_test  <- X_british_complete[british_test_indices, ]
y_train <- y_british_complete[-british_test_indices]
y_test  <- y_british_complete[british_test_indices]

ols_model <- lm(y_train ~ ., data = X_train)

y_train_pred_ols <- predict(ols_model, newdata = X_train)
train_mse_ols <- mean((y_train - y_train_pred_ols)^2)

# Calculate MSE for OLS model
y_pred_ols <- predict(ols_model, newdata = X_test)
mse_ols <- mean((y_test - y_pred_ols)^2)

R_sq_train <- 1 - sum((y_train - y_train_pred_ols)^2) / sum((y_train - mean(y_train))^2)
R_sq_test <- 1 - sum((y_test - y_pred_ols)^2) / sum((y_test - mean(y_test))^2)


saveRDS(ols_model$intercept, "/fastscratch/myscratch/xding/UKBB/british_gwas_clumped_snps/british_ols_model.rds")

saveRDS(ols_model$coef, "/fastscratch/myscratch/xding/UKBB/british_gwas_clumped_snps/british_ols_model.rds")

# phenotype_irish <- extracted_geno_data_irish_gwas[which(extracted_geno_data_irish_gwas$ethnicity == "Irish"),]
# X_irish <- phenotype_irish %>%
#   select(-IID, -`#FID`, -ethnicity, -LDL, -used.in.pca.calculation)

# y_irish <- phenotype_irish$LDL





