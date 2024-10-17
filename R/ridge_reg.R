
# install.packages("glmnet")
library(glmnet)
set.seed(123)

phenotype <- readRDS(file = "/fastscratch/myscratch/xding/UKBB/dat_LDL_British_Irish_pheotype_genotype_british_gwas.rds")

# View(head(british))

phenotype$sex <- as.numeric(phenotype$sex) - 1
# typeof(phenotype$sex)
phenotype_british <- phenotype[which(phenotype$ethnicity == "British"),]
phenotype_irish <- phenotype[which(phenotype$ethnicity == "Irish"),]

# summary(phenotype)
# missingness <- colSums(is.na(phenotype))
# is.na(phenotype)

X_british <- phenotype_british[,-c("#FID","IID","ethnicity","LDL","used.in.pca.calculation")]
# X_irish <- phenotype_irish[,-c("#FID","IID","ethnicity","LDL","used.in.pca.calculation")]

# X_british_matrix <- as.matrix(X_british)
y_british <- phenotype_british$LDL

complete_cases <- complete.cases(X_british, y_british)

# Subset the data to only include complete cases
X_british_complete <- X_british[complete_cases, ]
y_british_complete <- y_british[complete_cases]
dim(X_british_complete)


# X_british_matrix_sparse <- Matrix(X_british_matrix, sparse = TRUE)

n_samples <- nrow(X_british_complete)
british_test_indices <- sample(1:n_samples, size = round(0.1 * n_samples))

# Split the data into training and testing sets
X_train <- X_british_complete[-british_test_indices, ]
X_test  <- X_british_complete[british_test_indices, ]
y_train <- y_british_complete[-british_test_indices]
y_test  <- y_british_complete[british_test_indices]


saveRDS(X_train, "/fastscratch/myscratch/xding/UKBB/british_gwas_clumped_snps/british_X_train.rds")
saveRDS(X_test, "/fastscratch/myscratch/xding/UKBB/british_gwas_clumped_snps/british_X_test.rds")
saveRDS(y_train, "/fastscratch/myscratch/xding/UKBB/british_gwas_clumped_snps/british_y_train.rds")
saveRDS(y_test, "/fastscratch/myscratch/xding/UKBB/british_gwas_clumped_snps/british_y_test.rds")

X_train <- readRDS("/fastscratch/myscratch/xding/UKBB/british_gwas_clumped_snps/british_X_train.rds")
X_test <- readRDS("/fastscratch/myscratch/xding/UKBB/british_gwas_clumped_snps/british_X_test.rds")
y_train <- readRDS("/fastscratch/myscratch/xding/UKBB/british_gwas_clumped_snps/british_y_train.rds")
y_test <- readRDS("/fastscratch/myscratch/xding/UKBB/british_gwas_clumped_snps/british_y_test.rds")

dim(X_train)
dim(X_test)

### 1. Fit OLS model using lm()
ols_model <- lm(y_train ~ ., data = X_train)
summary(ols_model)
typeof(ols_model)
ols_model
saveRDS(ols_model, "/fastscratch/myscratch/xding/UKBB/british_gwas_clumped_snps/british_ols_model.rds")

# Residual standard error: 0.8077 on 301185 degrees of freedom
# Multiple R-squared:  0.1459,	Adjusted R-squared:  0.1366 
# F-statistic: 15.74 on 3268 and 301185 DF,  p-value: < 2.2e-16

y_train_pred_ols <- predict(ols_model, newdata = X_train)
train_mse_ols <- mean((y_train - y_train_pred_ols)^2)

# Calculate MSE for OLS model
y_pred_ols <- predict(ols_model, newdata = X_test)
mse_ols <- mean((y_test - y_pred_ols)^2)


train_mse_ols
# 0.6454492
mse_ols
# 0.6659574

R_sq_train <- 1 - sum((y_train - y_train_pred_ols)^2) / sum((y_train - mean(y_train))^2)
# 0.1458669
R_sq_test <- 1 - sum((y_test - y_pred_ols)^2) / sum((y_test - mean(y_test))^2)
# 0.1276008

### Ridge Models
ridge_model <- glmnet(X_train, y_train, alpha = 0)

