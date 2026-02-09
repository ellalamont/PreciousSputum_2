# First attempt on predicting outcome
# E. Lamont
# 10/29/25
# 12/11/25: Try with the 60 threshold and also maybe with all the samples


# https://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/



source("Import_data.R") 

# install.packages("caret")
library(caret)
# install.packages("glmnet")
library(glmnet)
library(pROC)


###########################################################
##################### ORGANIZE DATA #######################

# Keep only W0 sputum
W0SputumSamples60_pipeSummary <- GoodSamples60_pipeSummary %>% filter(Type == "Week 0 sputum")

P_metadata <- W0SputumSamples60_pipeSummary %>% select(SampleID2, Outcome)
P_metadata$Outcome <- factor(P_metadata$Outcome)


P_TPM <- GoodSamples60_tpmf %>% select(all_of(W0SputumSamples60_pipeSummary$SampleID2))

# Put everything in one dataframe
P_TPM_t <- P_TPM %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")
my_df <- inner_join(P_metadata, P_TPM_t, by = "SampleID2")

###########################################################
###################### NEW WAY - LASSO ####################

# do a 2/3 split? Other common is 0.8

# Separate the training and testing datasets
set.seed(42)
trainIndex <- createDataPartition(my_df$Outcome, p = 2/3, list = F)
train <- my_df[trainIndex,]
test <- my_df[-trainIndex,]

# Create matrices for glmnet
x_train <- train %>%
  select(-SampleID2, -Outcome) %>%
  as.matrix() # matrix of predictor variables
# y_train <- train$Outcome # the response or outcome variable, which is a binary variable.
y_train <- ifelse(train$Outcome == "Relapse", 1, 0)
x_test <- test %>%
  select(-SampleID2, -Outcome) %>%
  as.matrix()
# y_test <- test$Outcome 
y_test <- ifelse(test$Outcome == "Relapse", 1, 0)

# Check for NAs
sum(is.na(x_train)) # 0
sum(is.na(x_test)) # 0


# https://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/

# Find the best lambda using cross-validation
# lamba: a numeric value defining the amount of shrinkage. Should be specify by analyst.
set.seed(42)
cv.lasso <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial", type.measure = "deviance", nfolds = 5) # Needs to be deviance because the sample size is so small
# Warning messages:
# 1: In lognet(x, is.sparse, y, weights, offset, alpha, nobs,  ... :
#                one multinomial or binomial class has fewer than 8  observations; dangerous ground

plot(cv.lasso)
cv.lasso$lambda.min 
cv.lasso$lambda.1se # Recommended for small n?


coef(cv.lasso, cv.lasso$lambda.1se) # There are too many to see anything with this....


# Fit the final model on the training data
lasso.model <- glmnet(x_train, y_train, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.1se)

# 1/5/26
probs <- predict(lasso.model, newx = x_test, type = "response")
roc_obj <- roc(response = y_test, predictor = as.numeric(probs))
auc(roc_obj)
plot(roc_obj, print.auc = TRUE)


# Make prediction on test data
probabilities <- lasso.model %>% predict(newx = x_test)
predicted.classes <- ifelse(probabilities > 0.5, "relapse", "cure")

# Model accuracy rate
observed.classes <- test$Outcome
mean(predicted.classes == observed.classes) # 0, so not good!


###########################################################
#################### COMPUTE FULL MODEL ###################
# https://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/

# Fit the model
full.model <- glm(Outcome ~., data = x_train, family = binomial)
# Make predictions
probabilities <- full.model %>% predict(test.data, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg")
# Model accuracy
observed.classes <- test.data$diabetes
mean(predicted.classes == observed.classes)


###########################################################
################## W0 USE MONTE CARLO AUC #################
#1/5/26
# Done by ChatGPT!

set.seed(42)

n_repeats <- 200
auc_values <- numeric(n_repeats)

for (i in seq_len(n_repeats)) {
  
  idx <- createDataPartition(my_df$Outcome, p = 0.8, list = FALSE)
  
  train <- my_df[idx, ]
  test  <- my_df[-idx, ]
  
  x_train <- as.matrix(train %>% select(-SampleID2, -Outcome))
  y_train <- ifelse(train$Outcome == "Relapse", 1, 0)
  
  x_test <- as.matrix(test %>% select(-SampleID2, -Outcome))
  y_test <- ifelse(test$Outcome == "Relapse", 1, 0)
  
  # Tune lambda by deviance
  cvfit <- cv.glmnet(
    x_train,
    y_train,
    family = "binomial",
    alpha = 1
  )
  
  model <- glmnet(
    x_train,
    y_train,
    family = "binomial",
    alpha = 1,
    lambda = cvfit$lambda.1se
  )
  
  probs <- predict(model, x_test, type = "response")
  
  roc_obj <- roc(y_test, as.numeric(probs), quiet = TRUE)
  auc_values[i] <- auc(roc_obj)
}

mean(auc_values)
quantile(auc_values, c(0.025, 0.975))


# Permutation test - Randomly choses which samples are relapses.
set.seed(42)
n_perm <- 500
perm_auc <- numeric(n_perm)
for (i in seq_len(n_perm)) {
  
  idx <- createDataPartition(my_df$Outcome, p = 0.8, list = FALSE)
  
  train <- my_df[idx, ]
  test  <- my_df[-idx, ]
  
  x_train <- as.matrix(train %>% select(-SampleID2, -Outcome))
  y_train <- ifelse(train$Outcome == "Relapse", 1, 0)
  
  x_test <- as.matrix(test %>% select(-SampleID2, -Outcome))
  y_test <- ifelse(test$Outcome == "Relapse", 1, 0)
  
  # 🔀 permute labels ONLY in training
  y_perm <- sample(y_train)
  
  cvfit <- cv.glmnet(
    x_train,
    y_perm,
    family = "binomial",
    alpha = 1
  )
  
  model <- glmnet(
    x_train,
    y_perm,
    family = "binomial",
    alpha = 1,
    lambda = cvfit$lambda.1se
  )
  
  probs <- predict(model, x_test, type = "response")
  
  roc_obj <- roc(y_test, as.numeric(probs), quiet = TRUE)
  perm_auc[i] <- auc(roc_obj)
}
mean(perm_auc)

p_value <- mean(perm_auc >= mean(auc_values))
p_value

###########################################################
################## WEEK 2 - ALL SAMPLES ###################

# Keep only W2 sputum and put some limits on it
W2SputumSamples_pipeSummary <- my_pipeSummary %>% 
  filter(Type == "Week 2 sputum") %>%
  filter(Outcome != "Failure") %>%
  filter(Txn_Coverage_f > 10) %>%
  mutate(Outcome2 = case_when(Outcome == "Prob Relapse" ~ "Relapse", TRUE ~ Outcome))
# 28 cures, 14 relapses (including probable relapses)

P_metadata <- W2SputumSamples_pipeSummary %>% select(SampleID2, Outcome2)
P_metadata$Outcome2 <- factor(P_metadata$Outcome2)


P_TPM <- All_tpm_f %>% select(all_of(W2SputumSamples_pipeSummary$SampleID2))

# Put everything in one dataframe
P_TPM_t <- P_TPM %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")
my_df <- inner_join(P_metadata, P_TPM_t, by = "SampleID2")

# do a 2/3 split? Other common is 0.8

# Separate the training and testing datasets
set.seed(42)
trainIndex <- createDataPartition(my_df$Outcome2, p = 0.8, list = F)
train <- my_df[trainIndex,]
test <- my_df[-trainIndex,]

# Create matrices for glmnet
x_train <- as.matrix(train[, !names(train) %in% c("SampleID2", "Outcome2")]) # matrix of predictor variables
# y_train <- train$Outcome # the response or outcome variable, which is a binary variable.
y_train <- ifelse(train$Outcome2 == "Relapse", 1, 0)
x_test <- as.matrix(test[, !names(test) %in% c("SampleID2", "Outcome2")]) 
# y_test <- test$Outcome 
y_test <- ifelse(test$Outcome2 == "Relapse", 1, 0)

# Stratify cross-validation inside glmnet??
# set.seed(42)
# folds <- createFolds(y_train, k = 5, list = F)


# https://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/

# Find the best lambda using cross-validation
# lamba: a numeric value defining the amount of shrinkage. Should be specify by analyst.
set.seed(42)
cv.lasso <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial")
# Warning messages:
# 1: In lognet(x, is.sparse, y, weights, offset, alpha, nobs,  ... :
#                one multinomial or binomial class has fewer than 8  observations; dangerous ground
plot(cv.lasso)
cv.lasso$lambda.min 
cv.lasso$lambda.1se 


# coef(cv.lasso, cv.lasso$lambda.min) # There are too many to see anything with this....


# Fit the final model on the training data
lasso.model <- glmnet(x_train, y_train, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.min)
# Warning message:
#   In storage.mode(x) <- "double" : NAs introduced by coercion

# Make prediction on test data
probabilities <- lasso.model %>% predict(newx = x_test)
predicted.classes <- ifelse(probabilities > 0.5, "relapse", "cure")

# Model accuracy rate
observed.classes <- test$Outcome2
mean(predicted.classes == observed.classes) # 0, so not good!





