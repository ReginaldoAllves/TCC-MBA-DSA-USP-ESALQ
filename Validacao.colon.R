#### LIBRARIES

library('caret')             # Training, hyper parameter tuning and prediction
library('caretEnsemble')     # Ensemble of caret models
library('tidyverse')         # Data wrangling
library('skimr')             # Descriptive statistics
library('GGally')            # Descriptive statistics
library('randomForest')      # Required to build Random Forest models
library('gbm')               # Required to build Binomial Logistic Regression models
library('nnet')              # Required to build Neural Networks
library('klaR')              # Required to build Naive Bayes 
library('ggcorrplot')        # Required to draw a correlation matrix
library('pROC')              # Required to draw ROC curves  

#### DATASET

coad.all = read.table("./data/colon.data.txt", 
                      fill = TRUE, sep = "\t", header = TRUE, quote = "")

## Removal of Patient.ID

coad.all = dplyr::select(coad.all, -1)

## Removal of engineered features

coad.all = dplyr::select(coad.all,
                         -all_of(c(6,7)))

## ## For the minimal feature model

coad.all = dplyr::select(coad.all,
                         -all_of(c(1:5)))

## ## For the minimal feature model with frameshift

coad.all = dplyr::select(coad.all,
                         -all_of(c(1,2)),
                         -all_of(c(4,5)))

sapply(coad.all, anyNA)
any(duplicated(coad.all))

coad.data = coad.all

str(coad.data)

coad.data$msi = factor(coad.data$msi, 
                       levels = c(0,1), 
                       labels = c("MSS", "MSI"))


prop.table(table(coad.data$msi))

#set.seed(1)
#coad.bagMissing = preProcess(coad.data, method = "bagImpute")
#coad.data = predict(coad.bagMissing, newdata = coad.data)

coad.rangeModel <- preProcess(coad.data, method = "range")
coad.data <- predict(coad.rangeModel, newdata = coad.data)

## Reading classifiers

rf.all.model = readRDS('./classifiers/rf.all.model.rds')
nb.all.model = readRDS('./classifiers/nb.all.model.rds')
nn.all.model = readRDS('./classifiers/nn.all.model.rds')
lr.all.model = readRDS('./classifiers/lr.all.model.rds')

rf.min.model = readRDS('./classifiers/rf.min.model.rds')
nb.min.model = readRDS('./classifiers/nb.min.model.rds')
nn.min.model = readRDS('./classifiers/nn.min.model.rds')
lr.min.model = readRDS('./classifiers/lr.min.model.rds')


## Random Forest

confusionMatrix(reference = coad.data$msi, 
                data = coad.fittedRF, 
                mode = "everything", 
                positive = "MSI")

## Naive Bayes

confusionMatrix(reference = coad.data$msi, 
                data = coad.fittedNB, 
                mode = "everything", 
                positive = "MSI")

## Neural network

confusionMatrix(reference = coad.data$msi, 
                data = coad.fittedNN, 
                mode = "everything", 
                positive = "MSI")

## LR

confusionMatrix(reference = coad.data$msi, 
                data = coad.fittedLR, 
                mode = "everything", 
                positive = "MSI")


#### ROC curves

## Prediction with probabilitie

#### PREDICT

prob.all.fittedRF <- predict(rf.all.model, coad.data, type = "prob")
prob.all.fittedNB <- predict(nb.all.model, coad.data, type = "prob")
prob.all.fittedNN <- predict(nn.all.model, coad.data, type = "prob")
prob.all.fittedLR <- predict(lr.all.model, coad.data, type = "prob")

prob.min.fittedRF <- predict(rf.min.model, coad.data, type = "prob")
prob.min.fittedNB <- predict(nb.min.model, coad.data, type = "prob")
prob.min.fittedNN <- predict(nn.min.model, coad.data, type = "prob")
prob.min.fittedLR <- predict(lr.min.model, coad.data, type = "prob")

## ROC

roc.all.fittedRF <- roc(coad.data$msi, prob.all.fittedRF$MSS)
roc.all.fittedNB <- roc(coad.data$msi, prob.all.fittedNB$MSS)
roc.all.fittedNN <- roc(coad.data$msi, prob.all.fittedNN$MSS)
roc.all.fittedLR <- roc(coad.data$msi, prob.all.fittedLR$MSS)

roc.min.fittedRF <- roc(coad.data$msi, prob.min.fittedRF$MSS)
roc.min.fittedNB <- roc(coad.data$msi, prob.min.fittedNB$MSS)
roc.min.fittedNN <- roc(coad.data$msi, prob.min.fittedNN$MSS)
roc.min.fittedLR <- roc(coad.data$msi, prob.min.fittedLR$MSS)
