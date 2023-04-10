rm(list = ls())

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

ucec.all = read.table("./data/utero.data.txt", 
                      fill = TRUE, sep = "\t", header = TRUE, quote = "")

#### DATA WRANGLINH

## All molecular features (modelo completo)

ucec.all = dplyr::select(ucec.all, -Patient.ID)


## For the minimal feature model (modelo mínimo)

ucec.all = dplyr::select(ucec.all,
                         -Patient.ID,
                         -all_of(c(2:6)))

## Check for NAs

sapply(ucec.all, anyNA)

# Quantification of NAs

missingData <- ucec.all %>%
  summarise_all(funs(sum(is.na(.)))) %>% 
  gather("column") %>%
  rename(NumNAs = value) %>% 
  mutate(PrcNAs = NumNAs/nrow(ucec.all)) %>% 
  filter(NumNAs!=0) %>%
  arrange(desc(PrcNAs))

missingData

## Check for duplicated rows

any(duplicated(ucec.all))

#### SELECT DATA
  
ucec.data = ucec.all

#### DATA FORMATING

##  Converting numeric variables into factors

str(ucec.data)
ucec.data$msi = factor(ucec.data$msi, 
                       levels = c(0,1), 
                       labels = c("MSS", "MSI"))

## Near-Zero variables

exclude.cols = nearZeroVar(ucec.data, freqCut = 95/5, uniqueCut = 10)

#### DATASET SPLIT (75:25)

set.seed(1)
index = createDataPartition(ucec.data$msi, p = .75, list = FALSE)
trainingSet = ucec.data[index,]
testSet = ucec.data[-index,]

# Checking proportion of response variable

prop.table(table(ucec.data$msi))
prop.table(table(trainingSet$msi))
prop.table(table(testSet$msi))

#### PRE-PROCESSING

## Data normalization (0-1)

rangeModel <- preProcess(trainingSet, method = "range")
trainingSet <- predict(rangeModel, newdata = trainingSet)

#### PRE-PROCESS TEST DATASET

testSet <- predict(rangeModel, testSet)

#### TRAINING AND HYPERPARAMETER TUNING

### Control function

Ctrl <- trainControl(
  method = "repeatedcv",              # Repeated k-fold cross-validation
  number = 5,                         # Number of K
  repeats = 5,                        # Number of repeats
  savePredictions = "final",          # Save predictions of the final model only
  classProbs = TRUE,                  # Save probabilities
  summaryFunction = twoClassSummary   # Obtain AUC, Sensibility, and Specificity
)

### Models

## Random forest

set.seed(1)
rf.min.model <- train(msi ~ ., 
                  data = trainingSet, 
                  method = "rf", 
                  metric = "ROC", 
                  trControl = Ctrl, 
                  tuneLength = 10)

## Naive bayes

set.seed(1)
nb.min.model <- train(msi ~ ., 
                  data = trainingSet, 
                  method = "nb", 
                  metric = "ROC", 
                  trControl = Ctrl, 
                  tuneLength = 10)

## Neural network

set.seed(1)
nn.min.model <- train(msi ~ ., 
                  data = trainingSet, 
                  method = "nnet", 
                  metric = "ROC", 
                  trControl = Ctrl, 
                  tuneLength = 10)

## Logistic regression with step wise

set.seed(1)
lr.min.model <- train(msi ~ ., 
                   data = trainingSet, 
                   method = "glmStepAIC", 
                   metric = "ROC",
                   family = "binomial",
                   trControl = Ctrl, 
                   tuneLength = 10)

summary(lr.model)


## Plotting model performance

rfPerf.plot = plot(rf.model, main = "Model AUC with Random Forest")
nbPerf.plot = plot(nb.model, main = "Model AUC with Naive Bayes")
nnPerf.plot = plot(nn.model, main = "Model AUC with Neural Network")

summary(lr.model)

## Plotting variance importance

rfVarImp.plot = plot(varImp(rf.model), main="Variable Importance with Random Forest")
nbVarImp.plot = plot(varImp(nb.model), main="Variable Importance with Naive Bayes")
nnVarImp.plot = plot(varImp(nn.model), main="Variable Importance with Neural Network")
lr.plot = plot(?varImp(lr.model), main="Variable Importance with Logistic Regression")

#### PREDICT

fittedRF <- predict(rf.model, testSet)
fittedNB <- predict(nb.model, testSet)
fittedNN <- predict(nn.model, testSet)
fittedLR <- predict(lr.model, testSet)

#### CONFUSION MATRIX

## Random Forest

confusionMatrix(reference = testSet$msi, 
                  data = fittedRF, 
                  mode = "everything", 
                  positive = "MSI")

## Naive Bayes

confusionMatrix(reference = testSet$msi, 
                data = fittedNB, 
                mode = "everything", 
                positive = "MSI")


## Neural network

confusionMatrix(reference = testSet$msi, 
                data = fittedNN, 
                mode = "everything", 
                positive = "MSI")

## LR

confusionMatrix(reference = testSet$msi, 
                data = fittedLR, 
                mode = "everything", 
                positive = "MSI")


#### COMPARING PERFORMANCES

models_compare <- resamples(list(
  RandomForest=rf.model,
  NaiveBayes=nb.model,
  NeuralNetwork=nn.model,
  LogisticRegression=lr.model))


summary(models_compare)

scales <- list(x=list(relation="free"), y=list(relation="free"))

bwplot(models_compare, scales=scales)

#### AUC

## Prediction with probabilitie

#### PREDICT

prob.fittedRF <- predict(rf.model, testSet, type = "prob")
prob.fittedNB <- predict(nb.model, testSet, type = "prob")
prob.fittedNN <- predict(nn.model, testSet, type = "prob")
prob.fittedLR <- predict(lr.model, testSet, type = "prob")

## ROC

roc.fittedRF <- roc(testSet$msi, prob.fittedRF$MSS)
roc.fittedNB <- roc(testSet$msi, prob.fittedNB$MSS)
roc.fittedNN <- roc(testSet$msi, prob.fittedNN$MSS)
roc.fittedLR <- roc(testSet$msi, prob.fittedLR$MSS)