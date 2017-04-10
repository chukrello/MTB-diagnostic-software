library(caret)
library(glmnet)
library(xgboost)
library(plyr)
getwd()

df <- read.csv('regr_input_PZA.txt', header=F, sep='\t', na.strings=c('No result','M', 'F'))
df <- fread('regr_input_small.txt')
dim(df)
str(df[,330])
names(df)[330] <- 'INH'
df <- df[complete.cases(df),]
df <- df[1:1800,]
rownames(df) <- df[,1]
df[,1] <- NULL
apply(df[1:400,1:101], 2,function(x) chisq.test(table(df[,x], df$IZN))$p.value)


##Параметры тренинга модели

fitControl <- trainControl(## 7-fold CV
  method = "repeatedcv",
  number = 8,
  repeats = 15,
  classProbs = TRUE,
  savePredictions = T,
  summaryFunction = twoClassSummary)

##Разделение на тренинговую и тестовую выборку

inTraining <- createDataPartition(df$V330, p = .75, list = FALSE)
training <- df[ inTraining,]
testing  <- df[-inTraining,]

##Использование алгоритма. Random forest
set.seed(900)

rfmodel_1 <- train(V330~., data=training, method='rf', trControl=fitControl)
rfmodel_1$results
predict(rfmodel_1, newdata = testing, type='prob')

predictions <- predict(rfmodel_1, testing)
confusionMatrix(predictions, testing$IZN)

importance <- varImp(rfmodel_1, scale=FALSE)
importance


##Алгоритм - glmnet

glmnet_1 <- train(V330~.,data=training, method='glmnet', trControl=fitControl,family='binomial')

glmnet_1$results
predict(glmnet_1, newdata=testing, type='prob')
glmnet_1$finalModel
predictions <- predict(glmnet_1, testing)
confusionMatrix(predictions, testing$INH)

predictors(glmnet_1)
varImp(glmnet_1)
##Использование алгоритма - glm

glmodel_1 <- train(training$IZN ~.,data=training, method='glm', family='binomial',trControl=fitControl)

predictions=predict(glmodel_1, testing)
confusionMatrix(predictions, testing$Group)
predictors(glmodel_1)

##Использование алгоритма - xgboost

xgmodel_1 <- train(IZN~., data=training, method='xgbTree', trControl=fitControl)

xgmodel$results
predictions <- predict(xgmodel_1, testing)
confusionMatrix(predictions, testing$Drug)
predictors

## А теперь визуализируем

importance <- varImp(glmnet_1, scale=FALSE)

##Обобщённая функция для вычисления

make_prediction <- function(df, Drug, algorithm){
  fitControl <- trainControl(## 7-fold CV
    method = "repeatedcv",
    number = 8,
    repeats = 15,
    classProbs = TRUE,
    savePredictions = T,
    summaryFunction = twoClassSummary)

  inTraining <- createDataPartition(df1$Drug, p = .75, list = FALSE)
  training <- df[ inTraining,]
  testing  <- df[-inTraining,]
  
  resul <<- train(Drug~., data=df, trControl=fitControl, method=algorithm)
  predictions <- predict(resul, testing)
  best <<- predictors()
  return(confusionMatrix(predictions))}