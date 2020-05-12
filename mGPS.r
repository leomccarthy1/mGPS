
library(caret)
library(xgboost)
library(e1071)
mGPS <-  function(training = NULL,testing = NULL,classTarget,variables){
  
  #Check for training set
  
  if(is.null(training)){
    return(message("No training set given"))
    } else{


  training <- droplevels(training)
    
 #Train mGPS with 5-fold cross validation of training set for hyperparameter tuning. 
  message("Training mGPS...")

  set.seed(1234)
  folds <- createFolds(training[,classTarget], k = 5, returnTrain = T)
  
  
  
  trControlClass <-  trainControl(
    method = "cv",
    number = 5,  
    verboseIter = FALSE,
    returnData = FALSE,
    search = "grid",
    savePredictions = "final",
    classProbs = T,
    allowParallel = T,
    index = folds )
  
  
  
  trControl <-  trainControl(
    method = "cv",
    number = 5,  
    verboseIter = FALSE,
    returnData = FALSE,
    search = "grid",
    savePredictions = "final",
    allowParallel = T,
    index = folds)
  
  
  
  tune_grid <- expand.grid(
    nrounds = c(400,600),
    eta = c( 0.05, 0.1),
    max_depth = c(3,6,9),
    gamma = 0,
    colsample_bytree = c(0.6,0.8),
    min_child_weight = c(1),
    subsample = (0.7)
  )
  
  
  Xgb_region <- train(x = training[,variables],y = training[,"continent"],
                      method = "xgbTree",
                      trControl = trControlClass,
                      tuneGrid = tune_grid,
                      nthread = 1)
  
  l1_train <- data.frame(training[,c(variables)],Xgb_region[["pred"]][order(Xgb_region$pred$rowIndex),levels(training[,"continent"]) ])
  
  Xgb_class <- train(x = l1_train,y = training[,classTarget],
                     method = "xgbTree",
                     trControl = trControlClass,
                     tuneGrid = tune_grid,
                     nthread = 1)
  
  l2_train <- data.frame(l1_train,Xgb_class[["pred"]][order(Xgb_class$pred$rowIndex),levels(training[,classTarget]) ])
  
  
  Xgb_latitude <- train(x = l2_train ,y = training[,"latitude"],
                        method = "xgbTree",
                        trControl = trControl,
                        tuneGrid = tune_grid,
                        nthread = 1)
  
  l3_train <- data.frame(l2_train, "latPred" = Xgb_latitude[["pred"]][order(Xgb_latitude$pred$rowIndex),"pred" ])
  
  Xgb_longitude <- train(x = l3_train ,y = training[,"longitude"],
                         method = "xgbTree",
                         trControl = trControl,
                         tuneGrid = tune_grid,
                         nthread = 1)
  
}
  #check for test set, return trained model if no test set. 
  if(is.null(testing)){
    
    model <- function(test,variables){
      regProbs <- predict(Xgb_region, newdata = test[,variables],type ="prob")
      
      l1_test <- data.frame(test[,variables], regProbs)
      
      classPred <- predict(Xgb_class, newdata = l1_test)
      classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
      
      l2_test <-  data.frame(l1_test, classProbs) 
      latPred <- predict(Xgb_latitude, newdata = l2_test)
      
      l3_test <- data.frame(l2_test, latPred)
      longPred <- predict(Xgb_longitude, newdata = l3_test)
      return(list(classPred, latPred, longPred))
      
    }
    message("No test set...returning trained mGPS model function")
    return(list(Xgb_region,Xgb_class,Xgb_latitude,Xgb_longitude,"model" = model))
  }else{
  message("Generating predictions")
  #generate mGPS predictions for test set
  regProbs <- predict(Xgb_region, newdata = testing[,variables],type ="prob")
  
  l1_test <- data.frame(testing[,variables], regProbs)
  
  classPred <- predict(Xgb_class, newdata = l1_test)
  classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
  
  l2_test <-  data.frame(l1_test, classProbs) 
  latPred <- predict(Xgb_latitude, newdata = l2_test)
  
  l3_test <- data.frame(l2_test, latPred)
  longPred <- predict(Xgb_longitude, newdata = l3_test)
  
  
  
  
  
  
  return(list(classPred, latPred, longPred))
  
}
  
}








