
library(caret)
library(xgboost)
library(e1071)

##Data transformation 
data_normalise <- function(df) {

  return(df/rowSums(df))
  
}




##Feature selection algorithm
species_select <-
  function(x,
           y,
           remove_correlated = T,
           subsets = NULL,
           cores = 1) {
    doParallel::registerDoParallel(cores)
    
    
    
    if (remove_correlated == T) {
      correlated <- findCorrelation(
        cor(x, use = "complete.obs"),
        cutoff = 0.98,
        verbose = FALSE,
        names = FALSE,
        exact = FALSE
      
      )
      
      x <- x[,-c(correlated)]
      print(paste0("correlated features removed:",length(correlated)))
      
    }
   
    len <- ncol(x)
    if(is.null(subsets)){
      subsets <-
        c(floor(len/2), floor(len / 4), floor(len / 8), floor(len / 16), floor(len / 32),floor(len / 64))
    }
    
    rfe_ctrl <- rfeControl(
      functions = rfFuncs,
      method = "cv",
      number =  5,
      verbose = FALSE,
      allowParallel = TRUE
    )
    set.seed(123)
    featureElimination <- rfe(
      x = x,
      y = y,
      sizes = subsets,
      rfeControl = rfe_ctrl,
      tuneLength = 2
    )
    doParallel::registerDoParallel(1) 
    return(featureElimination)
    
    
  }

##Main mGPS algorithm 
mGPS <-
  function(training = NULL,
           testing = NULL,
           classTarget,
           hierarchy = c('continent','city','latitude','longitude'),
           variables,
           nthread = 1,
           coast = NULL) {
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
        nrounds = c(300,600),
        eta = c( 0.05, 0.1),
        max_depth = c(3,6,9),
        gamma = 0,
        colsample_bytree = c(0.6,0.8),
        min_child_weight = c(1),
        subsample = (0.7)
      )
      
      if(length(hierarchy) == 4){
      
      Xgb_region <- train(x = training[,variables],y = training[,hierarchy[1]],
                          method = "xgbTree",
                          trControl = trControlClass,
                          tuneGrid = tune_grid,
                          nthread = nthread)
      
      l1_train <- data.frame(training[,c(variables)],Xgb_region[["pred"]][order(Xgb_region$pred$rowIndex),levels(training[,hierarchy[1]]) ])
      
      }else{
        
        l1_train <- training[,c(variables)]
        
      }
      Xgb_class <- train(x = l1_train,y = training[,classTarget],
                         method = "xgbTree",
                         trControl = trControlClass,
                         tuneGrid = tune_grid,
                         nthread = nthread)
      
      l2_train <- data.frame(l1_train,Xgb_class[["pred"]][order(Xgb_class$pred$rowIndex),levels(training[,classTarget]) ])
      
      
      Xgb_latitude <- train(x = l2_train ,y = training[,"latitude"],
                            method = "xgbTree",
                            trControl = trControl,
                            tuneGrid = tune_grid,
                            nthread = nthread)
      
      l3_train <- data.frame(l2_train, "latPred" = Xgb_latitude[["pred"]][order(Xgb_latitude$pred$rowIndex),"pred" ])
      
      Xgb_longitude <- train(x = l3_train ,y = training[,"longitude"],
                             method = "xgbTree",
                             trControl = trControl,
                             tuneGrid = tune_grid,
                             nthread = nthread)
      
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
      if(length(hierarchy) == 4){
      regProbs <- predict(Xgb_region, newdata = testing[,variables],type ="prob")
      
      l1_test <- data.frame(testing[,variables], regProbs)
      }else{
        l1_test <- testing[,variables]
      }
      classPred <- predict(Xgb_class, newdata = l1_test)
      classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
      
      l2_test <-  data.frame(l1_test, classProbs) 
      latPred <- predict(Xgb_latitude, newdata = l2_test)
      
      l3_test <- data.frame(l2_test, latPred)
      longPred <- predict(Xgb_longitude, newdata = l3_test)
      
      #adjust out of bounds predictions
      longPred[longPred > 180] <- 180
      longPred[longPred < -180] <- -180
      latPred[latPred > 90] <- 90
      latPred[latPred < -90] <- -90
      #Pull to nearest coastline if provided
      find_coast <- function(long, lat) {
        distances_from_coastline <-
          sp::spDistsN1(coast, c(long, lat), longlat = TRUE)
        
        closest_point <-  which.min(distances_from_coastline)
        new_coords <- coast[closest_point,]
        
        return(new_coords)
        
      }
      if (!is.null(coast)) {
        toAdjust <-
          which(is.na(maps::map.where(database = "world", longPred, latPred)))
        
        adjusted <-
          mapply(find_coast, long = longPred[toAdjust], lat = latPred[toAdjust])
        
        
        longPred[toAdjust] <- adjusted[1,]
        latPred[toAdjust] <- adjusted[2,]
        
        
        }
      
      
      
      return(list(classPred, latPred, longPred))
      
    }
    
  }



