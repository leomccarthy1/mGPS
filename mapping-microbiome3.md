This is a comprehensive outline of the main modeling methods used in our
peper titled \"\" in order to assist in reproducibility of results.

Data Preperation
================

Load useful packages

``` {.r}
library(sp)
library(rworldmap)
```

    ## ### Welcome to rworldmap ###

    ## For a short introduction type :   vignette('rworldmap')

``` {.r}
library(caret)
```

    ## Loading required package: lattice

    ## Loading required package: ggplot2

``` {.r}
library(maps)
library(MASS)
library(randomForest)
```

    ## randomForest 4.6-14

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

``` {.r}
library(geosphere)
library(doParallel)
```

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: parallel

``` {.r}
library(caret)
library(plyr)
```

    ## 
    ## Attaching package: 'plyr'

    ## The following object is masked from 'package:maps':
    ## 
    ##     ozone

``` {.r}
library(maptools)
```

    ## Checking rgeos availability: TRUE

``` {.r}
library(rgeos)
```

    ## rgeos version: 0.5-2, (SVN revision 621)
    ##  GEOS runtime version: 3.7.2-CAPI-1.11.2 
    ##  Linking to sp version: 1.3-1 
    ##  Polygon checking: TRUE

``` {.r}
library(mapplots)
```

First we import the metasub datasets for metadata (city, latitude,
longitude etc) and bacterial abundance data. Merge these sets by sample
ID.

Control samples are be removed for obvious reasons.

``` {.r}
#Import data sets 
complete_meta <- read.csv(file = "complete_metadatav2.csv", header = TRUE)
bac_data <- read.csv(file= "refseq.krakenhll_species (1).csv", header = TRUE)

#Change NA's to 0 in bacterial abundance data 
bac_data[is.na(bac_data)] <- 0 

#merge bacterial and meta data
metasub_data <- merge(complete_meta,bac_data,by.x="uuid",by.y="uuid")

#Remove control samples
control_samples <- c( which(metasub_data$city %in% c("control", "other_control","neg_control","other","pos_control")), which(metasub_data$control_type %in% c("ctrl cities","negative_control","positive_control"))) 
metasub_data <- droplevels(metasub_data[-c(control_samples), ])
```

We will also remove some samples for which the origin isn't clear
e.g. the name of the origin does match the co-ordinates. Samples from
bouroughs within london will be re-labelled simply as London for
consistency.

There are a group of samples from small towns accross uk, we will remove
these as there are very few from each town and they will clutter up
results A LOT while ading little improvement to a model on a global
scale

``` {.r}
#remove sparse samples locations and dubiously labelled samples. 
remove_samples <- which(metasub_data$city %in% c("tsuruoka", "antarctica","lands_end", "birmingham", "belfast", "bradford", "bristol","bury","cairngorm", "edinborough", "jaywick", "newcastle","swansea","eastbourne", "st_louis", "eden"))
metasub_data <- metasub_data[-c(remove_samples), ]

#re-label london boroughs 
metasub_data$city[metasub_data$city %in% c("kensington","islington")] <- "london" 
metasub_data <- droplevels(metasub_data)
```

``` {.r}
#Correction of identified misslabelling of data 
metasub_data$latitude[metasub_data$city == "kyiv"] <- metasub_data$city_latitude[metasub_data$city == "kyiv"]
metasub_data$longitude[metasub_data$city == "kyiv"] <- metasub_data$city_longitude[metasub_data$city == "kyiv"]
metasub_data$continent[metasub_data$city == "porto"] <- "europe"
```

Not all of the samples have exact origin co-ordinates recorded, for the
purposes of global predictions we can substitute the central
co-ordinates of the city of origin for the missing exact co-ordiantes,
found in the city meta-data. On the global scale this shouldn't be too
much of an issue.

``` {.r}
#impute unknown exact co-ords for city co-ords
metasub_data[is.na(metasub_data$latitude),]$latitude <- metasub_data[is.na(metasub_data$latitude),]$city_latitude
metasub_data[is.na(metasub_data$longitude),]$longitude <- metasub_data[is.na(metasub_data$longitude),]$city_longitude

#correction to some incorrect city co-ords for a few london samples
metasub_data[metasub_data$city == "london",]$city_latitude <- 51.50853
metasub_data[metasub_data$city == "london",]$city_longitude <- -0.12574
```

Feature selection
=================

This dataset contains rather a large amount of predictor variables, many
of which will be random noise and detrimental to both model accuracy and
computation time.

To select an optimum set of variables for each of the target variables
of interest (city,latitude and longitude) we will use recursive feature
elimination with random forests. The caret function "rfe" uses
re-sampling and external validation to protect against selecting
features that lead us to overfit the training data. More information is
avalable at
<https://topepo.github.io/caret/recursive-feature-elimination.html>

``` {.r}
set.seed(123)

#Activate parallel processing 
registerDoParallel(parallel::detectCores() - 1)

#recursive feature elimination for the target variable "city"
rfe_ctrl <- rfeControl(functions = rfFuncs,
                   method = "cv",
                   number =  5,
                   verbose = FALSE,
                   allowParallel = TRUE
                   )
featureElimination <- rfe(x = metasub_data[, c(43:3799)],y = metasub_data$city,
                 sizes = c(100,200,300,400,500),
                 rfeControl = rfe_ctrl,
                 tuneLength = 0
                 )

optimumVarsCity <- featureElimination$optVariables[1:200]
```

Here we can see the effect of selecting smaller subsets of features on
preict city of origin.

``` {.r}
plot(featureElimination, type = c("g", "o"))
```

![](mapping-microbiome3_files/figure-markdown/unnamed-chunk-7-1.png)

Looking at the results it's clear using all the variables is not
optimal. 300 variables is the optimum subset from these results, however
the accuracy loss by using only 200 variables is minisuclue (\<0.001%)
and using only 200 variables would speed up compuation so we will use
this as our optimum subset of bacterial species for further analysis
when aiming to predict city of origin.

We can take a look at the most 20 important variables, visually.

``` {.r}
imp <- varImp(featureElimination)
dotchart(rev(imp[1:25, "Overall"]),labels= rev(row.names(imp)[1:25]),cex=1,pt.cex = 1.2,
         xlab="Mean Decrease In Accuracy")
```

![](mapping-microbiome3_files/figure-markdown/unnamed-chunk-8-1.png)

Modeling
========

Our model makes use of regressor chains (more info:
<http://cig.fi.upm.es/articles/2015/Borchani-2015-WDMKD.pdf>) in order
to preoduce multi output predictions for class (city), latittude and
longitude. This techniques involves chaining predictions for a series of
outputs as input features for the next level model. In this case 3
chained XGboost models are used with hyperparameter tuning carried out
by 5-fold cross validation of the trainng set using grid search methods.

Here this function here takes 4 arguments: a training set, a test set, a
target class variable, a set of variables to use when training the
model. The function outputs the class,latitude and longitude predictions
from applying this model to the input test set.

``` {.r}
MappingModel <-  function(training,testing,classTarget,variables){
   
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
   
      
      
      Xgb_class <- train(x = training[,variables],y = training[,classTarget],
         method = "xgbTree",
         trControl = trControlClass,
         tuneGrid = tune_grid)
      
      l2_train <- data.frame(training[,c(variables)],Xgb_class[["pred"]][order(Xgb_class$pred$rowIndex),levels(training[,classTarget]) ])
      
      
            Xgb_latitude <- train(x = l2_train ,y = training[,"latitude"],
               method = "xgbTree",
               trControl = trControl,
               tuneGrid = tune_grid)
            
            l3_train <- data.frame(l2_train, "latPred" = Xgb_latitude[["pred"]][order(Xgb_latitude$pred$rowIndex),"pred" ])
            
                   Xgb_longitude <- train(x = l3_train ,y = training[,"longitude"],
                     method = "xgbTree",
                     trControl = trControl,
                     tuneGrid = tune_grid)
            
      
     
      
                  classPred <- predict(Xgb_class, newdata = testing[,variables])
                  classProbs <- predict(Xgb_class, newdata = testing[,variables],type ="prob")
                  
                        l2_test <-  data.frame(testing[,variables], classProbs) 
                           latPred <- predict(Xgb_latitude, newdata = l2_test)
                     
                           l3_test <- data.frame(l2_test, latPred)
                              longPred <- predict(Xgb_longitude, newdata = l3_test)
                  
                  
                 
     
         
   
return(list(classPred, latPred, longPred))

}
```

Using our algorithm will will employ strativied 10-fold validation on
the metasub data set in order to generate test predictions for city,
latitude and longitude.

``` {.r}
#generate 5 stratified folds for test predictions.
set.seed(18)
trainFolds <-  createFolds(metasub_data$city, k = 5, returnTrain = T)


GeoPreds <- list()
registerDoParallel(7) 
#iteratively train the model on each of the 10 training folds and generate predictions using the coresponding test fold.
for (i in 1:5){

train <- metasub_data[trainFolds[[i]],]
test <- metasub_data[-trainFolds[[i]],]
   
testPreds <-MappingModel(training = train, testing = test, classTarget = "city",variables = optimumVarsCity)
GeoPreds[[i]] <- testPreds
   
}


#Combine these test predictions into one data set 
add_preds <- list()
for (i in 1:5){
   
   add_preds[[i]] <- cbind(metasub_data[-trainFolds[[i]],] , 
                                 "cityPred"= GeoPreds[[i]][[1]], 
                                 "latPred" = GeoPreds[[i]][[2]], 
                                 "longPred" = GeoPreds[[i]][[3]] )
         
   
   
}

MetasubDataPreds <- rbind.fill(add_preds)
MetasubDataPreds[MetasubDataPreds$longPred > 180,"longPred"] <- 180.000
```

The final stage of this model is to adjust any predicted co-ordiantes
that don't lie on land. This step makes the assumption that all future
tetsing data will be taken from land we beleive this is a valid approach
to take.

``` {.r}
##Final stage is to adjust any co-ordinates that are in the sea and pull them to the nearset land mass

#get world coastlines
coastlines <- cbind("x"  = SpatialLines2map(coastsCoarse)$x ,"y" =SpatialLines2map(coastsCoarse)$y)
   coastlines <- coastlines[complete.cases(coastlines),]
      coastlines <- coastlines[coastlines[,1] < 180 ,]


#Function to find the nearest land point of any GPS co-ordinates
find_coast <- function(long,lat){
    distances_from_coastline <-  spDistsN1(coastlines , c(long,lat), longlat = TRUE)
    
                                    closest_point <-  which.min(distances_from_coastline)
                                    new_coords <- coastlines[closest_point,]
    
                                       return(new_coords)
 
                                     }

      #Find points than need adjusting 
      GPS_where <- map.where(database = "world", MetasubDataPreds$longPred, MetasubDataPreds$latPred)
         toAdjust <- MetasubDataPreds[which(is.na(GPS_where)),]
         adjusted <- mapply(find_coast, long = toAdjust$longPred, lat = toAdjust$latPred )
      #Adjust points 
            MetasubDataPreds[which(is.na(GPS_where)), "latPred"] <- adjusted[2,]
            MetasubDataPreds[which(is.na(GPS_where)), "longPred"] <- adjusted[1,]
```

``` {.r}
#Print test results 
print(c(mean(MetasubDataPreds$cityPred ==MetasubDataPreds$city)
   ,RMSE(MetasubDataPreds$latPred, MetasubDataPreds$latitude),
       RMSE(MetasubDataPreds$longPred, MetasubDataPreds$longitude)))
```

    ## [1]  0.9262899  8.2106433 23.4692890

For 10-fold test predictions we have an accuracy of 93% for city
predictions and RMSE of 8.2 and 23.5 for latitude and longitude
respectively.

Now we have our predicions we can begin to visulize them in order to get
a feel for the effectiveness of the algorithm, as in this context RMSE
is hard to understand paractially. So we will use predicted co-ordinates
to a generate the distance of predictions from the true origin.

``` {.r}
#Using lat long predictions determine distance (km) of the prediction from true origin using haversine distance. 

for (i in 1:nrow(MetasubDataPreds)){
  MetasubDataPreds[i,"Distance_from_origin"] <- distm(c(MetasubDataPreds[i,"longPred"],MetasubDataPreds[i,"latPred"]), c(MetasubDataPreds[i,"longitude"],MetasubDataPreds[i,"latitude"]), fun = distHaversine)/1000
}


#Print distance from origin results 
print(c(mean(MetasubDataPreds$Distance_from_origin ),
median(MetasubDataPreds$Distance_from_origin ),
mean(MetasubDataPreds$Distance_from_origin < 500)))
```

    ## [1] 1013.4403826  256.7242975    0.6402948

The mean distance of the origin prediction from the tru origin is 1013km
and median 256.7 with 64% of samples within 500km

We can plot the predicted origin of each global sample. Coloured by
continent of origin. City prediction accuracy is showm by the pie
charts.

``` {.r}
#####world map showing by continent
map <- getMap(resolution = "coarse")

palette <-c( "darkorchid4","gold2","dodgerblue3","brown","orangered2","mediumspringgreen","deeppink2")
plot(map, xlim = c(-165,168),ylim = c(-70,60), col = "grey",border = "darkgrey", xlab = "Longitude", ylab = 'Latitude', bg = "lightskyblue1")

#find coord preds by region
for ( i in 1:length(levels(MetasubDataPreds$continent))){
  this_continent <- levels(MetasubDataPreds$continent)[i]
  find_lats <- MetasubDataPreds[MetasubDataPreds[,"continent"] == this_continent,][,"latPred"]
  find_longs <- MetasubDataPreds[MetasubDataPreds[,"continent"] == this_continent,][,"longPred"]
  
  #plot predicted co-ordinates
  points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.2)
  
  #plot city prediction accuravy by continent as pies
  correctly_pred <-  mean(MetasubDataPreds[MetasubDataPreds$continent == this_continent,"cityPred"]== 
                                 MetasubDataPreds[MetasubDataPreds$continent == this_continent,"city"]) 
  incorrectly_pred <- (1 - correctly_pred) 

  
  
  
  continent_lats <- c(55,69,8,40,-40,-10,-5)
  continent_longs <- c(125,0,60,-130,140,-80,5)
  
  add.pie(z = c(correctly_pred, incorrectly_pred), x = continent_longs[i], y = continent_lats[i]
             ,edges=200,
             radius=10,
             col=c(palette[i],"black") , labels = ""
  )
}

#Plot city sampling locations
map.axes(cex.axis = 0.8)
par(fig = c(0,0.4,0.0,0.5), new = T) 
plot(map,xlim = c(-165,168), ylim = c(-70,65), col = "grey", border = "darkgrey", bg ="lightskyblue1")
for ( i in 1:length(levels(MetasubDataPreds$continent))){
  this_continent <- levels(MetasubDataPreds$continent)[i]
  find_lats <- MetasubDataPreds[MetasubDataPreds$continent == this_continent,]$city_latitude
  find_longs <- MetasubDataPreds[MetasubDataPreds$continent == this_continent,]$city_longitude
  
  points(find_longs, find_lats, col = palette[i], pch = 17, cex = 1)
}

legend(-165,-15, c("East Asia","Eurpoe","Middle East",
                 "North America",
                 "Oceania",
                 "South America",
                 "Sub Saharan Afica"), pch = 17, col = palette, cex = 0.5, bg ="lightskyblue1")
box( col = 'black')
```

![](mapping-microbiome3_files/figure-markdown/unnamed-chunk-14-1.png)

Further code for plots of the predictions seperated by region can be
found in the file "global\_plots.R"

\#Predictions unseen cities

Here we will test our model on predicting samples from cities that are
not seen in the training data. For this stage we ommit the hyper
paramter tuning stage in order to save time with computation and only a
marginal loss in accuracy. An effective standard set of hyperparameters
are assigned instead.

``` {.r}
#Adjust model slightly to not inc lude hyperparameter tuning. 
MappingModel_fixedParams <-  function(training,testing,classTarget,variables){
   
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
        nrounds = c(400),
        eta = c( 0.05),
        max_depth = c(6),
        gamma = 0,
        colsample_bytree = c(0.6),
        min_child_weight = c(1),
        subsample = (0.7)
      )
   
      
      
      Xgb_class <- train(x = training[,variables],y = training[,classTarget],
         method = "xgbTree",
         trControl = trControlClass,
         tuneGrid = tune_grid)
      
      l2_train <- data.frame(training[,c(variables)],Xgb_class[["pred"]][order(Xgb_class$pred$rowIndex),levels(training[,classTarget]) ])
      
      
            Xgb_latitude <- train(x = l2_train ,y = training[,"latitude"],
               method = "xgbTree",
               trControl = trControl,
               tuneGrid = tune_grid)
            
            l3_train <- data.frame(l2_train, "latPred" = Xgb_latitude[["pred"]][order(Xgb_latitude$pred$rowIndex),"pred" ])
            
                   Xgb_longitude <- train(x = l3_train ,y = training[,"longitude"],
                     method = "xgbTree",
                     trControl = trControl,
                     tuneGrid = tune_grid)
            
      
     
      
                  classPred <- predict(Xgb_class, newdata = testing[,variables])
                  classProbs <- predict(Xgb_class, newdata = testing[,variables],type ="prob")
                  
                        l2_test <-  data.frame(testing[,variables], classProbs) 
                           latPred <- predict(Xgb_latitude, newdata = l2_test)
                     
                           l3_test <- data.frame(l2_test, latPred)
                              longPred <- predict(Xgb_longitude, newdata = l3_test)
                  
                  
                 
     
         
   
return(list(classPred, latPred, longPred))

}
```

First randomly select 10 cities from the data set. Then create 10
training folds that contain all data other than the data for each of
these cities in turn. Test the trained model by generating predictions
for samples from the city that was left out, do this for all 10 cities.

``` {.r}
#generate folds for test leave-out validation at city level predictions.
set.seed(9)
leaveOutCitys <- sample (levels(metasub_data$city), size = 10, replace =F)


LeaveOutFolds <- list()
for(i in 1:length(levels(metasub_data$city))){
   
   LeaveOutFolds[[i]] <- which(metasub_data$city == levels(metasub_data$city)[i])
   
}

 LeaveOutFolds 

GeoPredsLeaveOut <- list()
registerDoParallel(7) 
#iteratively train the model on each of the 10 training folds and generate predictions using the coresponding test fold.
for (i in 1:length(levels(metasub_data$city)) ){

train <- droplevels(metasub_data[-LeaveOutFolds[[i]],])
test <- metasub_data[LeaveOutFolds[[i]],]
   
testPreds <-MappingModel_fixedParams(training = train, testing = test, classTarget = "city",variables = optimumVarsCity)
GeoPredsLeaveOut[[i]] <- testPreds
   
}



GeoPredsLeaveOut
add_preds <- list()
for (i in 1:10){
   
   add_preds[[i]] <- cbind(metasub_data[LeaveOutFolds[[i]],] , 
                                 "cityPred"= GeoPredsLeaveOut[[i]][[1]], 
                                 "latPred" = GeoPredsLeaveOut[[i]][[2]], 
                                 "longPred" = GeoPredsLeaveOut[[i]][[3]] )
         
   
   
}

MetasubLeaveOutPreds <- droplevels(rbind.fill(add_preds))
```

Adjust to nearest coast as before

``` {.r}
GPS_where_lo <- map.where(database = "world", MetasubLeaveOutPreds$longPred, MetasubLeaveOutPreds$latPred)
         toAdjustLo <- MetasubLeaveOutPreds[which(is.na(GPS_where_lo)),]
         adjustedLo <- mapply(find_coast, long = toAdjustLo$longPred, lat = toAdjustLo$latPred )
      
            MetasubLeaveOutPreds[which(is.na(GPS_where_lo)), "latPred"] <- adjustedLo[2,]
            MetasubLeaveOutPreds[which(is.na(GPS_where_lo)), "longPred"] <- adjustedLo[1,]
```

Now find the distance from origin for these preditions.

``` {.r}
for (i in 1:nrow(MetasubLeaveOutPreds)){
  MetasubLeaveOutPreds[i,"Distance_from_origin"] <- distm(c(MetasubLeaveOutPreds[i,"longPred"],MetasubLeaveOutPreds[i,"latPred"]), c(MetasubLeaveOutPreds[i,"longitude"],MetasubLeaveOutPreds[i,"latitude"]), fun = distHaversine)/1000
}

print(c(mean(MetasubLeaveOutPreds$Distance_from_origin ),
median(MetasubLeaveOutPreds$Distance_from_origin ),
mean(MetasubLeaveOutPreds$Distance_from_origin < 1250) ) )
```

    ## [1] 3658.2795603 2661.1270019    0.2650794
