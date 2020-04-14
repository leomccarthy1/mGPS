

library(nominatim)
library(mapdataNE)
library(dplyr)
library(sp)
library(rworldmap)
library(caret)
library(rpart)
library(maps)
library(MASS)
library(randomForest)
library(geosphere)
library(doParallel)

environmental <- read.csv(file="Dataset_01_22_2018_enviro.csv",header=TRUE)
taxa <- read.csv(file="Dataset_01_22_2018_taxa.csv",header=TRUE)
SoilData <- merge(environmental,taxa,by.x="ID_Environmental",by.y="ID_Environmental")

#get country of each sample

SoilData$Country <- map.where(database="world", SoilData$Longitude, SoilData$Latitude)
SoilData[234:237,]$Country <- "Chile"
SoilData[c(2,6,9,15,22,67,34),]$Country <- "USA"

SoilData$Country <- factor(make.names(SoilData$Country))
colnames(SoilData)[3:4] <- c("longitude","latitude")

#Remove countries with insufficient samples
keep <- levels(SoilData$Country)[table(SoilData$Country) > 2]
SoilData <- droplevels(SoilData[SoilData$Country%in% keep, ])




###modeling ###

#Remove redundant features
registerDoParallel(parallel::detectCores() - 1)
rfe_ctrl <- rfeControl(functions = rfFuncs,
                       method = "cv",
                       number = 5,
                       verbose = FALSE,
                       allowParallel = T
)
featureEliminationSoil <- rfe(x = SoilData[,c(28:538)],y = SoilData$Country,
                              sizes = c(20,30,50,100,200),
                              rfeControl = rfe_ctrl,
                              tuneLength = 5
)


optimumVarsSoil <- featureEliminationSoil$optVariables


#Forest of optimal species subset
fullForest <- randomForest(x = SoilData[, optimumVarsSoil], y = SoilData$Country,importance = T)
#Extract and rank species importance
v <- varImp(fullForest, type = 1, scale = F)
species <- row.names(v[order(v$Overall,decreasing = T),0])
value <- v[order(v$Overall,decreasing = T),1]
#plot
png("ImpOTU.png", width = 12,height = 8, units = 'in', res = 600)
dotchart(rev(value[1:25]),labels= rev(species[1:25]),cex=1.2,pt.cex = 1.3,
         xlab="Mean decrease in accuracy", mgp = c(2.2,0,0))
dev.off()

write.csv(species,"soil_200_OTU.csv")

varImpPlot(fullForest, type = 1, scale = F)
# define model
mGPS <-  function(training,testing,classTarget,variables){
  
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
  
  
  Xgb_region <- train(x = training[,variables],y = training[,"Continent"],
                      method = "xgbTree",
                      trControl = trControlClass,
                      tuneGrid = tune_grid,
                      nthread = 1)
  
  l1_train <- data.frame(training[,c(variables)],Xgb_region[["pred"]][order(Xgb_region$pred$rowIndex),levels(training[,"Continent"]) ])
  
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


#5 fold cross validation
set.seed(18)
trainFolds <-  createFolds(SoilData$Country, k = 5, returnTrain = T)

GeoPredsSoil <- list()
registerDoParallel(7) 

for (i in 1:5){
  
  train <- SoilData[trainFolds[[i]],]
  test <- SoilData[-trainFolds[[i]],]
  
  testPreds <-mGPS(training = train, testing = test, classTarget = "Country",variables = optimumVarsSoil)
  GeoPredsSoil[[i]] <- testPreds
  
}



add_preds <- list()
for (i in 1:5){
  
  add_preds[[i]] <- cbind(SoilData[-trainFolds[[i]],] , 
                          "countryPred"= GeoPredsSoil[[i]][[1]], 
                          "latPred" = GeoPredsSoil[[i]][[2]], 
                          "longPred" = GeoPredsSoil[[i]][[3]] )
  
  
  
}

SoilDataPreds <- rbind.fill(add_preds)




#coast ajust#
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

GPS_where <- map.where(database = "world", SoilDataPreds$longPred, SoilDataPreds$latPred)
toAdjust <- SoilDataPreds[which(is.na(GPS_where)),]
adjusted <- mapply(find_coast, long = toAdjust$longPred, lat = toAdjust$latPred )

SoilDataPreds[which(is.na(GPS_where)), "latPred"] <- adjusted[2,]
SoilDataPreds[which(is.na(GPS_where)), "longPred"] <- adjusted[1,]




#plotting#
SoilData$Country
###global prediction points plot####

palette <-c( "deeppink2","darkorchid4","orangered2","gold2","brown","mediumspringgreen")
map <- getMap(resolution = "coarse")


png("SoilGlobal.png", width = 13,height = 8, units = 'in', res = 600)
plot(map,xlim = c(-160,160), col = "grey", border = "darkgrey", bg = "lightskyblue1", xlab = "", ylab = "")
title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
for (i in 1:length(levels(SoilDataPreds$Continent))){{
  this_continent <- levels(SoilDataPreds$Continent)[i]
  find_lats <- SoilDataPreds[SoilDataPreds[,"Continent"] == this_continent,]$latPred
  find_longs <- SoilDataPreds[SoilDataPreds[,"Continent"] == this_continent,]$longPred
  points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.4)
  }
legend(-160,0,legend=c("Africa", "Asia", "Oceania", "Europe", "North America","South America"),
         col=palette,pch = "+",cex=1.4, bg = "lightskyblue1")  
}
map.axes()
dev.off()







###distance from country of origin bar plot###

SoilDataPreds$GPSassign <- map.where(database = "world", SoilDataPreds$longPred, SoilDataPreds$latPred)
SoilDataPreds[is.na(SoilDataPreds$GPSassign),"GPSassign"] <- "In_sea"

for (i in 1:length((SoilDataPreds$Country))){ 
  country_lats <- map(database ="world",regions = SoilDataPreds[i,"Country"], plot = FALSE)$y
  country_longs <- map(database ="world",regions = SoilDataPreds[i,"Country"], plot = FALSE)$x
  pred_lat <- SoilDataPreds[i,]$latPred
  pred_long  <- SoilDataPreds[i,]$longPred
  distance <- c()
  if(SoilDataPreds[i,]$GPSassign == SoilDataPreds[i,]$Country){
    SoilDataPreds[i,"distance_from_country"] <- 0
  } else{
    for (n in 1:length(country_lats)){
      distance[n] <- distHaversine(c(pred_long ,pred_lat ),c(country_longs[n],country_lats[n]))/1000
      
    }
    SoilDataPreds[i,"distance_from_country"] <- min(distance, na.rm = TRUE)
  }
}


bar_df <- data.frame(row.names = c( "Overall",levels(SoilDataPreds$Country)))

for (i in 1: length(levels(SoilDataPreds$Country))){
  overall_prop <- mean(SoilDataPreds[,"distance_from_country"] < 500)
  bar_df[1,"0 - 500km"] <- overall_prop
  
  this_country <- levels(SoilDataPreds$Country)[i]
  prop <- mean(SoilDataPreds[SoilDataPreds$Country == this_country,][,"distance_from_country"] < 500)
  bar_df[i+1,"0 - 500km"] <- prop
}

for (i in 1: length(levels(SoilDataPreds$Country))){
  this_country <- levels(SoilDataPreds$Country)[i]
  prop <- mean(SoilDataPreds[SoilDataPreds$Country == this_country,][,"distance_from_country"] > 500 & SoilDataPreds[SoilDataPreds$Country == this_country,][,"distance_from_country"] < 1000)
  bar_df[i+1,"500 - 1000km"] <- prop
  
  overall_prop <- mean(SoilDataPreds[,"distance_from_country"] > 500 & SoilDataPreds[,"distance_from_country"] < 1000)
  bar_df[ 1,"500 - 1000km"] <- overall_prop
}

for (i in 1: length(levels(SoilDataPreds$Country))){
  this_country<- levels(SoilDataPreds$Country)[i]
  prop <- mean(SoilDataPreds[SoilDataPreds$Country == this_country,][,"distance_from_country"] > 1000 & SoilDataPreds[SoilDataPreds$Country == this_country,][,"distance_from_country"] < 2000)
  bar_df[i+1,"1000 - 2000km"] <- prop
  
  overall_prop <- mean(SoilDataPreds[,"distance_from_country"] > 1000 & SoilDataPreds[,"distance_from_country"] < 2000)
  bar_df[ 1,"1000 - 2000km"] <- overall_prop
}
for (i in 1: length(levels(SoilDataPreds$Country))){
  this_country <- levels(SoilDataPreds$Country)[i]
  prop <- mean(SoilDataPreds[SoilDataPreds$Country == this_country,][,"distance_from_country"] > 2000 & SoilDataPreds[SoilDataPreds$Country == this_country,][,"distance_from_country"] < 3000)
  bar_df[i+1,"2000 - 3000km"] <- prop
  
  overall_prop <- mean(SoilDataPreds[,"distance_from_country"] > 2000 & SoilDataPreds[,"distance_from_country"] < 3000)
  bar_df[ 1,"2000 - 3000km"] <- overall_prop
}
for (i in 1: length(levels(SoilDataPreds$Country))){
  this_country <- levels(SoilDataPreds$Country)[i]
  prop <- mean(SoilDataPreds[SoilDataPreds$Country == this_country,][,"distance_from_country"] > 3000 )
  bar_df[i+1,"> 3000km"] <- prop
  
  overall_prop <- mean(SoilDataPreds[,"distance_from_country"] > 3000)
  bar_df[ 1,"> 3000km"] <- overall_prop
}
size1 <- c()
for (i in 1: length(levels(SoilDataPreds$Country))){
  
  this_country <- levels(SoilDataPreds$Country)[i]
  size1[i] <- length(which(SoilDataPreds$Country == this_country))
}




par(xpd = T, mar = par()$mar + c(1,0,0,7))
bp <- barplot(t(bar_df*100), col=c("slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), names.arg=c("Overall",paste0(levels(SoilDataPreds$Country),"  (",size1,")")) ,args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.6,ylab = "Proportion of sample predictions %")
legend("topright",inset = c(-0.15,0.4), rev(c(colnames(bar_df))), fill = rev(c("slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 0.8)

par(mar=c(5, 4, 4, 2) + 0.1)

png("Soilbar.png", width = 13,height = 8, units = 'in', res = 600)
par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(bar_df*100), col=c("slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), 
              names.arg=c("Overall",paste0(levels(SoilDataPreds$Country),"  (",size1,")")) ,
              args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, 
              cex.names=.6,ylab = "", axisnames = F,axes = F, space =0)

axis(side =2, pos = 0)
mtext(text = c("Overall",paste0(levels(SoilDataPreds$Country),"  (",size1,")")) , side = 1, at = bp, line = 0, padj = 1, cex = 0.7)
title(ylab="Proportion of sample predictions %", mgp=c(0,0,0),cex.lab=1)
legend("topright",inset = c(-0.1,0.4), rev(c(colnames(bar_df))), fill = rev(c("slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 0.8)

par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()


