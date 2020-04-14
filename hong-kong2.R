  
load("~/Documents/Microbiome mapping/Markdown/GlobalMappingObj.RData")
library(sp)
library(rworldmap)
library(caret)
library(maps)
library(MASS)
library(randomForest)
library(geosphere)
library(DiagrammeR)
library(doParallel)
library(caret)
library(rgdal)
library(rgdal)
#Import data sets 
complete_meta <- read.csv(file = "complete_metadata.csv", header = TRUE)
bac_data <- read.csv(file= "refseq.krakenhll_species.csv", header = TRUE)
bac_data[is.na(bac_data)] <- 0 
metasub_data <- merge(complete_meta,bac_data,by.x="uuid",by.y="uuid")
control_samples <- c( which(metasub_data$city %in% c("control", "other_control","neg_control","other","pos_control")), which(metasub_data$control_type %in% c("ctrl cities","negative_control","positive_control"))) 
metasub_data <- droplevels(metasub_data[-c(control_samples), ])


#Extract the samples from hong kong for which the station of origin is known
  hong_kong <- hong_kong <- metasub_data[metasub_data$city == "hong_kong",]
  hong_kong <- hong_kong[-which(hong_kong$station == ""), ]
  
  keep <- levels(hong_kong$station)[table(hong_kong$station) > 2]
  hong_kong <- droplevels(hong_kong[hong_kong$station %in% keep, ] )


#variable selection for this dataset to remove redundant variables

  set.seed(123)
  
  
  registerDoParallel(parallel::detectCores() - 1)
  
  
  rfe_ctrl <- rfeControl(functions = rfFuncs,
                         method = "cv",
                         number =  5,
                         verbose = FALSE,
                         allowParallel = TRUE
  )
  

  featureEliminationHk <- rfe(x = hong_kong[, c(43:3799)],y = hong_kong$station,
                            sizes = c(100,200,300,400,500),
                            rfeControl = rfe_ctrl,
                            tuneLength = 1
  )
  
  featureEliminationHk
  
  optimumVarsHk <- featureEliminationHk$optVariables


  
#Define mapping model
  
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
  

#Apply our mapping model algorithm
  
  #generate 5 stratified folds for test predictions.
  set.seed(18)
  trainFoldsHk <-  createFolds(hong_kong$station, k = 5, returnTrain = T)
  
  
  GeoPredsHk <- list()
  registerDoParallel(parallel::detectCores() - 1) 
  #iteratively train the model on each of the 10 training folds and generate predictions using the coresponding test fold.
  for (i in 1:5){
    
    train <- hong_kong[trainFoldsHk[[i]],]
    test <- hong_kong[-trainFoldsHk[[i]],]
    
    testPreds <-mGPS(training = train, testing = test, classTarget = "station",variables = optimumVarsHk)
    GeoPredsHk[[i]] <- testPreds
    
  }
 





  add_preds <- list()
  for (i in 1:5){
    
    add_preds[[i]] <- cbind(hong_kong[-trainFoldsHk[[i]],] , 
                            "stationPred"= GeoPredsHk[[i]][[1]], 
                            "latPred" = GeoPredsHk[[i]][[2]], 
                            "longPred" = GeoPredsHk[[i]][[3]] )
    
    
    
  }
  
HkDataPreds <- rbind.fill(add_preds)


confusionMatrix(HkDataPreds$stationPred, HkDataPreds$station)

##Final stage is to adjust any co-ordinates that are in the sea and pull them to the nearset land mass

#get HK coastlines 
HKshp <- readOGR(dsn ="/Users/leomccarthy/Documents/Microbiome mapping/Boundaries and Rail/gadm36_HKG_shp", layer = "gadm36_HKG_1")
HKcoords <- cbind(SpatialPolygons2map(HKshp)$x, SpatialPolygons2map(HKshp)$y)
  HKcoords <- HKcoords[complete.cases(HKcoords),]

#Function to find the nearest land point of any GPS co-ordinates
find_coastHk <- function(long,lat){
  distances_from_coastline <-  spDistsN1(HKcoords , c(long,lat), longlat = TRUE)
  
  closest_point <-  which.min(distances_from_coastline)
  new_coords <- HKcoords[closest_point,]
  
  return(new_coords)
  
}

#Adjust for in-sea predictions
GPS_where <- map.where(database = "world", HkDataPreds$longPred,HkDataPreds$latPred)
toAdjust <- HkDataPreds[which(is.na(GPS_where)),]
adjusted <- mapply(find_coastHk, long = toAdjust$longPred, lat = toAdjust$latPred )

HkDataPreds[which(is.na(GPS_where)), "latPred"] <- adjusted[2,]
HkDataPreds[which(is.na(GPS_where)), "longPred"] <- adjusted[1,]





##Plots
#Get hk rail line


library(randomcoloR)
#color palette 
n <- length(levels(HkDataPreds$station))
palette <- randomColor(count = n, luminosity = "bright")

#plot predicted GPS locations
png("Hong_kong_test.png", width = 12, height = 8, units = 'in', res = 600)
plot(HKshp,xlim = c(114,114.2), ylim = c(22.17,22.47), col = "grey",border = "darkgrey", bg = "lightskyblue1", xlab = "", ylab = "")
title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
for (i in 1:length(levels(HkDataPreds$station))){
  this_station<- levels(HkDataPreds$station)[i]
  find_lats <- HkDataPreds[HkDataPreds[,"station"] == this_station,]$latPred
  find_longs <- HkDataPreds[HkDataPreds[,"station"] == this_station,]$longPred
  points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.5)
  
}
#plot station locations
for (i in 1:length(levels(HkDataPreds$station))){
  this_station<- levels(HkDataPreds$station)[i]
  
  station_lat <- mean(HkDataPreds[HkDataPreds$station == this_station,]$latitude)
  station_long <- mean(HkDataPreds[HkDataPreds$station == this_station,]$longitude)
  
  points(station_long, station_lat, col = "black", bg =palette[i] ,pch = 24, cex = 1.3)
}
#legend 
legend(114.3,22.45,legend=c(levels(HkDataPreds$station)), col = palette,pch = 17, cex = 0.8, bg = "lightskyblue1")
map.axes(cex.axis = 0.8)
map.scale( cex = 1)
#inset map 
par(fig = c(0,0.3,0.5,1), new = T) 
plot(map,xlim = c(110,116), ylim = c(22.15,22.45), col = "grey", bg ="lightskyblue1",border = "darkgrey")
points(HkDataPreds$city_longitude, HkDataPreds$city_latitude,pch = 16)
text(114.175,22.05, labels = "Hong Kong", pch = 2, cex = 0.8, font = 2)
rect(113.15, 21.4 ,115.2, 23.1, density = NULL, angle = 45,
     col = NA, border = NULL, pch = 2)
box( col = 'black')
dev.off()




## Map 2 - Pies 

find_island <- function(long,lat){
  ### generate west island coords
  #west_lats <- SpatialPolygons2map(HKshp[3,])$y
  #west_longs <- SpatialPolygons2map(HKshp[3,])$x
  
  north_lats <- SpatialPolygons2map(HKshp[c(4:10,12:14,16:18),])$y
  north_longs <- SpatialPolygons2map(HKshp[c(4:10,12:14,16:18),])$x
  
  south_lats <- SpatialPolygons2map(HKshp[c(11,1:2,15),])$y
  south_longs <- SpatialPolygons2map(HKshp[c(11,1:2,15),])$x
  
  west_lats <- SpatialPolygons2map(HKshp[3,])$y
  west_longs <- SpatialPolygons2map(HKshp[3,])$x
  
  if(lat > 22.4){
    island <- "North"
  }
  
  else if(point.in.polygon(long,lat,north_longs,north_lats) %in% c(1,2,3)){
    island <- "North"
  }
  else if(point.in.polygon(long,lat,south_longs,south_lats) %in% c(1,2,3)){
    island <- "South"
  }
  else if(point.in.polygon(long,lat,west_longs,west_lats) %in% c(1,2,3)){
    island <- "West"
  }
  else {
    island <- "None"
  }
  return(island)
}

HkDataPreds$Island <- factor(mapply(find_island, long = HkDataPreds$longitude, lat = HkDataPreds$latitude))
HkDataPreds$predIsland <- factor(mapply(find_island, long = HkDataPreds$longPred, lat =HkDataPreds$latPred))

library(mapplots)

png("Hong_kong_pie.png", width = 12, height = 8, units = 'in', res = 600)
plot(HKshp,xlim = c(114.15,114.15), ylim = c(22.17,22.47), col = "grey", border = "darkgrey",bg = "lightskyblue1", xlab = "",ylab = "")
title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
for (i in 1:length(levels(HkDataPreds$station))){
  this_station<- levels(HkDataPreds$station)[i]
  
  station_lat <- mean(HkDataPreds[HkDataPreds$station == this_station,]$latitude)
  station_long <- mean(HkDataPreds[HkDataPreds$station == this_station,]$longitude)
  
  
  correctly_predicted <- mean(HkDataPreds[HkDataPreds$station == this_station,]$stationPred == this_station ) 
  correct_island <- mean(HkDataPreds[HkDataPreds$station == this_station,]$predIsland == HkDataPreds[HkDataPreds$station == this_station,]$Island
                         & HkDataPreds[HkDataPreds$station == this_station,]$stationPred != this_station)
                           
  incorrectly_predicted <- (1 - (correctly_predicted + correct_island) ) 
  
  
  add.pie( x = station_long,y = station_lat, z = c(correctly_predicted,correct_island,incorrectly_predicted)
               ,edges=200,
               radius=0.012*(1-(1/length(HkDataPreds[HkDataPreds$station == this_station,]$stationPred)**0.5)),
               col=c("red","white","black"), labels = ""
  )
  
}




map.axes(cex.axis = 1)
map.scale( cex = 1)
par(fig = c(0.625,0.975,0.1,0.65), new = T) 
plot(HKshp,xlim = c(114.185,114.2), ylim = c(22.265,22.345), col = "grey", border = "darkgrey", bg ="lightskyblue1")
box( col = 'black')
for (i in 1:length(levels(HkDataPreds$station))){
  this_station<- levels(HkDataPreds$station)[i]
  
  station_lat <- mean(HkDataPreds[HkDataPreds$station == this_station,]$latitude)
  station_long <- mean(HkDataPreds[HkDataPreds$station == this_station,]$longitude)
  
  
  correctly_predicted <- mean(HkDataPreds[HkDataPreds$station == this_station,]$stationPred == this_station ) 
  correct_island <- mean(HkDataPreds[HkDataPreds$station == this_station,]$predIsland == HkDataPreds[HkDataPreds$station == this_station,]$Island
                         & HkDataPreds[HkDataPreds$station == this_station,]$stationPred != this_station)
  
  incorrectly_predicted <- (1 - (correctly_predicted + correct_island) ) 
  
  
  add.pie( x = station_long,y = station_lat, z = c(correctly_predicted,correct_island,incorrectly_predicted)
           ,edges=200,
           radius=0.012*(1-(1/length(HkDataPreds[HkDataPreds$station == this_station,]$stationPred)**0.25)),
           col=c("red","white","black"), labels = ""
  )
  
}

dev.off()












