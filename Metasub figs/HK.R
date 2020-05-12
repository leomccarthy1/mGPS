  

library(sp)
library(rworldmap)
library(rpart)
library(maps)
library(rgdal)
library(randomForest)
library(geosphere)
library(caret)
library(xgboost)
library(e1071)
library(randomcoloR)
library(doParallel)



#Import data sets 
metasub_data <-  read.csv(file = "DATA/msub_meta_taxa.csv", header = TRUE)


#Extract the samples from hong kong for which the station of origin is known
  hong_kong <- hong_kong <- metasub_data[metasub_data$city == "hong_kong",]
  hong_kong <- hong_kong[-which(hong_kong$station == ""), ]
  hong_kong <-  droplevels(hong_kong[hong_kong$latitude != hong_kong$city_latitude,])
  
  keep <- levels(hong_kong$station)[table(hong_kong$station) > 3]
  hong_kong <- droplevels(hong_kong[hong_kong$station %in% keep, ] )


#variable selection for this dataset to remove redundant variables

  
  registerDoParallel(parallel::detectCores() - 1)
  
  
  rfe_ctrl <- rfeControl(functions = rfFuncs,
                         method = "cv",
                         number =  5,
                         verbose = FALSE,
                         allowParallel = TRUE
  )
 
  set.seed(123)
  featureEliminationHk <- rfe(x = hong_kong[, c(43:3711)],y = hong_kong$station,
                            sizes = c(100,200,300,400,500),
                            rfeControl = rfe_ctrl,
                            tuneLength = 3
  )
  


  optimumVarsHk <- featureEliminationHk$optVariables[1:200]

  
  #Extract and rank species importance
  v <- varImp(featureEliminationHk$fit, type = 1, scale = F)
  v[,"taxa"] <- row.names(v)
  v <- v[order(v$Overall,decreasing = T),]
  top_species <- v[1:25,"taxa"]
  

  #plot
  
  par(font = 3)
  dotchart(rev(v[1:25,"Overall"])*100,labels= rev(top_species),cex=1.2,pt.cex = 1.3,
           xlab="Mean decrease in accuracy", mgp = c(2.2,0,0))
  
  
#local adaptation of mGPS model to ommit continent level.
  
mGPS_local <-  function(training,testing,classTarget,variables){
    
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
    
    
    #Xgb_region <- train(x = training[,variables],y = training[,"continent"],
                        #method = "xgbTree",
                        #trControl = trControlClass,
                        #tuneGrid = tune_grid,
                        #nthread = 1)
    
    #l1_train <- data.frame(training[,c(variables)],Xgb_region[["pred"]][order(Xgb_region$pred$rowIndex),levels(training[,"continent"]) ])
    
    l1_train <- training[,c(variables)]
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
    
    
    
    
    
    #regProbs <- predict(Xgb_region, newdata = testing[,variables],type ="prob")
    
    #l1_test <- data.frame(testing[,variables], regProbs)
    
    l1_test <- testing[,variables]
    
    classPred <- predict(Xgb_class, newdata = l1_test)
    classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
    
    l2_test <-  data.frame(l1_test, classProbs) 
    latPred <- predict(Xgb_latitude, newdata = l2_test)
    
    l3_test <- data.frame(l2_test, latPred)
    longPred <- predict(Xgb_longitude, newdata = l3_test)
    
    
    
    
    
    
    return(list(classPred, latPred, longPred))
    
}
  

#Apply our mGPS algorithm
  
  #generate 5 folds for test predictions.
  set.seed(18)
  trainFoldsHk <-  createFolds(hong_kong$station, k = 5, returnTrain = T)
  
  
  GeoPredsHk <- list()
  registerDoParallel(parallel::detectCores() - 1) 
  #iteratively train the model on each of the 5 training folds and generate predictions using the coresponding test fold.
  for (i in 1:5){
    
    train <- hong_kong[trainFoldsHk[[i]],]
    test <- hong_kong[-trainFoldsHk[[i]],]
    
    testPreds <-mGPS_local(training = train, testing = test, classTarget = "station",variables = optimumVarsHk)
    GeoPredsHk[[i]] <- testPreds
    
  }
 

  add_preds <- list()
  for (i in 1:5){
    
    add_preds[[i]] <- cbind(hong_kong[-trainFoldsHk[[i]],] , 
                            "stationPred"= GeoPredsHk[[i]][[1]], 
                            "latPred" = GeoPredsHk[[i]][[2]], 
                            "longPred" = GeoPredsHk[[i]][[3]] )
    
    
    
  }
  
library(plyr)
HkDataPreds <- rbind.fill(add_preds)

station_names <- c("Admiralty"     ,"Che Kung Temple"  ,  "Cheung Chau",       
"Cheung Sha Wan"   ,  "East Tsim Sha Tsui", "Fo Tan",            
"Fortress Hill"   ,   "Hung Hom"       ,    "Kowloon Tong" ,     
"Lai King",           "Mongkok East",       "North Point"     ,  
"Peng Chau",          "Prince Edward",      "Quarry Bay" ,       
"Residence 1",        "Residence 2",        "Residence 3",       
"Residence 4",        "Sai Kung",           "Sai Wan Ho Pier",   
"Sham Shui Po",       "Shek Kip Mei",       "Tai Mei Tuk",       
"Tai O",              "Tai Wai",            "Tsim Sha Tsui",     
"Tsing Yi",           "Tuen Mun Pier",      "Wan Chai",          
"Wu Kai Sha Pier",    "Yau Tong",           "Yung Shue Wan") 




##Final stage is to adjust any co-ordinates that are in the sea and pull them to the nearset land mass
#get HK coastlines 
HKshp <- readOGR(dsn ="Geodata/gadm36_HKG_shp", layer = "gadm36_HKG_1")
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




#Get distance from origin 
for (i in 1:nrow(HkDataPreds)){
  HkDataPreds[i,"Distance_from_origin"] <- distm(c(HkDataPreds[i,"longPred"],HkDataPreds[i,"latPred"]), c(HkDataPreds[i,"longitude"],HkDataPreds[i,"latitude"]), fun = distHaversine)/1000
}
median(HkDataPreds$Distance_from_origin)
postResample(HkDataPreds$latitude, HkDataPreds$latPred)

##Plots Fig S4 will be saved to working directory

library(randomcoloR)
#color palette 
n <- length(levels(HkDataPreds$station))
set.seed(15)
palette <- randomColor(count = n, luminosity = "bright")

#plot predicted GPS locations
map <- getMap(resolution = "high")
palette <- randomColor(count = n, luminosity = "bright")
png("FigS4.png", width = 12, height = 8, units = 'in', res = 600)
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
legend(114.3,22.45,legend=station_names, col = palette,pch = 17, cex = 0.8, bg = "lightskyblue1")
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




## Fig 2 - , will be saved to working directory

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
mean(HkDataPreds$Island == HkDataPreds$predIsland)
library(mapplots)

png("Fig2.png", width = 12, height = 8, units = 'in', res = 600)
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






# Fig S13, will be saved to working directory
levels(hong_kong$station) <- station_names
ag <-  aggregate(hong_kong[,top_species], by = list(hong_kong$station), FUN = median)
for (i in top_species){
  
  ag[,i] <- ag[,i] <- (ag[,i] - min(ag[,i]))/(max(ag[,i]) - min(ag[,i]))
}
data.long <- melt(ag)

library(tidyr)
library(hrbrthemes)
library(reshape2)
library(viridis)

ggplot(data = data.long, mapping = aes(x = Group.1,
                                       y =forcats::fct_rev(variable),
                                       fill = value)) +
  geom_tile()+
  scale_fill_viridis(limits = c(0,1)) +
  xlab("Station")+
  ylab("")+
  labs(fill="Relative abundance \n(normalised)")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1, size = 11),
        axis.text.y=element_text(size = 12, face = "italic"),
        
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 12))

ggsave(
  "Hk_FigS13.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 13,
  height = 8,
  
  dpi = 600,
  limitsize = TRUE)









