
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


#import bacterial and meta data
metasub_data <-  read.csv(file = "DATA/msub_meta_taxa.csv", header = TRUE)


##### get samples from new york with known long an lat ####

new_york <- metasub_data[metasub_data$city == "new_york_city",]
new_york  <-  droplevels(new_york[new_york$latitude != new_york$city_latitude,])


new_york_stations <- read.csv("/Users/leomccarthy/Documents/Microbiome-mapping/Geodata/NYstations.csv",header = TRUE)
  new_york_stations <- new_york_stations[!duplicated(new_york_stations$Station),]
  new_york_stations$Station <- make.names(new_york_stations$Station)


find_station<- function(long,lat){
    distances_from<-  spDistsN1(cbind(new_york_stations$Longitude, new_york_stations$Latitude) , c(long,lat), longlat = TRUE)
    
    closest_station <-  new_york_stations[which.min(distances_from),"Station"]
  
    
    return(closest_station)
    
  }



  
stations_found <-  mapply(find_station, long = new_york$longitude, lat = new_york$latitude )
new_york$station <- factor(stations_found)
new_york_stations <- droplevels(new_york_stations[ new_york_stations$Station %in% c(levels(new_york$station)) ,])


### remove stations with less than 3 samples ####

keep <- levels(new_york$station)[table(new_york$station) > 2]
new_york <-droplevels(new_york[new_york$station %in% keep, ])



#Remove redundant variables 


#Activate parallel processing 
registerDoParallel(parallel::detectCores() - 1)

#recursive feature elimination for the target variable "station"
rfe_ctrl <- rfeControl(functions = rfFuncs,
                       method = "cv",
                       number =  5,
                       verbose = FALSE,
                       allowParallel = TRUE
)
set.seed(123)
featureEliminationNY <- rfe(x = new_york[, c(43:3711)],y = new_york$station,
                          sizes = c(100,200,300,400,500),
                          rfeControl = rfe_ctrl,
                          tuneLength = 10
)

optimumVarsNY<- featureEliminationNY$optVariables


#Extract and rank species importance

v <- varImp(featureEliminationNY$fit, type = 1, scale = F)
v[,"taxa"] <- row.names(v)
v <- v[order(v$Overall,decreasing = T),]
top_species <- v[1:25,"taxa"]

#plot top 25 geo indicator species fig s13
par(font = 3)
dotchart(rev(v[1:25,"Overall"])*100,labels= rev(top_species),cex=1.2,pt.cex = 1.3,
         xlab="Mean decrease in accuracy", mgp = c(2.2,0,0))



## local adaptation of mGPS ommitting continent level

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
#generate 5 stratified folds for test predictions.
set.seed(18)
trainFolds <-  createFolds(new_york$station, k = 5, returnTrain = T)


GeoPredsNY <- list()
registerDoParallel(detectCores() - 1) 
#iteratively train the model on each of the 10 training folds and generate predictions using the coresponding test fold.
for (i in 1:5){
  
  train <- new_york[trainFolds[[i]],]
  test <- new_york[-trainFolds[[i]],]
  
  testPreds <-mGPS_local(training = train, testing = test, classTarget = "station",variables = optimumVarsNY)
  GeoPredsNY[[i]] <- testPreds
  
}



add_preds <- list()
for (i in 1:5){
  
  add_preds[[i]] <- cbind(new_york[-trainFolds[[i]],] , 
                          "stationPred"= GeoPredsNY[[i]][[1]], 
                          "latPred" = GeoPredsNY[[i]][[2]], 
                          "longPred" = GeoPredsNY[[i]][[3]] )
  
  
  
}

NYDataPreds <- rbind.fill(add_preds)


####plot ####

NY_stations_set <- new_york_stations[ new_york_stations$Station %in% c(levels(new_york$station)) ,]
NY_stations_set <- droplevels(NY_stations_set[order(NY_stations_set$Latitude),])
lats_of_stations <- NY_stations_set$Latitude


#### borough ####

NYBorough <- readOGR(dsn = "Geodata/Borough Boundaries (Water Areas Included)", layer = "geo_export_830a0aea-d1d1-43f6-8349-8c576c70d38e")


find_borough <- function(long,lat){
  
  
  
  if(point.in.polygon(long,lat, SpatialPolygons2map(NYBorough[1,])$x,SpatialPolygons2map(NYBorough[1,])$y) %in% c(1,2,3)){
    borough <- NYBorough[1,]$boro_name
  }
  
  else if(point.in.polygon(long,lat, SpatialPolygons2map(NYBorough[2,])$x,SpatialPolygons2map(NYBorough[2,])$y)%in% c(1,2,3)){
    borough <- NYBorough[2,]$boro_name
  
  }
  else if(point.in.polygon(long,lat, SpatialPolygons2map(NYBorough[3,])$x,SpatialPolygons2map(NYBorough[3,])$y) %in% c(1,2,3)){
    borough <- NYBorough[3,]$boro_name
  }
  else if(point.in.polygon(long,lat, SpatialPolygons2map(NYBorough[4,])$x,SpatialPolygons2map(NYBorough[4,])$y)%in% c(1,2,3)){
    borough <- NYBorough[4,]$boro_name
  }
  
  else if(point.in.polygon(long,lat, SpatialPolygons2map(NYBorough[5,])$x,SpatialPolygons2map(NYBorough[5,])$y)%in% c(1,2,3)){
    borough <- NYBorough[5,]$boro_name
  }
  
  return(borough)
}

NYDataPreds[,"Borough"] <- droplevels(mapply(find_borough, NYDataPreds$longitude,NYDataPreds$latitude))
NYDataPreds[,"PredBorough"] <- droplevels(mapply(find_borough, NYDataPreds$longPred,NYDataPreds$latPred))
mean(NYDataPreds$PredBorough == NYDataPreds$Borough)



#Ds]istance from origin
for (i in 1:nrow(NYDataPreds)){
  NYDataPreds[i,"Distance_from_origin"] <- 
    distm(c(NYDataPreds[i,"longPred"],NYDataPreds[i,"latPred"]), c(NYDataPreds[i,"longitude"],NYDataPreds[i,"latitude"]), fun = distHaversine)/1000
}
median(NYDataPreds[,"Distance_from_origin"] )
confusionMatrix(NYDataPreds$station, NYDataPreds$stationPred)


#Plotting Fig 3, will be saved to working directory 
library(raster)
pie_lats <- c(40.6819,40.737,40.7319)
pie_longs <- c(-73.967,-74.055,-73.84)
palette <-c( "gold2","darkorchid4","deeppink2")

png("Fig3.png", width = 12, height = 8, units = 'in', res = 600)
plot(NYshp ,xlim = c(-74,-73.85),ylim = c(40.67,40.8),col = "grey", xlab ="",ylab ="", bg= "lightskyblue1", lwd = 1.2,border = "grey40")
title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
lines(NYlines, col = "brown3", cex = 0.8)
for (i in 1:3){
  this_borough <- levels(NYDataPreds$Borough)[i]
  stations <- levels(droplevels(NYDataPreds[NYDataPreds$Borough == this_borough,"station"]))
  
  find_lats <- NYDataPreds[NYDataPreds[,"Borough"] == this_borough,]$latPred
  find_longs <- NYDataPreds[NYDataPreds[,"Borough"] == this_borough,]$longPred
  
  station_lat <- NY_stations_set[NY_stations_set[,"Station"] %in% stations,]$Latitude
  station_long <- NY_stations_set[NY_stations_set[,"Station"] %in% stations,]$Longitude
  
  points(find_longs,find_lats, col = palette[i], pch  = "+", cex = 1.3)
  points(station_long,station_lat, col = palette[i], pch  = 17, cex = 1.6)
  
  correct <- mean(NYDataPreds[NYDataPreds$Borough == this_borough,"PredBorough"] == this_borough)
  add.pie(z = c(correct,1-correct) , x = pie_longs[i], 
          y = pie_lats[i], radius = 0.01, col = c(palette[i],"black"), labels = "")

  
}
for (i in 1:3){
  this_borough <- levels(NYDataPreds$Borough)[i]
  label_lats <- c(-73.94,-74.028,-73.815)
  text(label_lats[i],pie_lats[i] ,labels = this_borough, col = palette[i], font = 2)
}
legend(-74.08,40.85,legend = c("Stations","Predicted sample locations"),pch = c(17,3),col ="darkgreen", cex = 0.7)
map.axes(cex.axis = 0.8)
map.scale( cex = 0.5)
par(fig = c(0,0.3,0.5,1), new = T) 
plot(NYshp,col = "grey", bg ="lightskyblue1",border = "grey40")
points(NYDataPreds$city_longitude, NYDataPreds$city_latitude, cex = 1.5, pch = 16)
text(-73.95,40.68, labels = "New York", pch = 1.35, font = 2)
rect(-74.1, 40.6 ,-73.8, 40.8, density = NULL, angle = 45,
     col = NA, border = NULL, pch = 2)
box( col = 'black')

dev.off()


# abund station plot 
levels(NYDataPreds$station)[levels(NYDataPreds$station)=="Brooklyn.Bridge.City.Hall.Chambers.St"] <- "City.Hall.Chambers.St"
levels(NYDataPreds$station)[levels(NYDataPreds$station)=="Sutphin.Blvd.Archer.Av...JFK" ] <- "Sutphin.Blvd.Archer.Av" 
ag <-  aggregate(NYDataPreds[,top_species], by = list(NYDataPreds$station), FUN = median)




for (i in top_species){
  ag[,i] <- (ag[,i] - min(ag[,i]))/(max(ag[,i]) - min(ag[,i]))
}
data.long <- melt(ag)

library(tidyr)
library(hrbrthemes)
library(reshape2)
library(viridis)


ggplot(data = data.long, mapping = aes(x = Group.1,
                                       y = forcats::fct_rev(variable),
                                       fill = value)) +
  geom_tile()+
  scale_fill_viridis(begin = 0, end = 1, limits = c(0,1)) +
  xlab("Station")+
  ylab("")+
  labs(fill="Relative abundance \n(normalised)")+
 theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1, size = 11),
       axis.text.y=element_text(size = 12, face = "italic"),
       legend.text = element_text(size = 11), 
       legend.title = element_text(size = 12))


ggsave(
  "NY_FigS13.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 13,
  height = 8,
  
  dpi = 600,
  limitsize = TRUE)








