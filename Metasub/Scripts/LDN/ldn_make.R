

library(caret)
library(xgboost)
library(e1071)
library(doParallel)

#Import data sets 
setwd(rprojroot::find_rstudio_root_file())
source('mGPS.r')
metasub_data <- metasub_data <-  read.csv(file = "Data/Metasub/msub_meta_taxa.csv", header = TRUE)
london_samples <-metasub_data[metasub_data$city == "london",]

#london station data
london_stations <- read.csv(file = "Data/Geo/London stations GPS.csv", header = TRUE)

#Remove samples with unknown exact Lat+long
london_samples <-  droplevels(london_samples[london_samples$latitude != london_samples$city_latitude,])


find_stationLDN<- function(long,lat){
  distances_from<-  spDistsN1(cbind(london_stations$Longitude, london_stations$Latitude) , c(long,lat), longlat = TRUE)
  closest_station <-  london_stations[which.min(distances_from),"Station"]
  return(closest_station)
  
}

stations_found <-  mapply(find_stationLDN, long =london_samples$longitude, lat = london_samples$latitude )
london_samples$station <- factor(make.names(stations_found))

keep<- levels(london_samples$station)[table(london_samples$station) > 3]
london_samples <- droplevels(london_samples[london_samples$station %in% keep, ])


#cluster stations into 15 spatial clsuters, "super-stations"
set.seed(123)
clust <- kmeans(london_samples[,c("longitude","latitude")], centers = 10)
london_samples[,"super_station"] <- factor(make.names(clust$cluster))




#recursive feature elimination for the target variable "station"
featureEliminationLDN <- species_select(x = london_samples[, c(43:3711)],y = london_samples$super_station,remove_correlated = F,subsets = c(100,200,300,400,500),cores = 8)
optimumVarsLDN<- featureEliminationLDN$optVariables[1:200]



#Extract and rank species importance
v <- varImp(featureEliminationLDN$fit, type = 1, scale = F)
v[,"taxa"] <- row.names(v)
v <- v[order(v$Overall,decreasing = T),]
dir.create('Metasub/Outputs/LDN', showWarnings = FALSE)
write.csv(v, file = "Metasub/Outputs/LDN/ldn_git.csv")

## Apply local mGPS model

#generate 5 folds for test predictions.
trainFolds <-  createFolds(london_samples$super_station, k = 5, returnTrain = T)


GeoPredsLDN <- list()
registerDoParallel(7) 
#iteratively train the model on each of the 10 training folds and generate predictions using the coresponding test fold.
for (i in 1:5){
  
  train <- london_samples[trainFolds[[i]],]
  test <- london_samples[-trainFolds[[i]],]
  
  testPreds <-mGPS(training = train, testing = test, classTarget = "super_station",variables = optimumVarsLDN,nthread = 8,hierarchy = c('super_station','latitude','longitude'))
  GeoPredsLDN[[i]] <- testPreds
  
}


add_preds <- list()
for (i in 1:5){
  
  add_preds[[i]] <- cbind(london_samples[-trainFolds[[i]],] , 
                          "super_stationPred"= GeoPredsLDN[[i]][[1]], 
                          "latPred" = GeoPredsLDN[[i]][[2]], 
                          "longPred" = GeoPredsLDN[[i]][[3]] )
  
  
  
}

LDNDataPreds <- rbind.fill(add_preds)

write.csv(LDNDataPreds,"Metasub/Outputs/LDN/ldn_results.csv")
























