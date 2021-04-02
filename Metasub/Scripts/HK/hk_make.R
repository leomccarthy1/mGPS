library(doParallel)
library(plyr)


#Import data sets 
setwd( rprojroot::find_rstudio_root_file())
source('mGPS.r')
metasub_data <-  read.csv(file = "Data/Metasub/msub_meta_taxa.csv", header = TRUE)


#Extract the samples from hong kong for which the station of origin is known
hong_kong <- hong_kong <- metasub_data[metasub_data$city == "hong_kong",]
hong_kong <- hong_kong[-which(hong_kong$station == ""), ]
hong_kong <-  droplevels(hong_kong[hong_kong$latitude != hong_kong$city_latitude,])
  
keep <- levels(hong_kong$station)[table(hong_kong$station) > 3]
hong_kong <- droplevels(hong_kong[hong_kong$station %in% keep, ] )


#variable selection for this dataset to remove redundant variables

  
featureEliminationHk <- species_select(x = hong_kong[,c(43:3711)],y = hong_kong$station,remove_correlated = F,subsets = c(1000,200,300,400,500),cores = 8)
optimumVarsHk <- featureEliminationHk$optVariables[1:200]



v <- varImp(featureEliminationHk$fit, type = 1, scale = F)
v[,"taxa"] <- row.names(v)
v <- v[order(v$Overall,decreasing = T),]
dir.create('Metasub/Outputs/HK', showWarnings = FALSE)
write.csv(v, file = "Metasub/Outputs/HK/hk_git.csv")



#Apply our mGPS algorithm
  
#generate 5 folds for test predictions.
set.seed(18)
trainFoldsHk <-  createFolds(hong_kong$station, k = 5, returnTrain = T)
  
GeoPredsHk <- list()
#iteratively train the model on each of the 5 training folds and generate predictions using the coresponding test fold.
HKshp <- rgdal::readOGR(dsn ="Data/Geo/gadm36_HKG_shp", layer = "gadm36_HKG_1", verbose = F)
HKcoords <- cbind(maps::SpatialPolygons2map(HKshp)$x, maps::SpatialPolygons2map(HKshp)$y)
HKcoords <- HKcoords[complete.cases(HKcoords),]

for (i in 1:5){
    
  train <- hong_kong[trainFoldsHk[[i]],]
  test <- hong_kong[-trainFoldsHk[[i]],]
    
  testPreds <-mGPS(training = train, testing = test, classTarget = "station",variables = optimumVarsHk,nthread = 8,hierarchy = c('station','latitude','longitude'), coast = HKcoords)
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

write.csv(HkDataPreds,"Metasub/Outputs/HK/hk_results.csv")


























