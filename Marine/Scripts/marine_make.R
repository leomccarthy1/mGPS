library(caret)
source("mGPS.r")
#Load Data
setwd( rprojroot::find_rstudio_root_file())
marine_taxa =  read.csv(file="Data/Marine/marine_taxa.csv",header=TRUE)
taxa <- names(marine_taxa[,c(2:2512)])
marine_taxa$Sea = factor(make.names(marine_taxa$Sea))



#Normalise data
marine_taxa[is.na(marine_taxa)] <- 0
marine_taxa[,taxa] <- marine_taxa[,taxa][,-which(colSums(marine_taxa[,taxa]) == 0)] 
marine_taxa[,taxa] <- data_normalise(marine_taxa[,taxa])


#Get GITs
featureElimination <- species_select(x = marine_taxa[,taxa],y = marine_taxa$Sea, remove_correlated = T, subsets = c(5000,500,300,200,100),cores = 7)
opt <- featureElimination$optVariables


# Get predictions using cross validatino
set.seed(123)
trainFolds <-  createFolds(marine_taxa$Sea, k = 5, returnTrain = T)

GeoPreds <- list()
doParallel::registerDoParallel(1) 
#iteratively train the model on each of the 5 training folds and generate predictions using the coresponding test fold.
for (i in 1:5){
  
  train <- marine_taxa[trainFolds[[i]],]
  test <- marine_taxa[-trainFolds[[i]],]
  
  testPreds <-mGPS(training = train, testing = test, classTarget = "Sea",variables = opt,nthread = 8,hierarchy = c('Sea','latitude','longitude'))
  GeoPreds[[i]] <- testPreds
  
}

add_preds <- list()
for (i in 1:5){
  
  add_preds[[i]] <- cbind(marine_taxa[-trainFolds[[i]],] , 
                          "seaPred"= GeoPreds[[i]][[1]], 
                          "latPred" = GeoPreds[[i]][[2]], 
                          "longPred" = GeoPreds[[i]][[3]] )
  
  
  
}

marine_preds <- plyr::rbind.fill(add_preds)


#Pull to nearest marine coastline
mins = c()
for (i in 1:nrow(marine_taxa)){
  dists <- sp::spDistsN1( as.matrix(marine_taxa[-i,c("latitude","longitude")]), c(marine_taxa[i,"latitude"],marine_taxa[i,"longitude"]), longlat = TRUE)
  mins <-  c(mins,min(dists))
  
}

seas <- rgdal::readOGR(dsn = "Data/Geo/ne_10m_geography_marine_polys", layer = "ne_10m_geography_marine_polys")
coastlines <- cbind("x"  =maps::SpatialPolygons2map(seas)$x ,"y" =maps::SpatialPolygons2map(seas)$y)
coastlines <- coastlines[complete.cases(coastlines),]
coastlines <- coastlines[coastlines[,1] < 180 ,]


find_coast <- function(long,lat){
  distances_from_coastline <-  sp::spDistsN1(coastlines , c(long,lat), longlat = TRUE)
  
  closest_point <-  which.min(distances_from_coastline)
  new_coords <- coastlines[closest_point,]
  
  return(new_coords)
  
}
data(maptools::wrld_simpl)

## Create a SpatialPoints object
set.seed(0)
points <- data.frame(marine_preds$longPred, marine_preds$latPred) 
pts <- sp::SpatialPoints(points, proj4string=CRS(proj4string(wrld_simpl)))
## Find which points fall over land
ii <- !is.na(over(pts, wrld_simpl)$FIPS)
toAdjust <- marine_preds[which(ii == TRUE),]
adjusted <- mapply(find_coast, long = toAdjust$longPred, lat = toAdjust$latPred )

marine_preds[which(ii == TRUE), "latPred"] <- adjusted[2,]
marine_preds[which(ii == TRUE), "longPred"] <- adjusted[1,]

write.csv(marine_preds[,c('Sample','Sea','longitude','latitude','seaPred','longPred','latPred')],"Marine/Outputs/marine_results.csv")
