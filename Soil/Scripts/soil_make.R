### Data load ###
source("mGPS.R")
library(caret)
setwd( rprojroot::find_rstudio_root_file())
environmental <- read.csv(file="Data/Soil/Dataset_01_22_2018_enviro.csv",header=TRUE)
taxa <- read.csv(file="Data/Soil/Dataset_01_22_2018_taxa.csv",header=TRUE)
SoilData <- merge(environmental,taxa,by.x="ID_Environmental",by.y="ID_Environmental")


### Data clean ###
SoilData$country <- maps::map.where(database="world", SoilData$Longitude, SoilData$Latitude)
SoilData[234:237,]$country <- "Chile"
SoilData[c(2,6,9,15,22,67,34),]$country <- "USA"

SoilData$country <- factor(make.names(SoilData$country))

#rename columns
colnames(SoilData)[3:5] <- c("longitude","latitude","continent")

#Remove countries with insufficient samples
keep <- levels(SoilData$country)[table(SoilData$country) > 2]
SoilData <- droplevels(SoilData[SoilData$country%in% keep, ])

#concert read data to relative abundance.
for (i in 28:538){
  SoilData[,i] <- (SoilData[,i]/10000)
}


### Get GIT's ###
featureElim <- species_select(x = SoilData[,c(28:538)],y = SoilData$country,remove_correlated = F,subsets = c(20,30,50,100,200),cores = 8)
optVars <- featureElim$optVariables

v <- varImp(featureElim$fit, type = 1, scale = F)
v[,"taxa"] <- row.names(v)
v <- v[order(v$Overall,decreasing = T),]
dir.create('Soil/Outputs', showWarnings = FALSE)
write.csv(v, file = "Soil/Outputs/soil_git.csv")


### Make predictions ###
coastlines <- cbind("x"  = maps::SpatialLines2map(rworldmap::coastsCoarse)$x ,"y" =maps::SpatialLines2map(rworldmap::coastsCoarse)$y)
coastlines <- coastlines[complete.cases(coastlines),]
coastlines <- coastlines[coastlines[,1] < 180 ,]
#5 fold cross validation
set.seed(18)
trainFolds <-  createFolds(SoilData$country, k = 5, returnTrain = T)

GeoPredsSoil <- list()
registerDoParallel(7) 

for (i in 1:5){
  
  train <- SoilData[trainFolds[[i]],]
  test <- SoilData[-trainFolds[[i]],]
  
  testPreds <- mGPS(training = train, testing = test, classTarget = "country",variables = optVars,nthread = 8,hierarchy = c('continent','country','latitude','longitude'), coast=coastlines)
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
write.csv(SoilDataPreds,"Soil/Outputs/soil_results.csv")
