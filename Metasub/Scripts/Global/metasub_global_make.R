source('mGPS.r')

#Import data sets 
setwd( rprojroot::find_rstudio_root_file())
complete_meta <- read.csv(file = "Data/Metasub/complete_metadata.csv", header = TRUE)
taxa_abund <-read.csv("Data/Metasub/metasub_taxa_abundance.csv", header = T)
taxa_abund <- unique(taxa_abund)

#merge bacterial and meta data
metasub_data <- merge(complete_meta,taxa_abund,by.x="uuid",by.y="uuid")

#Remove control samples
control_samples <- c( which(metasub_data$city %in% c("control", "other_control","neg_control","other","pos_control")), which(metasub_data$control_type %in% c("ctrl cities","negative_control","positive_control"))) 
metasub_data <- droplevels(metasub_data[-c(control_samples), ])


#re-label london boroughs 
metasub_data$city[metasub_data$city %in% c("kensington","islington")] <- "london" 
metasub_data <- droplevels(metasub_data)

#remove sparse samples locations and dubiously labelled samples. 

small_cities <-  names(which(summary(metasub_data$city) < 8))
remove_samples <- which(metasub_data$city %in%  c("antarctica", small_cities))
metasub_data <- droplevels(metasub_data[-c(remove_samples), ])


#Correction of identified misslabelling of data 
metasub_data$latitude[metasub_data$city == "kyiv"] <- metasub_data$city_latitude[metasub_data$city == "kyiv"]
metasub_data$longitude[metasub_data$city == "kyiv"] <- metasub_data$city_longitude[metasub_data$city == "kyiv"]
metasub_data$continent[metasub_data$city == "porto"] <- "europe"

metasub_data[is.na(metasub_data$latitude),]$latitude <- metasub_data[is.na(metasub_data$latitude),]$city_latitude
metasub_data[is.na(metasub_data$longitude),]$longitude <- metasub_data[is.na(metasub_data$longitude),]$city_longitude

#correction to some incorrect city co-ords for a few london samples
metasub_data[metasub_data$city == "london",]$city_latitude <- 51.50853
metasub_data[metasub_data$city == "london",]$city_longitude <- -0.12574


### Find GITs ####
featureElim <- species_select(x = metasub_data[, c(43:3711)],y = metasub_data$city,remove_correlated = F, c(50,100,200,300,500,1500),cores = 8)
optVars <- featureElimination$optVariables


v <- varImp(featureElim$fit, type = 1, scale = F)
v[,"taxa"] <- row.names(v)
v <- v[order(v$Overall,decreasing = T),]
dir.create('Metasub/Outputs', showWarnings = FALSE)
write.csv(v, file = "Metasub/Outputs/metasub_global_git.csv")
write.csv(data.frame("n_vars" = featureElimination$results$Variables, "accuracy" = featureElimination$results$Accuracy), file = "Metasub/Outputs/metasub_git_subsets.csv")

### Make predictions using mGPS ###
coastlines <- cbind("x"  = maps::SpatialLines2map(rworldmap::coastsCoarse)$x ,"y" =maps::SpatialLines2map(rworldmap::coastsCoarse)$y)
coastlines <- coastlines[complete.cases(coastlines),]
coastlines <- coastlines[coastlines[,1] < 180 ,]
#generate 5 stratified folds for test predictions.
set.seed(18)
trainFolds <-  createFolds(metasub_data$city, k = 5, returnTrain = T)



GeoPreds <- list()

#iteratively train the model on each of the 5 training folds and generate predictions using the coresponding test fold.
for (i in 1:5){
  
  train <- metasub_data[trainFolds[[i]],]
  test <- metasub_data[-trainFolds[[i]],]
  
  testPreds <-mGPS(training = train, testing = test, classTarget = "city",variables = optVars,nthread = 8,hierarchy = c('continent','city','latitude','longitude'), coast=coastlines)
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
write.csv(MetasubDataPreds,"Metasub/Outputs/metasub_global_results.csv")




