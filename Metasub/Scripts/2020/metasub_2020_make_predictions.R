
library(caret)
source("mGPS.r")

#Load 2020 replication data sets (geo/meta data and taxa read data)
rep_taxa <- read.csv(file="Data/Metasub/metasub_covid-19.cap2__capalyzer-v0_1_0__kraken2_taxa.1603842648-1762.read_counts (1).csv",header=TRUE)
rep_meta <- read.csv(file="Data/Metasub/MetaSUB_COVID-19_metadata.csv",header=TRUE, stringsAsFactors = FALSE)
metasub <- read.csv(file = "Data/Metasub/metasub_data_full.csv", header = TRUE)[,-1]
msub_taxa <-  names(metasub[,43:ncol(metasub)])

#Clean ataxa data
for (i in 2:ncol(rep_taxa)){
  names(rep_taxa)[i] <- strsplit(colnames(rep_taxa)[i],"s__")[[1]][2]
  
}
rep_taxa <- rep_taxa[!is.na(names(rep_taxa))]
rep_taxa[is.na(rep_taxa)] <- 0
taxa <- names(rep_taxa[,2:ncol(rep_taxa)])
rep_taxa[,taxa] <- data_normalise(rep_taxa[,taxa])
taxa <- taxa[taxa %in% msub_taxa]

#remove samples with very incorrect looking lat longs
rep_meta <- rep_meta[rep_meta$latitude != 0.1 & rep_meta$longitude != 0.1,]

rep_meta$city <- as.factor(make.names(rep_meta$city))

na <- c("United States","Canada")
sa <- c("Brazil","Argentina","Chile")
ea <- c("Hong Kong","Republic of Korea","japan","Japan","China")
me <- c("Pakistan","Israel")
eu <- c("Germany","Poland","United Kingdom","Sweden","Switzerland")

rep_meta[rep_meta$country %in% na,'continent'] <- "North America"
rep_meta[rep_meta$country %in% sa,'continent'] <- "South America"
rep_meta[rep_meta$country %in% ea,'continent'] <- "East Asia"
rep_meta[rep_meta$country %in% me,'continent'] <- "Middle East"
rep_meta[rep_meta$country %in% eu,'continent'] <- "Europe"
rep_meta$continent <- as.factor(make.names(rep_meta$continent))

rep_meta <- droplevels(rep_meta[rep_meta$city != "X",])
rep_meta[rep_meta$city == "Seattle..WA","city"] <- as.factor("Seattle")
for(i in levels(rep_meta$city)){
  
  rep_meta[(is.na(rep_meta$latitude)) & (rep_meta$city == i),"latitude"] <- mean(rep_meta[rep_meta$city == i,"latitude"], na.rm = T)
  rep_meta[(is.na(rep_meta$longitude)) & (rep_meta$city == i),"longitude"] <- mean(rep_meta[rep_meta$city == i,"longitude"], na.rm = T)
}
rep_meta <- rep_meta[!is.na(rep_meta$latitude),]
#krakow lats and longs are wrong way round in origionl data 
rep_meta[rep_meta$city == "Krakow",c("latitude","longitude")] <- rep_meta[rep_meta$city == "Krakow",c("longitude","latitude")]

#Merge geodata and taxa data for each sample
rep_metasub <-droplevels( merge(rep_meta[,c("ID","city","latitude","longitude","continent")],rep_taxa,by.x = "ID", by.y = "X"))


##### Get optimum taxa (GITs) #####
taxa <- taxa[!(taxa %in% names(which(colSums(rep_metasub[,taxa]) == 0)))]
featureElim <- species_select(x = rep_metasub[,taxa], rep_metasub$city,remove_correlated = T, subsets = c(500,250,200,100),cores = 8)
optVars <- featureElim$optVariables[1:250]
v <- varImp(featureElim$fit, type = 1, scale = F)
v[,"taxa"] <- row.names(v)
v <- v[order(v$Overall,decreasing = T),]
dir.create('Metasub/Outputs/2020', showWarnings = FALSE)
write.csv(v, file = "Metasub/Outputs/2020/metasub_2020_git.csv")


####Get predictions###
coastlines <- cbind("x"  = maps::SpatialLines2map(rworldmap::coastsCoarse)$x ,"y" =maps::SpatialLines2map(rworldmap::coastsCoarse)$y)
coastlines <- coastlines[complete.cases(coastlines),]
coastlines <- coastlines[coastlines[,1] < 180 ,]


set.seed(123)
trainFolds <-  createFolds(rep_metasub$city, k = 5, returnTrain = T)

GeoPreds <- list()
doParallel::registerDoParallel(1) 
#iteratively train the model on each of the 5 training folds and generate predictions using the coresponding test fold.
for (i in 1:5){
  
  train <- rep_metasub[trainFolds[[i]],]
  test <- rep_metasub[-trainFolds[[i]],]
  
  testPreds <-mGPS(training = train, testing = test, classTarget = "city",variables = optVars,nthread = 8,hierarchy = c('continent','city','longitude','latitude'), coast=coastlines)
  GeoPreds[[i]] <- testPreds
  
}


add_preds <- list()
for (i in 1:5){
  
  add_preds[[i]] <- cbind(rep_metasub[-trainFolds[[i]],] , 
                          "cityPred"= GeoPreds[[i]][[1]], 
                          "latPred" = GeoPreds[[i]][[2]], 
                          "longPred" = GeoPreds[[i]][[3]] )
  
  
  
}

replication_preds <- plyr::rbind.fill(add_preds)




####seave results to file#####
write.csv(replication_preds[,c('ID','continent','city','longitude','latitude','cityPred','longPred','latPred')],"Metasub/Outputs/2020/metasub_2020_results.csv")





