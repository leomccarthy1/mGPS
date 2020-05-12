


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
metasub_data <- metasub_data <-  read.csv(file = "DATA/msub_meta_taxa.csv", header = TRUE)
london_samples <-metasub_data[metasub_data$city == "london",]

#london station data
london_stations <- read.csv(file = "Geodata/London stations GPS.csv", header = TRUE)

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




#Activate parallel processing 
registerDoParallel(parallel::detectCores() - 1)

#recursive feature elimination for the target variable "station"
rfe_ctrl <- rfeControl(functions = rfFuncs,
                       method = "cv",
                       number =  5,
                       verbose = FALSE,
                       allowParallel = TRUE,
                       seeds = seeds
)

seeds <- list(c(1,1,1,1,1,1),c(2,2,2,2,2,2),c(3,3,3,3,3,3),c(4,4,4,4,4,4),c(5,5,5,5,5,5) , 6)
featureEliminationLDN <- rfe(x = london_samples[, c(43:3711)],y = london_samples$super_station,
                            sizes = c(100,200,300,400,500),
                            rfeControl = rfe_ctrl,
                            tuneLength = 2
)

optimumVarsLDN<- featureEliminationLDN$optVariables[1:200]



#Extract and rank species importance
v <- varImp(featureEliminationLDN$fit, type = 1, scale = F)
v[,"taxa"] <- row.names(v)
v <- v[order(v$Overall,decreasing = T),]
top_species <- v[1:25,"taxa"]

#plot for fig s13 of top 25 geo indicator species 
png("ImpSpeciesLDN.png", width = 12,height = 8, units = 'in', res = 600)
par(font = 3)
dotchart(rev(v[1:25,"Overall"])*100,labels= rev(top_species),cex=1.2,pt.cex = 1.3,
         xlab="Mean decrease in accuracy", mgp = c(2.2,0,0))
dev.off()

## Apply local mGPS model
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

#generate 5 folds for test predictions.
trainFolds <-  createFolds(london_samples$super_station, k = 5, returnTrain = T)


GeoPredsLDN <- list()
registerDoParallel(7) 
#iteratively train the model on each of the 10 training folds and generate predictions using the coresponding test fold.
for (i in 1:5){
  
  train <- london_samples[trainFolds[[i]],]
  test <- london_samples[-trainFolds[[i]],]
  
  testPreds <-mGPS_local(training = train, testing = test, classTarget = "super_station",variables = optimumVarsLDN)
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




#Plotting


ldn_stations_set <- london_stations[ london_stations$Station %in% c(levels(LDNDataPreds$station)) ,]
ldn_stations_set <-droplevels( ldn_stations_set[order(ldn_stations_set$Latitude),])

lats_of_stations <- aggregate(LDNDataPreds$latitude, by = list(LDNDataPreds$super_station), FUN = mean)
lats_of_stations$x
aggregate(LDNDataPreds$latitude, by = list(LDNDataPreds$super_station), FUN = mean)

### generate a colour platte,one colour per station, colours get lighter as the latitude of the station increases#### 
cols <- colorRampPalette(c("darkgreen", "orange"))
n <- length(levels(LDNDataPreds$super_station))
palette <- cols(n)[as.numeric(cut(lats_of_stations$x,breaks = n))]
ldn_line <- readOGR(dsn ="Geodata/LondonRail",layer = "Tracks")


palette <- palette <-c( "gold2","brown","dodgerblue3","darkorchid4","orangered2","olivedrab2","deeppink2","mediumspringgreen", "gray21","cyan2")
london_stations[,"Station"] <- make.names(london_stations[,"Station"])
install.packages("rworldxtra")

png("LdnPreds.png", width = 12, height = 8, units = 'in', res = 600)
map <- getMap(resolution = "high")
plot(map,xlim = c(-.5,.2), ylim = c(51.42,51.62), col = "grey", bg = "lightskyblue1",border = "darkgrey")
title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
lines(ldn_line, col = "brown3")


for (i in 1:nlevels(LDNDataPreds$super_station)){
  this_station<- levels(LDNDataPreds$super_station)[i]
  find_lats <- LDNDataPreds[LDNDataPreds[,"super_station"] == this_station,]$latPred
  find_longs <- LDNDataPreds[LDNDataPreds[,"super_station"] == this_station,]$longPred
  
  st <- levels(droplevels(LDNDataPreds[LDNDataPreds[,"super_station"] == this_station,]$station))
  points(ldn_stations_set[ldn_stations_set[,"Station"] %in% st,]$Longitude, ldn_stations_set[ldn_stations_set[,"Station"] %in% st,]$Latitude, col = palette[i],bg = palette[i], pch = 17, cex = 1.5)
  #points(LDNDataPreds[LDNDataPreds[,"station"] == this_station,]$longitude, LDNDataPreds[LDNDataPreds[,"station"] == this_station,]$latitude, col = palette[i], pch = 19, cex = 0.6)
  points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.5)
  
}
#legend(-0.465,51.45, legend = c("Stations","Sample predictions"), col = c(palette),pch = c(17,3),box.lty= 1, cex = 1.2, bg = "lightsteelblue1")
axis(1); axis(2)
par(fig = c(0,0.3,0.5,1), new = T) 
plot(map,xlim = c(-1,1), ylim = c(50,53), col = "grey", bg ="lightskyblue1",border = "grey40")
points(LDNDataPreds$city_longitude, LDNDataPreds$city_latitude, cex = 1, pch = 16)
text(LDNDataPreds$city_longitude+1.2,LDNDataPreds$city_latitude, labels = "London", pch = 1.35, font = 2, cex = 0.8)
rect(-0.12574-0.5, 51.50853-0.3 ,-0.12574+0.5, 51.50853+0.3, density = NULL, angle = 45,
     col = NA, border = NULL, pch = 2)
box( col = 'black')

dev.off()



for (i in 1:nrow(LDNDataPreds)){
  LDNDataPreds[i,"Distance_from_origin"] <- distm(c(LDNDataPreds[i,"longPred"],LDNDataPreds[i,"latPred"]), c(LDNDataPreds[i,"longitude"],LDNDataPreds[i,"latitude"]), fun = distHaversine)/1000
}
median(LDNDataPreds$Distance_from_origin)
#Abundance by station plot for fig s13, will be saved to working directory 
levels(LDNDataPreds$super_station) <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10","R11","R12","R13","R14","R15")
ag <-  aggregate(LDNDataPreds[,top_species], by = list(LDNDataPreds$super_station), FUN = median)


for (i in top_species){
  ag[,i] <- (ag[,i] - min(ag[,i]))/(max(ag[,i]) - min(ag[,i]))
}


library(tidyr)
library(hrbrthemes)
library(reshape2)

data.long <- melt(ag)

library(viridis)
ggplot(data = data.long, mapping = aes(x = Group.1,
                                       y = forcats::fct_rev(variable),
                                       fill = value)) +
  geom_tile()+
  scale_fill_viridis(limits = c(0,1)) +
  xlab("Region")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1, size = 11),
        axis.text.y=element_text(size = 12, face = "italic"),
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 12),
        axis.title.x = element_text(size = 12))+
  labs(fill="Relative abundance \n(normalised)")
  #theme(axis.text.x = element_text( hjust = 1))
  

ggsave(
  "LDN_FigS13.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 14,
  height = 8,
  
  dpi = 600,
  limitsize = TRUE,)



























