library(sp)
library(mapdataNE)
library(rworldmap)
library(nominatim)
library(geonames)
library(caret)
library(rpart)
library(maps)
library(MASS)
library(randomForest)
library(geosphere)
library(revgeo)
library(ggmap)
library(ggplot2)
library(plyr)
library(party)
library(Boruta)
library(dplyr)
library(randomForestSRC)

#######reading and tidying data######
city <- read.csv(file = "complete_metadata.csv", header = TRUE)
####remove cases for which city isnt known, this is a quick way to remove controls######
city <- city[-c(4488:4493), ] 
city <- city[complete.cases(city[ , 19:20]),]
bacdata <- read.csv(file= "refseq.krakenhll_species (1).csv", header = TRUE)
bacdata[is.na(bacdata)] <- 0   ##turning NAs to 0 for abundance
#### merge the 2 data sets for ease 
total2 <- merge(city,bacdata,by.x="uuid",by.y="uuid")
#####set lat long as city lat long for those where exact isnt known
total2[is.na(total2[,"latitude"]), "latitude"] <- total2[is.na(total2[,"latitude"]), "city_latitude"]
total2[is.na(total2[,"longitude"]), "longitude"] <- total2[is.na(total2[,"longitude"]), "city_longitude"]
total2 <- total2[complete.cases(total2[ ,7:8]),]
 #### adjusting known incorrect values in data set
total2$latitude[total2$city == "kyiv"] <- total2$city_latitude[total2$city == "kyiv"]
total2$longitude[total2$city == "kyiv"] <- total2$city_longitude[total2$city == "kyiv"]
total2$continent[total2$city == "porto"] <- "europe"
total2 <- droplevels(total2)
total2$city <- as.factor(total2$city)
levels(starting_set$city) <- c("Auckland" ,"Baltimore" ,"Barcelona", "Berlin","Bogota", 
                         "Brisbane","Denver","Doha" ,"Fairbanks","Hamilton",
                         "Hanoi","Hong_Kong","Ilorin","Kuala_Lumpur","Kyiv",
                         "Lisbon", "London" , "Minneapolis", "Naples"  ,"New_York_City",
                         "Offa" ,"Oslo" ,"Paris","Porto","Rio_de_Janeiro",
                         "Sacramento","San_Francisco" ,"Santiago" , "Sao_Paulo" ,"Sendai",
                         "Seoul","Singapore", "Sofia","Stockholm","Taipei",
                         "Tokyo","Vienna","Yamaguchi","Zurich")

#### generate feautre importance based on city classification random forest######
features_forest_city <- randomForest(city~., data = total2[,c(5,40:ncol(total2))], ntree = 300, importance = T)
feature_importances <- randomForest::importance(features_forest_city, type = 1, scale = F)

selected_species <- order(-feature_importances)
starting_set <-total2[,c(1:39,selected_species[1:500]+39)]
total2 <- starting_set




png("Feature_importance_distribution.png", width = 13, height = 8, units = 'in', res = 600)
plot(randomForest::importance(features_forest_city, scale = F, type = 1),ylab = "Species importance (Mean decrease in accuracy)", xlab = "Species Index")
dev.off()

png("Feature_importance.png", width = 13, height = 8, units = 'in', res = 600)
imp <- feature_importances[order(-feature_importances), , drop = FALSE]
dotchart(rev(imp[1:25]),labels=rev(row.names(imp[1:25, ,drop = F])),cex=1,pt.cex = 1.2,
         xlab="Mean Decrease In Accuracy")
dev.off()

#### build random fores on full set to predict continent and get probabilities of each continent, add these to data set 
original_continent_forest <- randomForest(continent~., data = total2[,c(27,40:ncol(total2))])
continent.probs <- predict(original_continent_forest, type = "prob")
total2 <- cbind(total2,continent.probs)
##build random forest to predict city, find probabilities of each city and add to data set#### 
total2$city <- droplevels(total2$city)
original_city_forest <- randomForest(city~., data =total2[,c(5,40:ncol(total2))])
original_city_forest 
city.probs <-  predict(original_city_forest, type = "prob")
total2 <-cbind(total2, city.probs)



#### final latitude and longitude random forests 
original_latitude_forest <- randomForest(latitude~., data = total2[,c(7,40:ncol(total2))])
pred_oob_lats <- predict(original_latitude_forest)
total2$predicted_lat <- predict(original_latitude_forest)

original_longitude_forest <- randomForest(longitude~., data = total2[,c(8,40:ncol(total2))])
pred_oob_longs <- predict(original_longitude_forest)
total2$predicted_long <- predict(original_longitude_forest)



### find distance from origin using haversine formula ####
distance1 <- c()
lat.pred <- total2$predicted_lat
long.pred <- total2$predicted_long
t1 <- total2[,"city_longitude"]
t2 <- total2[,"city_latitude"]
for (i in 1:nrow(total2)){
  distance1[i] <- distm(c(long.pred[i],lat.pred[i]), c(t1[i],t2[i]), fun = distHaversine)/1000
  
}
mean(distance1)

total2[,"Distances_in_km"] <- distance1


#####pull too coast#####

coastlines <- SpatialLines2map(coastsCoarse)
coastlines <- cbind(coastlines$x, coastlines$y)
coastlines <- coastlines[complete.cases(coastlines),]
coastlines <- coastlines[coastlines[,1] < 180 ,]
find_coast <- 
  function(long,lat){
    distances_from_coastline <-  distm(cbind(coastlines[,1],coastlines[,2]), cbind(rep(long,length(coastlines[,1])),rep(lat,length(coastlines[,2]))), fun = distHaversine)
    closest_point <-  which.min(distances_from_coastline)
    new_coords <- coastlines[closest_point,]
    
    return(new_coords)
  }


pred_countrires1 <- map.where(database = "mapdataNE::world10", total2$predicted_long, total2$predicted_lat)
pred_countrires1[is.na(pred_countrires1)] <-"not found"
total2[,"predicted_country"] <- factor(make.names(pred_countrires1))

in_sea1 <- total2[total2$predicted_country == "not.found",]
found1 <- mapply(find_coast, long = in_sea1$predicted_long, lat = in_sea1$predicted_lat)
total2[total2$predicted_country == "not.found", "predicted_lat"] <- found1[2,]
total2[total2$predicted_country == "not.found", "predicted_long"] <- found1[1,]




##### training and testing #####



#### generate test sample of size 20% of each city population

test_sample <- c()
for (i in 1:length(levels(starting_set$city))){
  city_to_sample <- levels(starting_set$city)[i]
  set_to_sample <- starting_set[starting_set[,"city"] == city_to_sample,]
  test_sample <- c(test_sample,sample(as.numeric(row.names(set_to_sample)), trunc(0.25 * nrow(set_to_sample))))
  
}
test_set <- starting_set[c(test_sample) ,]
train_set <- starting_set[-c(test_sample) ,]



#### build random fores on full set to predict continent and get probabilities of each continent, add these to data set 
train_set$continent <- droplevels(train_set$continent)
train_continent_forest <- randomForest(continent~., data = train_set[,c(27,40:ncol(train_set))], ntree = 300)
train_continent.probs <- predict(train_continent_forest, type = "prob")
train_set <- cbind(train_set,train_continent.probs)

#build random forest to predict city, find probabilities of each city and add to data set 
train_set$city <- droplevels(train_set$city)
train_city_forest <- randomForest(city~., data =train_set[,c(5,40:ncol(train_set))], ntree =300)
train_city.probs <-  predict(train_city_forest, type = "prob")
train_set <-cbind(train_set, train_city.probs)


#### final latitude and longitude random forests 
train_latitude_forest <- randomForest(latitude~., data = train_set[,c(7,40:ncol(train_set))], ntree = 300)
train_longitude_forest <- randomForest(longitude~., data = train_set[,c(8,40:ncol(train_set))], ntree=300)




#### now apply trained models to test set 



Model_1 <- function(bacdata1){
   
    
  cont_prob <- predict(train_continent_forest, newdata = bacdata1, type = "prob")
  bacdata1 <- cbind(bacdata1,cont_prob)
  
  city_prob <- predict(train_city_forest, newdata = bacdata1, type = "prob")
  bacdata1 <- cbind(bacdata1,city_prob)
  
  
  
  lat <- predict(train_latitude_forest, newdata = bacdata1)
  long <- predict(train_longitude_forest, newdata = bacdata1)
  
  
  
  
return(list("lat" = lat,"long" = long))
}
test_results <- Model_1(test_set[,40:ncol(test_set)])
test_set <- droplevels(test_set)

distance2 <- c()
lat.pred <- test_results$lat
long.pred <- test_results$long
t1 <- test_set[,"city_longitude"]
t2 <- test_set[,"city_latitude"]
for (i in 1:nrow(test_set)){
  distance2[i] <- distm(c(long.pred[i],lat.pred[i]), c(t1[i],t2[i]), fun = distHaversine)/1000
}
mean(distance2)

test_set[,"Distances_in_km"] <- distance2

bar_df2 <- data.frame(row.names = c("Overall",levels(test_set$city)))

for (i in 1: length(levels(test_set$city))){
  this_city <- levels(test_set$city)[i]
  prop <- mean(test_set[test_set$city == this_city,][,"Distances_in_km"] < 100)
  bar_df2[i+1,"0 - 100km"] <- prop
  
  overall_prop <- mean(test_set[,"Distances_in_km"] < 100)
  bar_df2[ 1,"0 - 100km"] <- overall_prop
}
for (i in 1: length(levels(test_set$city))){
  this_city <- levels(test_set$city)[i]
  prop <- mean(test_set[test_set$city == this_city,][,"Distances_in_km"] > 100 & test_set[test_set$city == this_city,][,"Distances_in_km"] < 500)
  bar_df2[i+1,"100 - 500km"] <- prop
  
  overall_prop <- mean(test_set[,"Distances_in_km"] > 100 & test_set[,"Distances_in_km"] < 500)
  bar_df2[ 1,"100 - 500km"] <- overall_prop
}



for (i in 1: length(levels(test_set$city))){
  this_city <- levels(test_set$city)[i]
  prop <- mean(test_set[test_set$city == this_city,][,"Distances_in_km"] > 500 & test_set[test_set$city == this_city,][,"Distances_in_km"] < 1000)
  bar_df2[i+1,"500 - 1000km"] <- prop
  
  overall_prop <- mean(test_set[,"Distances_in_km"] > 500 & test_set[,"Distances_in_km"] < 1000)
  bar_df2[ 1,"500 - 1000km"] <- overall_prop
}

for (i in 1: length(levels(test_set$city))){
  this_city <- levels(test_set$city)[i]
  prop <- mean(test_set[test_set$city == this_city,][,"Distances_in_km"] > 1000 & test_set[test_set$city == this_city,][,"Distances_in_km"] < 2000)
  bar_df2[i+1,"1000 - 2000km"] <- prop
  
  overall_prop <- mean(test_set[,"Distances_in_km"] > 1000 & test_set[,"Distances_in_km"] < 2000)
  bar_df2[ 1,"1000 - 2000km"] <- overall_prop
}
for (i in 1: length(levels(test_set$city))){
  this_city <- levels(test_set$city)[i]
  prop <- mean(test_set[test_set$city == this_city,][,"Distances_in_km"] > 2000 & test_set[test_set$city == this_city,][,"Distances_in_km"] < 3000)
  bar_df2[i+1,"2000 - 3000km"] <- prop
  
  overall_prop <- mean(test_set[,"Distances_in_km"] > 2000 & test_set[,"Distances_in_km"] < 3000)
  bar_df2[ 1,"2000 - 3000km"] <- overall_prop
}
for (i in 1: length(levels(test_set$city))){
  this_city <- levels(test_set$city)[i]
  prop <- mean(test_set[test_set$city == this_city,][,"Distances_in_km"] > 3000 )
  bar_df2[i+1,">3000km"] <- prop
  
  overall_prop <- mean(test_set[,"Distances_in_km"] > 3000)
  bar_df2[ 1,">3000km"] <- overall_prop
}
size2 <- c()
for (i in 1: length(levels(test_set$city))){
  
  this_city <- levels(test_set$city)[i]
  size2[i] <- length(which(test_set$city == this_city))
}

levels(test_set$city) <- levels(total2$city)
png("Test_bar.png", width = 12,height = 7, units = 'in', res = 600)
par(xpd = T, mar = par()$mar + c(1,0,0,7))
bp2 <- barplot(t(bar_df2), col=c("lightyellow1","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), names.arg=c("Overall",paste0(levels(test_set$city),"  (",size2,")")) ,args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.55)

legend("topright",inset = c(-0.15,0.4), rev(c(colnames(bar_df2))), fill = rev(c("lightyellow1","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 0.8)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()






