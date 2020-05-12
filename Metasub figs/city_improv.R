



metasub_data<- read.csv("DATA/msub_meta_taxa.csv")
optimumVarsCity <- read.csv("DATA/opt_species_metasub.csv")
optimumVarsCity <- as.character(optimumVarsCity[,"x"])
#load mGPS model
source(mGPS.R)

#Generate 10 folds
set.seed(123)
folds <-  createFolds(metasub_data$city, k = 10, returnTrain = F)

GeoPreds <- list()
library(doParallel)

#Crate 5 sets of predictions on the same test fold, increasing the number of training folds used. 
registerDoParallel(detectCores()-1)
for (i in c(4,5,6,7,8)){
  training <- metasub_data[unlist(folds[1:i], use.names = F),]
  testing <- metasub_data[unlist(folds[9:10],use.names = F),]
  
  GeoPreds[[i]] <- mGPS(training = training, testing = testing, classTarget = "city",variables = optimumVarsCity)
  
}

test_results <- list()
for (i in c(4,5,6,7,8)){
  test_results[[i]] <- cbind(metasub_data[unlist(folds[9:10],use.names = F),c("city","latitude","longitude")], 
                             "cityPred"= GeoPreds[[i]][[1]], 
                             "latPred" = GeoPreds[[i]][[2]], 
                             "longPred" = GeoPreds[[i]][[3]] )
  
}

#Calculate distance from origin of predictions
for (i in c(4,5,6,7,8) ){
  for(p in 1:length(test_results[[i]]$city)){
  
  test_results[[i]][p,"dist"] <- distm(c(test_results[[i]][p,"longitude"],test_results[[i]][p,"latitude"]), c(test_results[[i]][p,"longPred"],test_results[[i]][p,"latPred"]), fun = distHaversine)/1000
}
}



#Make frame for plotting 

accuracy_improv <- data.frame(row.names = c("40","50","60","70","80")) 
accuracy_improv$prop_dist <- c(mean(test_results[[4]][,"dist"] < 500),
                               mean(test_results[[5]][,"dist"] < 500),
                               mean(test_results[[6]][,"dist"] < 500),
                               mean(test_results[[7]][,"dist"] < 500), 
                               mean(test_results[[8]][,"dist"] < 500))

accuracy_improv$med_dist <-  c(median(test_results[[4]][,"dist"]),
                                 median(test_results[[5]][,"dist"] ),
                                 median(test_results[[6]][,"dist"] ),
                                 median(test_results[[7]][,"dist"] ), 
                                 median(test_results[[8]][,"dist"] ))

  



#Plot#
par(mar = c(5,5,2,5))
acur.bar <- barplot(t(accuracy_improv),ylab = "Predictions within 500km (%)", xlab = )

plot(c(40,50,60,70,80),accuracy_improv$prop_dist*100,type = "o", ylab = "Predictions within 500km (%)", xlab = "Training data size (%)")

par(new=TRUE)
plot(c(40,50,60,70, 80),axes = F,ylab = NA,xlab =NA,accuracy_improv$med_dist,col = "red", type = "o")
axis(side = 4,col.axis = "red", col = "red", cex.axis = 1)
mtext(side = 4,line = 3, "Mean distance from origin (km)", col = "red")





