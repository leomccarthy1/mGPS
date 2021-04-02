# mGPS
**Microbiome biodiversity typify fine-scale biogeography**

*Leo McCarthy and Eran Elhaik*


This repositry contains all code used for data cleaning, modelling and results contained in the mGPS analysis. Outlining the method used for devoping our multi-output model for hierarchical location predictions i.e. (continent/country/city) and location co-ordinates(latitude + longitude) using bacterial abundance data as input. 

## Data

### Taxa
All dat and geodata required to reproduce the analysis is contained within the `Data` folder and needs to be unzipped before running the analysis. 
* `Dataset_01_22_2018_enviro.csv` - Soil origin metadata
* `Dataset_01_22_2018_taxa.csv` - Soil OTU taxa read data
* `complete_metadata.csv` - MetaSUB enrionnmental and geo data
* `metasub_taxa_abundance.csv` - MetaSUB taxa abundance data
* `opt_taxa_metasub.csv` - Optimal 200 species for global MetaSUB predictions
* `opt_taxa_soil.csv` - Optimal 200 OTU taxa for global soil preidctions



## Usage 

* `mGPS()` function takes several arguments:   
  - `training` -- Taxa data of samples used to train the model. 
  - `testing` -- taxa abundance data of samples for which predictions are to be generated  
  - `classTarget` -- granularity for geographic class prediction either country,city or trasit station etc. 
  - `variables` -- a vector containing names of species or taxa to be used as variables for prediction. This needs definining even if all taxa/species in the training data frame are to be used, so that geographic information are not mistakenly used as predictors. 
  - `hierarchy` -- The geographic hierarchy for predictions i.e. 
  - `nthread` -- number of threads to utilise 
  - `coast` -- (optional) data.frame of co-ordinates for predictions to be bound by
  
If no test set is given then a trained model is returned that takes a test set as the input. 
  
For this implementation of mGPS, hyperparameter tuning is carried out at every level of the chained model using a small grid search applied to the training set provided, the same validation fold splits are used at every level. Predictions are then generated using the test data set provided. 

```R
mGPS(training = train, 
     testing = test, 
     classTarget = "city",
     variables = optVars,
     nthread = 8,
     hierarchy = c('continent','city','latitude','longitude'), 
     coast=coastlines)
```

* `species_select()` - function for selecting optimal Geographically informative taxa:
  - `x` -- predictior variables
  - `y` -- target variable 
  - `remove_correlated` -- should correlated predictor variables be removed (>98% correlation). 
  - `subsets` -- The variable subset sizes to try 
  - `cores` -- number of cores to utilise 

```R
species_select(x = metasub_data[, taxa],
               y = metasub_data$city,
               remove_correlated = F,
               c(50,100,200,300,500),
               cores = 8)
```


 
## Results and figures

Each data set used in the analysis has its own folder here. Each dataset has a `{dataset}_make.R` file for generating predictions using `mGPS` and cross validation. Also there is a `{dataset}_plots.Rmd` file for each dataset in the analysis which produces plots and tables found in the manuscript.

Each `make` and `plots` file for the corresponding datasets should be run in turn. After cloning the repo the `plots` files can be run without first running the corresponding `make` fiule as the predictions are stored in this repo. 
  
  










