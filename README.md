# Mapping bacteria - mGPS
Machine learning approach for geographic predictions of environmental samples based on bacterial abundance. 

This repositry contains all code used for data cleaning, modelling and results contained in the mGPS analysis. Outlining the method used for devoping our multi-output model for hierarchical location predictions i.e. (continent/country/city) and location co-ordinates(latitude + longitude) using bacterial abundance data as input. 

## Data

All data required is contained within the `DATA` folder and can also be downloaded
<a id="raw-url" href="https://github.com/leomccarthy1/mGPS/blob/master/DATA/mGPS_data.zip?raw=true">Here</a>
* `Dataset_01_22_2018_enviro.csv` - Soil environmental and geo metadata
* `Dataset_01_22_2018_taxa.csv` - Soil OTU taxa abundance data
* `complete_metadata.csv` - MetaSUB enrionnmental and geo data
* `refseq.krakenhll_species.csv` - MetaSUB species abundance data
* `metasub_data_full.csv` - (cleaned) MetaSUB geo and species data combined by sample ID
* `metasub_200species.csv` - Optimal 200 species for global MetaSUB predictions
* `soil_200_OTU.csv` - Optimal 200 OTU taxa for global soil preidctions


## Usage 

* `mGPS.R` Contains code for application of the main mGPS algorithm the function takes several arguments:   
  - `training` -- Taxa data of samples used to train the model. The chaining structure of the model requires that data for taxa abundance, continent, `classTarget` and co-ordinates are provided for each sample. However only abundance data will be used as predictor variables.  
  - `testing` -- taxa abundance data of samples for which predictions are to be generated  
  - `classTarget` -- granularity for geographic class prediction i.e. country,city or trasit station etc. 
  - `variables` -- a vecttor containing names of species or taca to be used as variables for prediction. This needs definining even if all taxa/species in the training data frame are to be used, so that geographic information are not mistakenly used as predictors. 
  
If no test set is given then a trained model is returned that takes a test set as the input. 
  
For this implementation of mGPS, hyperparameter tuning is carried out at every level of the chained model using a small grid search applied to the training set provided, the same validation fold splits are used at every level. Predictions are then generated using the test data set provided. This also allows for predictions to be generated using different taxa to the ones used in the mGPS paper. 


* `metasub_global.(R)md` is an R notbook explicitly outlining the steps taken from data cleaning through to the mian modelling process and generating results for predictions on a global scale using the MetaSUB data set. This is the best place to start in understanding the modeling workflow implemented here and also contains code for plots for [] manuscript. 

* `soil_global.R` contains code used for applying algorithm to soil data, generating (outer) cross validation predictions and creating plots for manuscript. Includes data cleaning and feature selection




