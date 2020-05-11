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

* `mGPS.R` contains code for application of the main mGPS algorithm the function takes several arguments:   
  - `training` -- bacterial data used to train our model  
  - `testing` -- bacterial data for which predictions are generated  
  - `classTarget` -- granularity for geographic class prediction i.e. country,city or trasit station etc
  - `variables` -- a vecttor containing names of vartiables to be used as features for prediction. This needs definining even if all variables are to be used.
  
  If no training set is given then predictions are generated using our pretrained model. If no test set is given then a trained model is returned that takes a test set as the input. 
  
  For predictions using a pre-trained model the "trained" folder must be contained in teh working directory.
  
  Predictions for the paper were generated using a nested 5-fold corss validation procedure, where hyperparameter tuning was carried out using inner cross validation on each training fold independently. 

* `metasub_global.(R)md` is an R notbook explicitly outlining the steps taken from data cleaning through to the mian modelling process and generating results for predictions on a global scale using the MetaSUB data set. This is the best place to start in understanding the modeling workflow implemented here and also contains code for plots for [] manuscript. 

* `soil_global.R` contains code used for applying algorithm to soil data, generating (outer) cross validation predictions and creating plots for manuscript. Includes data cleaning and feature selection




