# mGPS
Machine learning approach for geographic predictions of environmental samples based on bacterial abundance. 

This repositry contains all code used for data cleaning, modelling and results contained in the mGPS analysis. Outlining the method used for devoping our multi-output model for hierarchical location predictions i.e. (continent/country/city) and location co-ordinates(latitude + longitude) using bacterial abundance data as input. 

## Data

### Taxa
All data required is contained within the `DATA` folder and can also be downloaded
<a id="raw-url" href="https://www.kaggle.com/tmdb/tmdb-movie-metadata">Here</a>
* `Dataset_01_22_2018_enviro.csv` - Soil origin metadata
* `Dataset_01_22_2018_taxa.csv` - Soil OTU taxa read data
* `complete_metadata.csv` - MetaSUB enrionnmental and geo data
* `metasub_taxa_abundance.csv` - MetaSUB taxa abundance data
* `opt_taxa_metasub.csv` - Optimal 200 species for global MetaSUB predictions
* `opt_taxa_soil.csv` - Optimal 200 OTU taxa for global soil preidctions

### Geodata 
All the geodata required for geographic plots in the mGPS paper. <a id="raw-url" href="https://github.com/leomccarthy1/mGPS/blob/master/Geodata/Geodata.zip?raw=true">Download</a>

## Usage 

* `mGPS.R` Contains code for application of the main mGPS algorithm the function takes several arguments:   
  - `training` -- Taxa data of samples used to train the model. The chaining structure of the model requires that data for taxa abundance, continent, `classTarget` and co-ordinates are provided for each sample. However only abundance data will be used as predictor variables.  
  - `testing` -- taxa abundance data of samples for which predictions are to be generated  
  - `classTarget` -- granularity for geographic class prediction either country,city or trasit station etc. 
  - `variables` -- a vector containing names of species or taxa to be used as variables for prediction. This needs definining even if all taxa/species in the training data frame are to be used, so that geographic information are not mistakenly used as predictors. This does however mean that training can take some time and parallisation through the <a id="raw-url" href="https://cran.r-project.org/web/packages/doParallel/doParallel.pdf">doParallel</a> package is recommended. 
  
If no test set is given then a trained model is returned that takes a test set as the input. 
  
For this implementation of mGPS, hyperparameter tuning is carried out at every level of the chained model using a small grid search applied to the training set provided, the same validation fold splits are used at every level. Predictions are then generated using the test data set provided. This also allows for predictions to be generated using different taxa to the ones used in the mGPS paper. 
 
 
## Cross validation, results and figures
* `metasub_figs`
  - `metasub_global.(R)md` is an R notbook explicitly outlining the dat cleansing and generation of results and figures for the mGPS paper. This is the best place to start in understanding the modeling workflow implemented here and also contains code for plots for [Fig1, FigS1, FigS2, FigS3,FigS9, FigS10,] manuscript. 
  - `HK.R` Hong Kong results and figures[Fig2,FigS4]
  - `NY.R` New York results and figures [Fig3,FigS5]
  - `LDN.R` London results and figures
  - `city_imrprov.R` - Fig S7
* `soil_biome_figs`
  - `soil_biome.Rmd` Soil microbiome results and figures. [Fig4, FigS8, Fig14,Figs15]
  
  










