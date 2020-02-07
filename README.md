# Mapping bacteria - mGPS
Machine learning approach for geographic predictions of environmental samples based on bacterial abundance. 

This repositry contains all code used for data cleaning, modelling and alaysing/plotting the results contained in the mGPS paper. The paper in question outlines the method used for devoping our multi-output model for prediction of location class target (continent/country/city) and location co-ordinates(latitude + longitude) using bacterial abundance data as input. The approach used is a variation on chained regrssors/classifiers (https://www.cs.waikato.ac.nz/~eibe/pubs/chains.pdf) (https://link.springer.com/article/10.1007/s10994-016-5546-z). 

## Data

# Usage 

* `mGPS.R` contains terse code for output of predictions given traning and testing bacterial data.   
  - `trainging` -- bacterial data used to train our model  
  - `testing` -- bacterial data for which predictions are generated  
  - `classTarget` -- label for geographic class prediction i.e. continent,country,city etc
  - `vars` -- a vecttor containing names of vartiables to be used as features for prediction. This needs definining even if all            variables are to be used

* `metasub_global.md` (with Rmd of the same name) is an R notbook explicitly outlining the steps taken from data cleaning through to the mian modelling process and generating results for predictions on a global scale using the metaSub data set. This is definately the best place to start in understanding the modeling workflow implemented here and aims to ease the process of reproducing our results


* `Mplots.R` code used for generating MetaSUB-data plots used in mGPS paper

* `soil_global.R` contains code used for applying our algorithm to soil data, generating (outer) cross validation predictions. 

* `Splots.R` code used for generating soil-data plots used in mGPS paper



