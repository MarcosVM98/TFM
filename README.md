# TFM
This repository contains the contents related to the project of my TFM called "Automated wildfire season detection at a global scale:  Application for the development of a predictive system of fireactivity.". This repository allows full reproducibility of the results we obtained in our project as all the stochastic procedures have a fixed random seed.

It contains the following directories:

-"data" contains all the data used to develop the project, including the one that was generated during the study. It contains 4 directories and a .Rdata object:
	
	+"Fire" contains the burned area, biomes and spatial coordinates data
	+"CPC" contains climate teleconnection indices data
	+"Correlation" contains the results of the correlation study related to unimodal fire seasons
	+"ModelData" contains the results of the models' performance in several validation metrics
	+"fireSeason_def_2.Rdata" is a dataframe containing the coordinates, biome, cluster, and fire season of each pixel

-"scripts" contains the R scripts that have the main functions used in the notebooks:
	
	+"biomes.R" is related to the obtention of the biome_dataframe_masked_df.Rdata (the biomes dataframe)
	+"clustering_functions.r" contains the functions related to clustering
	+"correlation_functions.R" contains the functions related to the correlation study
	+"modelling_functions_new.R" contains the functions related to the development of the predictive models
	+"worldmap.Rdata" is the layout used in the spatial plots

-"notebooks" contains the notebooks where code is executed. There are several notebooks:
	
	+"Definitive_Clustering_v2.ipynb" contains the clustering performance
	+"Correlation_Per75_v2.ipynb" contains the correlation study using the original burned area time series and related to the unimodal fire seasons
	+"Correlation_Per75_with_deltas_v2.ipynb" contains the correlation study using the delta burned area time series and related to the unimodal fire seasons
	+"Correlation_bimodal_main_v2.ipynb" contains the correlation study using the original and the delta burned area time series and related to bimodal's main fire seasons
	+"Correlation_bimodal_secondary_v2.ipynb" contains the correlation study using the original and the delta burned area time series and related to bimodal's secondary fire seasons
	+"Modelling_with_deltas_v2.ipynb" contains the development of the models
	+"Figures.ipynb" contains the code to generate the figures that summarize the project results

-"Figures" contains the figures that summarize the project results in .pdf format
