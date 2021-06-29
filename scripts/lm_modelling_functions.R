# Packages required
library(tidyr)
library(dplyr)
require(sp)
require(magrittr)
library(caret)
library(randomForest)


# Data required
load('fireSeasonPer75_def.Rdata', verbose = T)
load("data/ba_mon_time_series_masked.Rdata", verbose = T)
load("data/ba_mon_clim_masked_df.Rdata", verbose = T)
load('corrDfAnnualClus.Rdata', verbose = T)


#' @title Burned area time series calculation
#' @description Obtaining the annual burned area time series as the sum of the burned area during the fire season in each year of all the points in the cluster
#' @param m Subset of the burned area dataframe containing only the burned area data in the points of the cluster during the months of the fire season
#' @param meses Vector containing the months of the fire season of that cluster
#' @param ind.meses Vector containing the positions in the burned area dataframe of the months of the fire seasons during the period 2001-2020
#' @param ind.coords Vector containing the position of the points of the cluster in masked_coords dataframe
#' @return vector with the burned area time series
get.ba.serie <- function(m, meses, ind.meses, ind.coords){    
    if (length(meses) == 1){
        ba.serie = m[,1]
        for (j in 2:length(ind.coords)){                    
            ba.serie = ba.serie + m[,j]
        }
    } else {
        ba.serie = c()
        for (i in 1:(length(ind.meses)/length(meses))){
            ba.serie = c(ba.serie, sum(m[((i-1)*(length(meses))+1):(i*length(meses)),1]))
        }
        for (j in 2:length(ind.coords)){
            l = c()
            for (i in 1:(length(ind.meses)/length(meses))){
                l = c(l, sum(m[((i-1)*(length(meses))+1):(i*length(meses)),j]))
            }
            ba.serie = ba.serie + l
        }
    }        
    return (ba.serie)
}


#' @title Climate index time series calculation
#' @description Obtaining the annual climate index time series of the cpcs given. The value of each year is the average of the values of the index during the fire seasons months
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ind.coords Vector containing the position of the points of the cluster in masked_coords dataframe
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different data frame
#' @param meses Vector containing the months of the fire season of that cluster
#' @param ind.meses Vector containing the positions in the masked_ba_series.log array of the months of the fire seasons during the period 2001-2020
#' @param df dataframe whose new columns will be the climate indexes time series
#' @return dataframe whose new columns will be the climate indexes time series
get.cpc.serie <- function(fireSeasons, ind.coords, list.cpcs, meses, ind.meses, df){            
    # If the fire season takes place in only one year 
    if (fireSeasons[ind.coords,]$start.1[1] <= fireSeasons[ind.coords,]$end.1[1]){
        ind.years = as.numeric(substr(dates[ind.meses], 1, 4))

        if (length(meses) == 1){
            for (i in 1:length(list.cpcs)){
                cpc.serie = list.cpcs[[i]][which(list.cpcs[[i]]$Year %in% ind.years), c(meses + 1)]
                df = cbind(df, cpc.serie)
            }                    
        } else {
            for (i in 1:length(list.cpcs)){                    
                n = list.cpcs[[i]][which(list.cpcs[[i]]$Year %in% ind.years), c(meses + 1)]
                cpc.serie = apply(n, 1, mean)
                df = cbind(df, cpc.serie)
            }
        }                

    # If the fire season takes place in two different years
    } else {                
        ind.years = as.numeric(substr(dates[ind.meses], 1, 4))

        for (i in 1:length(list.cpcs)){
            n = list.cpcs[[i]][which(list.cpcs[[i]]$Year %in% ind.years), c(meses + 1)]
            n = as.vector(t(n))

            cpc.serie = c()
            for (j in 1:(length(ind.meses)/length(meses))){
                cpc.serie = c(cpc.serie, mean(n[(fireSeasons[ind.coords,]$end.1[1] + 1 + (j-1)*length(meses)):
                                (fireSeasons[ind.coords,]$end.1[1] + j*length(meses))]))
            }
            df = cbind(df, cpc.serie)
        }
    }                        
    return (df)
}


#' @title Climate indexes and burned area time series calculation for a given cluster
#' @description Obtaining the annual climate index time series of the cpcs given and the burned area time serie. The value of each year of the cimate index time serie is the average of the values of the index during the fire seasons months. The annual burned area time series is obtained as the sum of the burned area during the fire season in each year of all the points in the cluster
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different data frame
#' @param biome Number between 1 and 13 representing a biome
#' @param cluster Number of the cluster whose burned area time series are we going to obtain
#' @param pvalue Threshold for considering significant each correlation
#' @return list containing in the element 'df' a dataframe whose first column will be the burned area time serie and the other ones will be the climate indexes time series. The second element is a vector containing the climate indexes with significant correlation
lm.clus.data.preparation <- function(fireSeasons, ba.series, corr.df, list.cpcs, biome, cluster, pvalue = 0.05){
    
    # Positions of the cluster data points in the dataframes
    ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == cluster)

    # Indexes with significant correlation with this cluster
    cpcs = c()
    for (i in 1:length(list.cpcs)){
        if (!is.na(corr.df[ind.coords, 2*i+2][1]) && corr.df[ind.coords, 2*i+2][1] < pvalue){
            cpcs = c(cpcs, i)
        }
    }

    cat('Correlated indexes:', cpcs, '\n')

    # Positions of the months of the fire season in 'dates' dataframe
    # If the fire season takes place in only one year 
    if (fireSeasons[ind.coords,]$start.1[1] <= fireSeasons[ind.coords,]$end.1[1]){
        meses = seq(fireSeasons[ind.coords,]$start.1[1], fireSeasons[ind.coords,]$end.1[1])
        ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)

        ind.meses = ind.meses[1:(length(meses) * floor(length(ind.meses)/length(meses)))]

    # If the fire season takes place in two different years
    } else {                
        meses = c(seq(1, fireSeasons[ind.coords,]$end.1[1]), seq(fireSeasons[ind.coords,]$start.1[1], 12))                
        ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)

        if (fireSeasons[ind.coords,]$end.1[1] <= 4){
            ind.meses = ind.meses[(fireSeasons[ind.coords,]$end.1[1]+1):(length(ind.meses))]
        } else {
            ind.meses = ind.meses[(fireSeasons[ind.coords,]$end.1[1]+1):(length(ind.meses)
                                                                         -4-(13-fireSeasons[ind.coords,]$start.1[1]))]
        }
    }

    # Subset containing the data in the points of the cluster during the fire season months
    m = ba.series[ind.meses, ind.coords]

    # Obtaining burned area time series as the sum of the burned area during the fire season of all the points in the cluster
    ba.serie <- get.ba.serie(m, meses, ind.meses, ind.coords)                
    df = data.frame('ba' = ba.serie)

    # Obtaining the climate indexes time series
    df <- get.cpc.serie(fireSeasons, ind.coords, list.cpcs, meses, ind.meses, df)
    colnames(df) = c('ba', (1:length(list.cpcs)))
            
    return (list(df = df, cpcs = cpcs))
}


#' @title Model's validation
#' @description Obtains the value of some metrics for the results obtained by a model
#' @param obs Vector with the real values
#' @param pred Vector with the predicted values
#' @returns dataframe containing the value of the RMSE, R2 coefficient and variance and 90 percentile ratios
lm.clus.model.validation <- function(obs, pred){
    
    rss = sum((pred - obs) ^ 2)  ## residual sum of squares
    tss = sum((obs - mean(obs)) ^ 2)  ## total sum of squares
    R2 = 1 - rss/tss
    RMSE = sqrt(rss/length(pred))
    RVar = var(pred)/var(obs)
    Rp90 = quantile(pred, prob = 0.9) / quantile(obs, prob = 0.9)
    
    results = data.frame(RMSE = RMSE, R2 = R2, RVar = RVar, Rp90 = Rp90)
    
    return (results)    
}




#' @title Linear model obtention and validation
#' @description Obtaining the linear model and evaluating the results. It also plots the real and the predicted values of the burned area, the burned area versus the climate index with significant correlation if there is only one and climate index time serie
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different dataframe
#' @param biome Number between 1 and 13 representing a biome
#' @param cluster Number of the cluster whose burned area time series are we going to obtain
#' @param pvalue Threshold for considering significant each correlation
#' @return the linear model previously calculated
lm.clus <- function(fireSeasons, ba.series, corr.df, list.cpcs, biome, cluster, pvalue = 0.05){
    # Obtaining the time series
    data = lm.clus.data.preparation(fireSeasons, ba.series, corr.df, list.cpcs, biome, cluster, pvalue)
    df = data$df
    cpcs = data$cpcs
    
    # Training the linear model
    ctrl <- trainControl(method = "LOOCV")# Leave-one-out CV
    mod <- train(ba ~ ., data = df[,c(1, c(cpcs + 1))], method = "lm", trControl = ctrl)
    
    # Validating the model
    results = lm.clus.model.validation(mod$pred$obs, mod$pred$pred)
    
    # Doing the plots
    par(mfrow=c(2,2))
    
    # Observed vs predicted plot
    plot(mod$pred$obs, col = 'green', main = 'Observed vs predicted', ylab = 'Burned area')
    points(mod$pred$pred, col = 'red')
    
    # Climate time series plots
    for (i in 1:length(cpcs)){
        plot(df[,cpcs[i]+1], main = paste('Climate index', toString(cpcs[i])))
    }
    
    # Climate index vs burned area
    if (length(cpcs) == 1){
        ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == cluster)
        plot(df[,cpcs[i]+1], df$ba, main = paste('Climate index vs ba with corr = ', 
                                          toString(round(corr.df[ind.coords, 2*cpcs + 1][1], digits = 3))), 
             ylab = 'Burned area')
    }
    
    # Printing the validation results
    print(results)
    
    return (mod)
}