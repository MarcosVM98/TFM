# Functions to build and test the linear models


# Packages required
library(tidyr)
library(dplyr)
require(sp)
require(magrittr)
library(caret)
library(Matrix)


#' @title Burned area time series calculation
#' @description Obtaining the annual burned area time series as the sum of the burned area during the fire season in each year of all the points in the cluster
#' @param m Subset of the masked_ba_series.log dataframe containing only the burned area data in the points of the cluster during the months of the fire season
#' @param months Vector containing the months of the fire season of that cluster
#' @param ind.months Vector containing the positions in the masked_ba_series.log dataframe of the months of the fire seasons during the period 2001-2020
#' @param ind.coords Vector containing the position of the points of the cluster in masked_coords dataframe
#' @return vector with the burned area time series
get.ba.serie.detrended <- function(m, months, ind.months, ind.coords){    
    if (length(months) == 1){
        ba.serie = m[,1]
        for (j in 2:length(ind.coords)){                    
            ba.serie = ba.serie + m[,j]
        }
    } else {
        ba.serie = c()
        for (i in 1:(length(ind.months)/length(months))){
            ba.serie = c(ba.serie, sum(m[((i-1)*(length(months))+1):(i*length(months)),1]))
        }
        for (j in 2:length(ind.coords)){
            l = c()
            for (i in 1:(length(ind.months)/length(months))){
                l = c(l, sum(m[((i-1)*(length(months))+1):(i*length(months)),j]))
            }
            ba.serie = ba.serie + l
        }
        ind = 1:length(ba.serie)
        dff = data.frame(ind, ba.serie)
        ba.serie = lm(ba.serie ~ ind, data = dff)$residuals
    }        
    return (ba.serie)
}



#' @title Climate index time series calculation
#' @description Obtaining the annual climate index time series of the cpcs given. The value of each year is the average of the value of the index during the fire seasons months
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param dates Array containing the dates of each burned area observation
#' @param ind.coords Vector containing the position of the points of the cluster in masked_coords dataframe
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different data frame
#' @param ind.start Vector containing the positions in the masked_ba_series.log array of the starting month of the fire season
#' @param df dataframe whose new columns will be the climate indices time series
#' @param lt ...
#' @param duration Number of months considered each year for the climate time series calculation
#' @return dataframe whose new columns will be the climate indices time series
get.cpc.serie.detrended <- function(fireSeasons, dates, ind.coords, list.cpcs, ind.start, df, lt, duration = 1){
    
    # Update the indices
    ind.start = ind.start - lt
    ind.months = ind.start
    if (duration > 1){
        ind.months.0 = ind.months
        for (i in 2:duration){
            ind.months = c(ind.months, ind.months.0 - i + 1)
        }  
    }
    ind.months = sort(ind.months)
    
    if (ind.months[1] < 1){
        ind.years = as.numeric(substr(dates[ind.months[which(ind.months > 0)]], 1, 4))
        ind.years = c(min(ind.years)-1, ind.years)
        ind.months = ind.months + 12
    } else {
        ind.years = as.numeric(substr(dates[ind.months], 1, 4))
    }
    ind.months = unique(ind.months)
    ind.years = unique(ind.years)
    
    for (i in 1:length(list.cpcs)){
        n = list.cpcs[[i]][which(list.cpcs[[i]]$Year %in% ind.years), -1]
        n = as.vector(t(n))
        n = n[ind.months]

        cpc.serie = c()
        for (i in 1:(length(ind.months)/duration)){
            cpc.serie = c(cpc.serie, mean(n[(duration*(i-1)+1):(duration*i)]))
        }
        
        ind = 1:length(cpc.serie)
        dff = data.frame(ind, cpc.serie)
        cpc.serie = lm(cpc.serie ~ ind, data = dff)$residuals
        df = cbind(df, cpc.serie)
    }
                     
    return (df)
}




#' @title Climate indices and burned area time series calculation for a given cluster
#' @description Obtaining the annual climate index time series of the cpcs given and the burned area time serie. The value of each year of the cimate index time series are the residuals of the average of the values of the index during the fire seasons months. The annual burned area time series is obtained as the rediduals of the sum of the burned area during the fire season in each year of all the points in the cluster
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param coords Array containing the x and y coordinates of each point of the ba.series dataframe
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between burned area during fire season months in each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indices in the period that we are studying. Data of each index must be in a different data frame
#' @param biome Number between 1 and 13 representing a biome
#' @param cluster Number of the cluster whose burned area time series are we going to obtain
#' @param lt ...
#' @param duration Number of months considered each year for the climate time series calculation
#' @param pvalue Threshold for considering significant each correlation. Default to 0.05
#' @return A copy of corr.df dataframe with a new column with the correlation
clus.data.preparation <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, lt, duration = 1, pvalue = 0.05){   
         
    ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == cluster)
    ind.start = which(as.numeric(substr(dates, 6, 7)) %in% fireSeasons[ind.coords,]$start.1[1])
            
    if (fireSeasons[ind.coords,]$start.1[1] <= fireSeasons[ind.coords,]$end.1[1]){
        meses = seq(fireSeasons[ind.coords,]$start.1[1], fireSeasons[ind.coords,]$end.1[1])
        ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)
        if (fireSeasons[ind.coords,]$start.1[1] <= 4 & fireSeasons[ind.coords,]$end.1[1] > 4){# because the last month with burned area data is April
            ind.start = ind.start[-length(ind.start)]
        }
    } else {
        meses = c(seq(1, fireSeasons[ind.coords,]$end.1[1]), seq(fireSeasons[ind.coords,]$start.1[1], 12))                
        ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)

        if (fireSeasons[ind.coords,]$end.1[1] <= 4){# because the last month with burned area data is April
            ind.meses = ind.meses[(fireSeasons[ind.coords,]$end.1[1]+1):(length(ind.meses))]
            if (fireSeasons[ind.coords,]$start.1[1] <= 4){
                ind.start = ind.start[-length(ind.start)]
            }
        } else {
            ind.meses = ind.meses[(fireSeasons[ind.coords,]$end.1[1]+1):(length(ind.meses)-4-(13-fireSeasons[ind.coords,]$start.1[1]))]
            ind.start = ind.start[-length(ind.start)]
        }
    }
    ind.meses = ind.meses[1:(length(meses) * floor(length(ind.meses)/length(meses)))]

    m = ba.series[ind.meses, ind.coords]                

    # Obtaining burned area time series as the sum of the burned area during the fire season of all the points in the cluster
    ba.serie <- get.ba.serie.detrended(m, meses, ind.meses, ind.coords)  
    df = as.data.frame(ba.serie)

    # Obtaining climate indices time series
    for (i in 1:length(list.cpcs)){
        df0 <- get.cpc.serie.detrended(fireSeasons, dates, ind.coords, list(list.cpcs[[i]]), ind.start, df, lt, duration)
        if (corr.df[ind.coords[1],2*i+2] < pvalue){
            colnames(df0)[dim(df0)[2]] = paste('P', i, 'LT', lt, 'D', duration, sep='')
            df = df0
        }
    }
    
    return (df)
}



#' @title Model validation
#' @description Obtains the value of some metrics for the results obtained by a model
#' @param obs Vector of real values
#' @param pred Vector of predicted values
#' @param as.df Decide whether the results are a vector or a data frame. Default to FALSE
#' @return dataframe or vector containing the value of the RMSE, bias, correlation, variance ratio, total accuracy and tercile accuracy
model.validation <- function(obs, pred, as.df = F){
    
    rss = sum((pred - obs) ^ 2)  ## residual sum of squares
    RMSE = 100 * sqrt(rss/length(pred)) / sd(obs)
    bias = 100*(mean(pred) - mean(obs)) / sd(obs)
    cor = cor.test(obs, pred, method = 'pearson')
    RVar = var(pred) / var(obs)
    
    t1 = quantile(obs, prob = 1/3)
    t2 = quantile(obs, prob = 2/3)
    
    obs_ter = obs
    obs_ter[which(obs < t1)] <- 1
    obs_ter[which(obs >= t1)] <- 2
    obs_ter[which(obs >= t2)] <- 3

    t1 = quantile(pred, prob = 1/3)
    t2 = quantile(pred, prob = 2/3)
    pred_ter = pred
    pred_ter[which(pred < t1)] <- 1
    pred_ter[which(pred >= t1)] <- 2
    pred_ter[which(pred >= t2)] <- 3
    
    acc.total = sum(pred_ter == obs_ter) / length(obs_ter)
    acc = c()
    for (i in 1:3){
        acc[i] = sum(pred_ter[which(obs_ter == i)] == i) / sum(obs_ter == i)
    }
    
    if (as.df == T){
        results = data.frame(RMSE = RMSE, bias = bias, RVar = RVar, cor.pvalue = cor$p.value, cor = cor$estimate, acc = acc.total, acc.t1 = acc[1], acc.t2 = acc[2], acc.t3 = acc[3])
    } else {
        results = c(RMSE, bias, RVar, cor$p.value, cor$estimate, acc.total, acc[1], acc[2], acc[3])
    }
    
    return (results)    
}




#' @title Linear model obtention and validation for a particular cluster
#' @description Obtaining the linear model and evaluating the results
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists, and shape of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between burned area during fire season months in each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different dataframe
#' @param biome Number between 1 and 13 representing a biome
#' @param cluster Number of the cluster whose burned area time series are we going to obtain
#' @param lt ...
#' @param duration Number of months considered each year for the climate time series calculation
#' @param pvalue Threshold for considering significant each correlation. Default to 0.05
#' @param all Should all the lower possible values of lt and duration be considered?
#' @return list containing the linear model previously calculated, the results of the validation of the model and a index indicating how many correlated indices with this cluster are (ind = -1 is there are not)
lm.obtention <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, lt, duration = 1, pvalue = 0.05, all = F){
    # Obtaining the time series
    if (all == F){
        df = clus.data.preparation(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, lt, duration, pvalue)
    } else {
        for (i in 1:lt){
            for (j in 1:duration){                
                if (i == 1 & j == 1){
                    df = clus.data.preparation(fireSeasons, ba.series, dates, corr.df[[1]][[1]], list.cpcs, biome, cluster, 1, 1, pvalue)
                } else {
                    df.new = clus.data.preparation(fireSeasons, ba.series, dates, corr.df[[i]][[j]], list.cpcs, biome, cluster, i, j, pvalue)                    
                    if (dim(df.new)[2] > 1){
                        if (dim(df.new)[2] == 2){
                            df = cbind(df, df.new)
                            df = df[,-(dim(df)[2]-1)]
                        } else {
                            df = cbind(df, df.new[,-1])
                        }
                    }
                }
            }
        }

        # Second selection of predictors: avoid having a rank-defficient matrix
        while (rankMatrix(as.matrix(df)) < dim(df)[2]){
            ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == cluster)[1]
            k.max = 2
            P = as.integer(substr(colnames(df)[k.max], 2, 2))
            LT = as.integer(substr(colnames(df)[k.max], 5, 5))
            D = as.integer(substr(colnames(df)[k.max], 7, 7))
            pVALUE = corr.df[[LT]][[D]][ind.coords, 2*D+2]
            for (k in 3:dim(df)[2]){
                P = as.integer(substr(colnames(df)[k], 2, 2))
                LT = as.integer(substr(colnames(df)[k], 5, 5))
                D = as.integer(substr(colnames(df)[k], 7, 7))
                if (corr.df[[LT]][[D]][ind.coords, 2*D+2] > pVALUE){
                    pVALUE = corr.df[[LT]][[D]][ind.coords, 2*D+2]
                    k.max = k
                }
            }
            df = df[,-k.max]
        }
    }
    
    ind = dim(df)[2] - 1
    if (dim(df)[2] == 1){
        ind = -1
        return (list(mod = NA, results = NA, ind = ind))
    }
    
    # Training the linear model
    ctrl <- trainControl(method = "LOOCV")# Leave-one-out CV
    mod <- train(ba.serie ~ ., data = df, method = "lm", trControl = ctrl)
    
    # Validating the model
    results = model.validation(mod$pred$obs, mod$pred$pred)
    
    return (list(mod = mod, results = results, ind = ind))
}




#' @title Linear model obtention and validation for all the clusters
#' @description Obtaining the linear models and evaluating the results. Clusters with no correlated indices will not have a model
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists, and shape of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param corr.df Dataframe with the same shape as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indices in the period that we are studying. Data of each index must be in a different dataframe
#' @param lt ...
#' @param duration Number of months considered each year for the climate time series calculation
#' @param mode Type of fire season whose correlation are we going to calculate. It could be 'unimodal' (default), 'bimodal1' (main fire season of bimodal fire seasons) or 'bimodal2' (secondary fire season of bimodal fire seasons)
#' @param pvalue Threshold for considering significant each correlation. Default to 0.05
#' @param all Should all the lower possible values of lt and duration be considered?
#' @return list containing the linear models previously calculated and the results of the validation of each model
lm.all <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, lt, duration = 1, mode = 'unimodal', pvalue = 0.05, all = F){
    
    if (mode == 'unimodal'){
        shape = 1
    } else if (mode == 'bimodal1'){
        shape = 2
    } else if (mode == 'bimodal2'){
        shape = 2
        fireSeasons[which(fireSeasons$form == shape),]$start.1 = fireSeasons[which(fireSeasons$form == shape),]$start.2
        fireSeasons[which(fireSeasons$form == shape),]$end.1 = fireSeasons[which(fireSeasons$form == shape),]$end.2
    }
    
    lm = list()# To store the models
    results = data.frame(biome = 0, cluster = 0, lm.Npred = 0, lm.RMSE = 0, lm.bias = 0, lm.RVar = 0, lm.cor.pvalue = 0, lm.cor = 0, lm.acc = 0, lm.acc.t1 = 0, lm.acc.t2 = 0, lm.acc.t3 = 0)# To store the quality
    
    for (biome in 1:13){
        lm.biome = list()
        clusters = sort(unique(fireSeasons[which(fireSeasons$BIOME == biome),]$cl))
        for (cl in 1:length(clusters)){
                
            cluster = clusters[cl]         
            ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == cluster)
            
            # If the cluster has a different shape, we move to the next cluster
            if (fireSeasons[ind.coords,]$form[1] != shape){
                lm.biome[[cl]] = NA
                add = c(biome, cl, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                results = rbind(results, add)
                next
            }
            
            # Obtains the linear model
            r = lm.obtention(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, lt, duration, pvalue, all)
            
            if (r$ind == -1){
                lm.biome[[cl]] = NA
                add = c(biome, cl, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                results = rbind(results, add)
                next
            }
            
            lm.biome[[cl]] <- r$mod

            # Obtaining the quality
            add = r$results
            add = c(biome, cl, r$ind, add)
            
            results = rbind(results, add)
                                    
        }
        lm[[biome]] = lm.biome
    }
    return (list(lm = lm, results = results[-1,]))
}
