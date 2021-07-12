# Functions to build and test the models


# Packages required
library(tidyr)
library(dplyr)
require(sp)
require(magrittr)
library(caret)
library(randomForest)
library(tree)


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
#' @param t persistence index. Default to 0
#' @return dataframe whose new columns will be the climate indexes time series
get.cpc.serie <- function(fireSeasons, dates, ind.coords, list.cpcs, meses, ind.meses, df, t = 0){
    
    if (t > 0){# persistence
        ind.meses = ind.meses - t# update the indexes
        meses = (meses - t) %% 12# update the months whose climate index time are we going to obtain
        meses[which(meses == 0)] = 12
        if (ind.meses[1] < 1){# if ind.meses has non-positive values, we need data for the year which is before the first one of dates array 
            ind.years = as.numeric(substr(dates[ind.meses[which(ind.meses > 0)]], 1, 4))
            ind.years = c(min(ind.years)-1, ind.years)
        } else {
            ind.years = as.numeric(substr(dates[ind.meses], 1, 4))
        }
        
        # Recalculate the months of beginning and end of the fire season
        if (1 %in% meses & 12 %in% meses){
            dif = diff(sort(meses))
            k = 1
            while (dif[k] == 1 & k < length(dif)){
                k = k + 1
            }
            fireSeasons[ind.coords,]$start.1[1] = sort(meses)[k+1]
            fireSeasons[ind.coords,]$end.1[1] = sort(meses)[k]
        } else {
            fireSeasons[ind.coords,]$start.1[1] = min(meses)
            fireSeasons[ind.coords,]$end.1[1] = max(meses)
        }
    } else {        
        ind.years = as.numeric(substr(dates[ind.meses], 1, 4))
    }
                                                    
    # If the fire season takes place in only one year 
    if (fireSeasons[ind.coords,]$start.1[1] <= fireSeasons[ind.coords,]$end.1[1]){

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
#' @param t persistence index. Default to 0
#' @return list containing in the element 'df' a dataframe whose first column will be the burned area time serie and the other ones will be the climate indexes time series. The second element is a vector containing the climate indexes with significant correlation
clus.data.preparation <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue = 0.05, t = 0, useDeltas = F){
    
    # Positions of the cluster data points in the dataframes
    ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == cluster)

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

    cpcs.all = c()
    j = 0
    while (j <= t){# One iteration per value between 0 and the persistence t        
        # Indexes with significant correlation with this cluster
        cpcs = c()
        for (i in (j*length(list.cpcs)+1):((j+1)*length(list.cpcs))){
            if (!is.na(corr.df[ind.coords, 2*i+2][1]) && corr.df[ind.coords, 2*i+2][1] < pvalue){
                cpcs = c(cpcs, i)
            }
        }
        # Obtaining the climate indexes time series
        df <- get.cpc.serie(fireSeasons, dates, ind.coords, list.cpcs, meses, ind.meses, df, t = j)        
        
        cpcs.all = c(cpcs.all, cpcs)
        j = j + 1
    }
    
    if (is.null(names(list.cpcs))){
        colnames(df) = c('ba', (1:((t+1)*length(list.cpcs))))# con esto funcionaba
    } else {
        if (t == 0){
            colnames(df) = c('ba', names(list.cpcs))
        } else {
            names = c('ba')
            for (i in 0:t){
                for (j in 1:length(list.cpcs)){
                    names = c(names, paste(names(list.cpcs)[j], '.', toString(i), sep = ''))
                }
            }
            colnames(df) = names
        }
    }
    
    if (useDeltas == T){
        df = apply(df, 2, diff)
    }
            
    return (list(df = df, cpcs = cpcs.all))
}


#' @title Model's validation
#' @description Obtains the value of some metrics for the results obtained by a model
#' @param obs Vector with real values
#' @param pred Vector with predicted values
#' @returns dataframe containing the value of the RMSE, bias, correlation, variance ratio, total accuracy and tercile accuracy
model.validation <- function(obs, pred, as.df = F){
    
    rss = sum((pred - obs) ^ 2)  ## residual sum of squares
    RMSE = abs(100*sqrt(rss/length(pred)) / mean(obs))
    bias = 100 * (mean(pred) / mean(obs)) / mean(obs)
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




#' @title Linear model obtention and validation
#' @description Obtaining the linear model and evaluating the results. It also plots the real and the predicted values of the burned area, the burned area versus the climate index with significant correlation if there is only one and climate index time serie
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different dataframe
#' @param biome Number between 1 and 13 representing a biome
#' @param cluster Number of the cluster whose burned area time series are we going to obtain
#' @param pvalue Threshold for considering significant each correlation
#' @param t persistence index. Default to 0
#' @return the linear model previously calculated
lm.clus.plot <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue = 0.05, t = 0, useDeltas = F){
    # Obtaining the time series
    data = clus.data.preparation(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue, t, useDeltas)
    df = data$df
    cpcs = data$cpcs
    
    # Training the linear model
    ctrl <- trainControl(method = "LOOCV")# Leave-one-out CV
    mod <- train(ba ~ ., data = df[,c(1, c(cpcs + 1))], method = "lm", trControl = ctrl)
    
    # Validating the model
    results = model.validation(mod$pred$obs, mod$pred$pred, as.df = T)
    
    # Doing the plots
    par(mfrow=c(2,2))
    
    # Observed vs predicted plot
    plot(mod$pred$obs, col = 'green', main = 'Observed vs predicted', ylab = 'Burned area', type = 'l')
    lines(mod$pred$pred, col = 'red')
    
    # Climate time series plots
    #for (i in 1:length(cpcs)){
    #    plot(df[,cpcs[i]+1], main = paste('Climate index', toString(cpcs[i])))
    #}
    
    # Climate index vs burned area
    #if (length(cpcs) == 1){
    #    ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == cluster)
    #    plot(df[,cpcs[i]+1], df$ba, main = paste('Climate index vs ba with corr = ', 
    #                                      toString(round(corr.df[ind.coords, 2*cpcs + 1][1], digits = 3))), 
    #         ylab = 'Burned area')
    #}
    
    # Printing the validation results
    print(results)
    
    return (mod)
}



#' @title Linear model obtention and validation for a particular cluster
#' @description Obtaining the linear model and evaluating the results
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different dataframe
#' @param biome Number between 1 and 13 representing a biome
#' @param cluster Number of the cluster whose burned area time series are we going to obtain
#' @param pvalue Threshold for considering significant each correlation
#' @param t persistence index. Default to 0
#' @return list containing the linear model previously calculated, the results of the validation of the model and a index indicating how many correlated indexes with this cluster are (ind = -1 is there are not)
lm.obtention <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue = 0.05, t = 0, useDeltas = F){
    # Obtaining the time series
    data = clus.data.preparation(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue, t, useDeltas)
    df = data$df
    cpcs = data$cpcs
    
    if (length(cpcs) == 0){
        ind = -1
        return (list(mod = NA, results = NA, ind = ind))
    } else {
        ind = length(cpcs)
    }
    
    # Training the linear model
    ctrl <- trainControl(method = "LOOCV")# Leave-one-out CV
    mod <- train(ba ~ ., data = df[,c(1, c(cpcs + 1))], method = "lm", trControl = ctrl)
    
    # Validating the model
    results = model.validation(mod$pred$obs, mod$pred$pred)
    
    return (list(mod = mod, results = results, ind = ind))
}



#' @title Linear model obtention and validation for all the clusters
#' @description Obtaining the linear models and evaluating the results. Clusters with no correlated indexes will not have a model
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different dataframe
#' @param mode Type of fire season whose correlation are we going to calculate. It could be 'unimodal' (by default), 'bimodal1' (main fire season of bimodal fire seasons) or 'bimodal2' (secondary fire season of bimodal fire seasons)
#' @param pvalue Threshold for considering significant each correlation
#' @param t persistence index. Default to 0
#' @return list containing the linear models previously calculated and the results of the validation of each model
lm.all <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, mode = 'unimodal', pvalue = 0.05, t = 0, useDeltas = F){
    
    if (mode == 'unimodal'){
        form = 1
    } else if (mode == 'bimodal1'){
        form = 2
    } else if (mode == 'bimodal2'){
        form = 2
        fireSeasons[which(fireSeasons$form == form),]$start.1 = fireSeasons[which(fireSeasons$form == form),]$start.2
        fireSeasons[which(fireSeasons$form == form),]$end.1 = fireSeasons[which(fireSeasons$form == form),]$end.2
    }
    
    lm = list()# To store the models
    results = data.frame(biome = 0, cluster = 0, lm.Npred = 0, lm.RMSE = 0, lm.bias = 0, lm.RVar = 0, lm.cor.pvalue = 0, lm.cor = 0, lm.acc = 0, lm.acc.t1 = 0, lm.acc.t2 = 0, lm.acc.t3 = 0)# To store the quality
    
    for (biome in 1:13){
        lm.biome = list()
        clusters = sort(unique(fireSeasons[which(fireSeasons$BIOME == biome),]$cl))
        for (cl in 1:length(clusters)){
                
            cluster = clusters[cl]         
            ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == cluster)
            
            # If the cluster has a different form, we move to the next cluster
            if (fireSeasons[ind.coords,]$form[1] != form){
                lm.biome[[cl]] = NA
                add = c(biome, cl, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                results = rbind(results, add)
                next
            }
            
            r = lm.obtention(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue, t, useDeltas)
            
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



#' @title Random forest model obtention and validation for a particular cluster
#' @description Obtaining the random forest model and evaluating the results
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different dataframe
#' @param biome Number between 1 and 13 representing a biome
#' @param cluster Number of the cluster whose burned area time series are we going to obtain
#' @param pvalue Threshold for considering significant each correlation
#' @param t persistence index. Default to 0
#' @return list containing the random forest model previously calculated, the results of the validation of the model and a index indicating how many trees are in the model
rf.obtention.plot <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue = 0.05, t = 0, useDeltas = F){
    # Obtaining the time series
    data = clus.data.preparation(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue, t, useDeltas)
    df = data$df        
    
    
    # Training the random forest model
    ctrl <- trainControl(method = "LOOCV")# Leave-one-out CV
    tunegrid <- expand.grid(.mtry=ceiling(length(list.cpcs)/3))# Grid for optimizing the number of predictors
    set.seed(23)
    
    # Training the first model
    ntree.def = 10
    model.def <- train(ba~.,
       data = df,
       method = 'rf',
       tuneGrid = tunegrid,
       trControl = ctrl,
       ntree = 25)
    
    results = cbind(data.frame(ntree = 25), model.validation(model.def$pred$obs, model.def$pred$pred, as.df = T))

    # Training new models with different ntree parameters and storing the best one
    for (ntree in seq(50, 200, 25)){

        fit <- train(ba~.,
                   data = df,
                   method = 'rf',
                   tuneGrid = tunegrid,
                   trControl = ctrl,
                   ntree = ntree)
        
        add = cbind(data.frame(ntree = ntree), model.validation(fit$pred$obs, fit$pred$pred, as.df = T))                    
        results = rbind(results, add)
    }
    
    par(mfrow=c(2,2))
    plot(results$ntree, results$RMSE, main = 'RMSE')
    #plot(results$ntree, results$bias, main = 'bias')
    plot(results$ntree, results$RVar, main = 'RVar')
    plot(results$ntree, results$cor.pvalue, main = 'cor p-value')
    plot(results$ntree, results$acc, main = 'acc')
    
    return (results)
}




#' @title Random forest model obtention and validation for a particular cluster
#' @description Obtaining the random forest model and evaluating the results
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different dataframe
#' @param biome Number between 1 and 13 representing a biome
#' @param cluster Number of the cluster whose burned area time series are we going to obtain
#' @param pvalue Threshold for considering significant each correlation
#' @param t persistence index. Default to 0
#' @return list containing the random forest model previously calculated, the results of the validation of the model and a index indicating how many trees are in the model
rf.obtention <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue = 0.05, t = 0, useDeltas = F){
    # Obtaining the time series
    data = clus.data.preparation(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue, t, useDeltas)
    df = data$df        
    
    # Training the random forest model
    ctrl <- trainControl(method = "LOOCV")# Leave-one-out CV
    tunegrid <- expand.grid(.mtry=ceiling(length(list.cpcs)/3))# Grid for optimizing the number of predictors
    set.seed(23)
    
    # Training the first model
    ntree.def = 10
    model.def <- train(ba~.,
       data = df,
       method = 'rf',
       tuneGrid = tunegrid,
       trControl = ctrl,
       ntree = 10)

    # Training new models with different ntree parameters and storing the best one
    for (ntree in seq(20, 200, 20)){

      fit <- train(ba~.,
                   data = df,
                   method = 'rf',
                   tuneGrid = tunegrid,
                   trControl = ctrl,
                   ntree = ntree)
      if (fit$results$RMSE < model.def$results$RMSE){
          model.def = fit
          ntree.def = ntree
      }
    }
    
    # Validating the model
    results = model.validation(model.def$pred$obs, model.def$pred$pred)
    
    return (list(mod = model.def, results = results, ntree = ntree.def))
}



#' @title Random forest models obtention and validation for all the clusters
#' @description Obtaining the random forest models and evaluating the results
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different dataframe
#' @param mode Type of fire season whose correlation are we going to calculate. It could be 'unimodal' (by default), 'bimodal1' (main fire season of bimodal fire seasons) or 'bimodal2' (secondary fire season of bimodal fire seasons)
#' @param pvalue Threshold for considering significant each correlation
#' @param t persistence index. Default to 0
#' @return list containing the random forests models previously calculated and the results of the validation of each model
rf.all <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, mode = 'unimodal', pvalue = 0.05, t = 0, useDeltas = F){
    
    if (mode == 'unimodal'){
        form = 1
    } else if (mode == 'bimodal1'){
        form = 2
    } else if (mode == 'bimodal2'){
        form = 2
        fireSeasons[which(fireSeasons$form == form),]$start.1 = fireSeasons[which(fireSeasons$form == form),]$start.2
        fireSeasons[which(fireSeasons$form == form),]$end.1 = fireSeasons[which(fireSeasons$form == form),]$end.2
    }
    
    rf = list()# To store the models
    results = data.frame(biome = 0, cluster = 0, rf.Ntree = 0, rf.RMSE = 0, rf.bias = 0, rf.RVar = 0, rf.cor.pvalue = 0, rf.cor = 0, rf.acc = 0, rf.acc.t1 = 0, rf.acc.t2 = 0, rf.acc.t3 = 0)# To store the quality
    
    for (biome in 1:13){
        rf.biome = list()
        clusters = sort(unique(fireSeasons[which(fireSeasons$BIOME == biome),]$cl))
        for (cl in 1:length(clusters)){
                
            cluster = clusters[cl]         
            ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == cluster)
            
            # If the cluster has a different form, we move to the next cluster
            if (fireSeasons[ind.coords,]$form[1] != form){
                rf.biome[[cl]] = NA
                add = c(biome, cl, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                results = rbind(results, add)
                next
            }
            
            r = rf.obtention(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue, t, useDeltas)
            
            rf.biome[[cl]] <- r$mod

            # Obtaining the quality
            add = r$results
            add = c(biome, cl, r$ntree, add)
            
            results = rbind(results, add)
                                    
        }
        rf[[biome]] = rf.biome
    }
    return (list(rf = rf, results = results[-1,]))
}





#' @title Knn model obtention and validation for a particular cluster
#' @description Obtaining the knn model and evaluating the results
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different dataframe
#' @param biome Number between 1 and 13 representing a biome
#' @param cluster Number of the cluster whose burned area time series are we going to obtain
#' @param pvalue Threshold for considering significant each correlation
#' @param t persistence index. Default to 0
#' @return list containing the knn model previously calculated and the results of the validation of the model
knn.obtention <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue = 0.05, t = 0, useDeltas = F){
    # Obtaining the time series
    data = clus.data.preparation(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue, t, useDeltas)
    df = data$df        
    
    # Training the random forest model
    ctrl <- trainControl(method = "LOOCV")# Leave-one-out CV
    tunegrid <- expand.grid(k = 1)
    
    set.seed(23)
    # Training the model
    model.def <- train(ba~.,
       data = df,
       method = 'knn',
       tuneGrid = tunegrid,
       trControl = ctrl,
       preProcess = c("center","scale"))
    
    # Validating the model
    results = model.validation(model.def$pred$obs, model.def$pred$pred)
    
    return (list(mod = model.def, results = results))
}



#' @title Knn models obtention and validation for all the clusters
#' @description Obtaining the knn models and evaluating the results
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different dataframe
#' @param mode Type of fire season whose correlation are we going to calculate. It could be 'unimodal' (by default), 'bimodal1' (main fire season of bimodal fire seasons) or 'bimodal2' (secondary fire season of bimodal fire seasons)
#' @param pvalue Threshold for considering significant each correlation
#' @param t persistence index. Default to 0
#' @return list containing the knn models previously calculated and the results of the validation of each model
knn.all <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, mode = 'unimodal', pvalue = 0.05, t = 0, useDeltas = F){
    
    if (mode == 'unimodal'){
        form = 1
    } else if (mode == 'bimodal1'){
        form = 2
    } else if (mode == 'bimodal2'){
        form = 2
        fireSeasons[which(fireSeasons$form == form),]$start.1 = fireSeasons[which(fireSeasons$form == form),]$start.2
        fireSeasons[which(fireSeasons$form == form),]$end.1 = fireSeasons[which(fireSeasons$form == form),]$end.2
    }
    
    knn = list()# To store the models
    results = data.frame(biome = 0, cluster = 0, knn.RMSE = 0, knn.bias = 0, knn.RVar = 0, knn.cor.pvalue = 0, knn.cor = 0, knn.acc = 0, knn.acc.t1 = 0, knn.acc.t2 = 0, knn.acc.t3 = 0)# To store the quality
    
    for (biome in 1:13){
        knn.biome = list()
        clusters = sort(unique(fireSeasons[which(fireSeasons$BIOME == biome),]$cl))
        for (cl in 1:length(clusters)){
                
            cluster = clusters[cl]         
            ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == cluster)
            
            # If the cluster has a different form, we move to the next cluster
            if (fireSeasons[ind.coords,]$form[1] != form){
                knn.biome[[cl]] = NA
                add = c(biome, cl, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                results = rbind(results, add)
                next
            }
            
            r = knn.obtention(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue, t, useDeltas)
            
            knn.biome[[cl]] <- r$mod

            # Obtaining the quality
            add = r$results
            add = c(biome, cl, add)
            
            results = rbind(results, add)
                                    
        }
        knn[[biome]] = knn.biome
    }
    return (list(knn = knn, results = results[-1,]))
}





#' @title Regression tree model obtention and validation for a particular cluster
#' @description Obtaining the regression tree model and evaluating the results. 
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param dates Array containing the dates of each burned area observation
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different dataframe
#' @param biome Number between 1 and 13 representing a biome
#' @param cluster Number of the cluster whose burned area time series are we going to obtain
#' @param pvalue Threshold for considering significant each correlation
#' @param t persistence index. Default to 0
#' @return list containing the tree previously calculated and the dataframe used for obtaining and validating the model
tree.obtention <- function(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue = 0.05, t = 0){
    # Obtaining the time series
    data = clus.data.preparation(fireSeasons, ba.series, dates, corr.df, list.cpcs, biome, cluster, pvalue, t)
    df = data$df        
    
    set.seed(23)
    
    # Initial tree
    i0 = 1
    t0 = -1
    j0 = -1
    tree.old <- tree(ba ~ ., data = df[-i0,])
    
    # Finding the best values for the parameters
    for (i in 1:length(df[,1])){# leave one out
        df0 = df[-i,]
        for (j in c(2,3,5,7)){# minsize parameter
            for (t in c(0.01, 0.05, 0.075, 0.25, 0.5)){# mindev parameter
                control.pars <- tree.control(nobs = nrow(df0), mindev = t, minsize = j)
                tree.new <- tree(ba ~ ., data = df0, control = control.pars)
                if ((predict(tree.new, newdata = df[i,])- df[i,1])^2 < (predict(tree.old, newdata = df[i0,])- df[i0,1])^2){
                    tree.old <- tree.new
                    i0 = i
                    j0 = j
                    t0 = t
                }
            }
        }
    }
    plot(tree.old)
    text(tree.old)
    
    results = model.validation(df[,1], predict(tree.old, df[,-1]), as.df = T)
    
    print(results)
    cat('Minsize =', j0, 'Mindev =', t0)
    
    return (list('tree'=tree.old, 'df'=df))
}