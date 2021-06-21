# Functions used for obtaining the models for the prediction of the burned area in each cluster. Only linear models and random forests are considered

# Packages required
library(tidyr)
library(dplyr)
require(sp)
require(magrittr)
library('caret')
library('randomForest')


# Data required
load('fireSeasonPer75_def.Rdata', verbose = T)
load("data/ba_mon_time_series_masked.Rdata", verbose = T)
load("data/ba_mon_clim_masked_df.Rdata", verbose = T)
load('corrDfAnnualClus.Rdata', verbose = T)


masked_ba_series.log = log1p(masked_ba_series)



#' @title Linear models obtaining
#' @description Obtains a linear model for the prediction of the sum of the burned area during the fire season's months using the average of some climate indexes in these months. A model for each cluster is calculated. Only climate indexes with significant correlation with the burned area in the cluster are used as predictors. Clusters with no correlated indexes or with a different fire seasons's form from the one specified in "mode" argument will not have a model.
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different data frame
#' @param mode Type of fire season that we consider. It could be 'unimodal' (by default), 'bimodal1' (main fire season of bimodal fire seasons) or 'bimodal2' (secondary fire season of bimodal fire seasons). Clusters with different type of fire season will not have a model.
#' @return A list containing another list with all the linear models and a dataframe with information about the quality of each model.
lm.clus <- function(fireSeasons, corr.df, list.cpcs, mode = 'unimodal'){
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
    results = data.frame(biome = 0, cluster = 0, Npred = 0, RMSE = 0, R2 = 0, MAE = 0, RVar = 0, Rp90 = 0)# To store the quality
    ctrl <- trainControl(method = "LOOCV")# leave-one-out CV

    for (biome in 1:13){
        lm.biome = list()
        clusters = sort(unique(fireSeasons[which(fireSeasons$BIOME == biome),]$cl))
        n.clusters = length(clusters)
        for (cl in 1:n.clusters){
                
            clus = clusters[cl]         
            ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == clus)
            # If the cluster has a different form, we move to the next cluster
            if (fireSeasons[ind.coords,]$form[1] != form){
                lm.biome[[cl]] = NA
                add = c(biome, cl, NA, NA, NA, NA, NA, NA)
                results = rbind(results, add)
                next
            }
            
            # Discard the indexes with no significant correlation
            cpcs = c()
            for (i in 3:dim(corr.df)[2]){
                if (!is.na(corr.df[ind.coords, i][1])){
                    cpcs = c(cpcs, i)
                }
            }
            # If there are no indexes with significant correlation, we move to the next cluster
            if (length(cpcs) == 0){
                lm.biome[[cl]] = NA
                add = c(biome, cl, 0, NA, NA, NA, NA, NA)
                results = rbind(results, add)
                next
            }            
            
             # If the fire season takes place in only one year 
            if (fireSeasons[ind.coords,]$start.1[1] <= fireSeasons[ind.coords,]$end.1[1]){
                meses = seq(fireSeasons[ind.coords,]$start.1[1], fireSeasons[ind.coords,]$end.1[1])
                ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)
                   
                ind.meses = ind.meses[1:(length(meses) * floor(length(ind.meses)/length(meses)))]
                
                m = masked_ba_series.log[ind.meses, ind.coords]
                
                # Obtaining burned area time series as the sum of the burned area during the fire season of all the points in the cluster
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
                
                df = data.frame('ba' = ba.serie)
                ind.years = as.numeric(substr(dates[ind.meses], 1, 4))
                
                # Obtaining the climate indexes time series
                if (length(meses) == 1){
                    for (i in 1:length(cpcs)){
                        cpc.serie = list.cpcs[[cpcs[i]-2]][which(list.cpcs[[cpcs[i]-2]]$Year %in% ind.years), c(meses + 1)]
                        df = cbind(df, cpc.serie)
                    }                    
                } else {
                    for (i in 1:length(cpcs)){                    
                        n = list.cpcs[[cpcs[i]-2]][which(list.cpcs[[cpcs[i]-2]]$Year %in% ind.years), c(meses + 1)]
                        cpc.serie = apply(n, 1, mean)
                        df = cbind(df, cpc.serie)
                    }
                }
                colnames(df) = c('ba', cpcs-2)
                
                # Training the linear model
                mod <- train(ba ~ ., data = df, method = "lm", trControl = ctrl)
                lm.biome[[cl]] <- mod
                
                add = c(biome, cl, length(cpcs), mod$results$RMSE, mod$results$Rsquared, mod$results$MAE,
                       var(mod$pred$pred)/var(mod$pred$obs), quantile(mod$pred$pred, prob = 0.9)/quantile(mod$pred$obs,
                        prob = 0.9))
                results = rbind(results, add)
            
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
                
                m = masked_ba_series.log[ind.meses, ind.coords]
                
                # Obtaining burned area time series as the sum of the burned area during the fire season of all the points in the cluster
                ba.serie = c()
                for (i in 1:(length(ind.meses)/length(meses))){
                    ba.serie = c(ba.serie, sum(m[((i-1)*(length(meses)) + 1):(i*length(meses)),1]))
                }
                for (j in 2:length(ind.coords)){
                    l = c()
                    for (i in 1:(length(ind.meses)/length(meses))){
                        l = c(l, sum(m[((i-1)*(length(meses)) + 1):(i*length(meses)),j]))
                    }
                    ba.serie = ba.serie + l
                }
                
                df = data.frame('ba' = ba.serie)
                ind.years = as.numeric(substr(dates[ind.meses], 1, 4))
                
                # Obtaining the climate indexes time series
                for (i in 1:length(cpcs)){
                    n = list.cpcs[[cpcs[i]-2]][which(list.cpcs[[cpcs[i]-2]]$Year %in% ind.years), c(meses + 1)]
                    n = as.vector(t(n))

                    cpc.serie = c()
                    for (j in 1:(length(ind.meses)/length(meses))){
                        cpc.serie = c(cpc.serie, mean(n[(fireSeasons[ind.coords,]$end.1[1] + 1 + (j-1)*length(meses)):
                                        (fireSeasons[ind.coords,]$end.1[1] + j*length(meses))]))
                    }
                    df = cbind(df, cpc.serie)
                }
                
                colnames(df) = c('ba', cpcs-2)
                
                # Training the linear model
                mod <- train(ba ~ ., data = df, method = "lm", trControl = ctrl)
                lm.biome[[cl]] <- mod
                
                add = c(biome, cl, length(cpcs), mod$results$RMSE, mod$results$Rsquared, mod$results$MAE,
                       var(mod$pred$pred)/var(mod$pred$obs), quantile(mod$pred$pred, prob = 0.9)/quantile(mod$pred$obs,
                        prob = 0.9))
                results = rbind(results, add)
            }                        
        }
        lm[[biome]] = lm.biome
    }
    return (list(lm = lm, results = results[-1,]))
}




#' @title Random forest models obtaining
#' @description Obtains a random forest model for the prediction of the sum of the burned area during the fire season's months using the average of some climate indexes in these months. A model for each cluster is calculated. Clusters with a different fire seasons's form from the one specified in "mode" argument will not have a model.
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different data frame
#' @param mode Type of fire season that we consider. It could be 'unimodal' (by default), 'bimodal1' (main fire season of bimodal fire seasons) or 'bimodal2' (secondary fire season of bimodal fire seasons). Clusters with different type of fire season will not have a model.
#' @return A list containing another list with all the random forest models and a dataframe with information about the quality of each model.
rf.clus <- function(fireSeasons, list.cpcs, mode = 'unimodal'){
    if (mode == 'unimodal'){
        form = 1
    } else if (mode == 'bimodal1'){
        form = 2
    } else if (mode == 'bimodal2'){
        form = 2
        fireSeasons[which(fireSeasons$form == form),]$start.1 = fireSeasons[which(fireSeasons$form == form),]$start.2
        fireSeasons[which(fireSeasons$form == form),]$end.1 = fireSeasons[which(fireSeasons$form == form),]$end.2
    }
    set.seed(23)
    rf = list()# To store the model's information
    results = data.frame(biome = 0, cluster = 0, mtry = 0, ntree = 0, RMSE = 0, R2 = 0, MAE = 0, RVar = 0, Rp90 = 0)# To store the quality information
    ctrl <- trainControl(method = "LOOCV")# Leave-one-out CV
    tunegrid <- expand.grid(.mtry=c(1:length(list.cpcs)))# Grid for optimizing the number of predictors
    
    for (biome in 1:13){
        rf.biome = list()
        clusters = sort(unique(fireSeasons[which(fireSeasons$BIOME == biome),]$cl))
        n.clusters = length(clusters)
        for (cl in 1:n.clusters){
                
            clus = clusters[cl]         
            ind.coords = which(fireSeasons$BIOME == biome & fireSeasons$cl == clus)
            # Clusters with different type of fire season are not considered
            if (fireSeasons[ind.coords,]$form[1] != form){
                rf.biome[[cl]] = NA
                add = c(biome, cl, NA, NA, NA, NA, NA, NA, NA)
                results = rbind(results, add)
                next
            }                      
            
            # If the fire season takes place in only one year
            if (fireSeasons[ind.coords,]$start.1[1] <= fireSeasons[ind.coords,]$end.1[1]){
                meses = seq(fireSeasons[ind.coords,]$start.1[1], fireSeasons[ind.coords,]$end.1[1])
                ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)
                   
                ind.meses = ind.meses[1:(length(meses) * floor(length(ind.meses)/length(meses)))]
                
                m = masked_ba_series.log[ind.meses, ind.coords]
                
                # Obtaining burned area time series as the sum of the burned area during the fire season of all the points in the cluster
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
                
                df = data.frame('ba' = ba.serie)
                ind.years = as.numeric(substr(dates[ind.meses], 1, 4))
                
                # Obtaining climate indexes time series
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
                colnames(df) = c('ba', 1:length(list.cpcs))
                
                # Training the first model
                model.def = fit <- train(ba~.,
                   data = df,
                   method = 'rf',
                   tuneGrid = tunegrid,
                   trControl = ctrl,
                   ntree = 250)

                # Training new models with different ntree parameters and storing the best one
                for (ntree in c(500,1000,1500,2000,2500,3000)){
                  fit <- train(ba~.,
                               data = df,
                               method = 'rf',
                               tuneGrid = tunegrid,
                               trControl = ctrl,
                               ntree = ntree)
                  if (fit$results$RMSE[as.numeric(fit$bestTune)] < model.def$results$RMSE[as.numeric(model.def$bestTune)]){
                      model.def = fit
                  }
                }
                ind = as.numeric(model.def$bestTune)
                rv = var(model.def$pred[which(model.def$pred$mtry == ind),]$pred) / var(model.def$pred[which(model.def$pred$mtry == ind),]$obs)# var ratio
                rp90 = quantile(model.def$pred[which(model.def$pred$mtry == ind),]$pred, prob = 0.9) / quantile(model.def$pred[which(model.def$pred$mtry == ind),]$obs, prob = 0.9)# p90 ratio
                add = c(biome, cl, model.def$results$mtry[ind], model.def$finalModel$ntree, model.def$results$RMSE[ind],
                        model.def$results$Rsquared[ind], model.def$results$MAE[ind], rv, rp90)
                results = rbind(results, add)
                rf.biome[[cl]] = model.def
            
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
                
                m = masked_ba_series.log[ind.meses, ind.coords]
                
                # Obtaining burned area time series as the sum of the burned area during the fire season of all the points in the cluster
                ba.serie = c()
                for (i in 1:(length(ind.meses)/length(meses))){
                    ba.serie = c(ba.serie, sum(m[((i-1)*(length(meses)) + 1):(i*length(meses)),1]))
                }
                for (j in 2:length(ind.coords)){
                    l = c()
                    for (i in 1:(length(ind.meses)/length(meses))){
                        l = c(l, sum(m[((i-1)*(length(meses)) + 1):(i*length(meses)),j]))
                    }
                    ba.serie = ba.serie + l
                }
                
                df = data.frame('ba' = ba.serie)
                ind.years = as.numeric(substr(dates[ind.meses], 1, 4))
                
                # Obtaining climate indexes time series
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
                
                colnames(df) = c('ba', 1:length(list.cpcs))
                
                # Training the first model
                model.def = fit <- train(ba~.,
                   data = df,
                   method = 'rf',
                   tuneGrid = tunegrid,
                   trControl = ctrl,
                   ntree = 250)

                # Training new models with different ntree parameters and storing the best one
                for (ntree in c(500,1000,1500,2000,2500,3000)){

                  fit <- train(ba~.,
                               data = df,
                               method = 'rf',
                               tuneGrid = tunegrid,
                               trControl = ctrl,
                               ntree = ntree)
                  if (fit$results$RMSE[as.numeric(fit$bestTune)] < model.def$results$RMSE[as.numeric(model.def$bestTune)]){
                      model.def = fit
                  }
                }
                ind = as.numeric(model.def$bestTune)
                rp90 = quantile(model.def$pred[which(model.def$pred$mtry == ind),]$pred, prob = 0.9) / quantile(model.def$pred[which(model.def$pred$mtry == ind),]$obs, prob = 0.9)# p90 ratio
                rv = var(model.def$pred[which(model.def$pred$mtry == ind),]$pred) / var(model.def$pred[which(model.def$pred$mtry == ind),]$obs)# var ratio
                add = c(biome, cl, model.def$results$mtry[ind], model.def$finalModel$ntree, model.def$results$RMSE[ind],
                        model.def$results$Rsquared[ind], model.def$results$MAE[ind], rv, rp90)
                results = rbind(results, add)
                rf.biome[[cl]] = model.def
            }                        
        }
        rf[[biome]] = rf.biome
    }
    return (list(rf = rf, results = results[-1,]))
}