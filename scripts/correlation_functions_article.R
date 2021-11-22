# Functions used for obtaining the correlation between burned area during the fire season and climate indices.

# Packages required
library(tidyr)
library(dplyr)
require(sp)
require(magrittr)
require(RColorBrewer)

# Color palette for the plots
group.colors <- brewer.pal(11, 'Spectral')

# Data required
load("scripts/worldmap.Rdata", verbose = T)



#' @title Burned area time series calculation
#' @description Obtaining the annual burned area time series as the sum of the burned area during the fire season in each year of all the points in the cluster
#' @param m Subset of the masked_ba_series.log dataframe containing only the burned area data in the points of the cluster during the months of the fire season
#' @param meses Vector containing the months of the fire season of that cluster
#' @param ind.meses Vector containing the positions in the masked_ba_series.log dataframe of the months of the fire seasons during the period 2001-2020
#' @param ind.coords Vector containing the position of the points of the cluster in masked_coords dataframe
#' @return vector with the burned area time series
get.ba.serie.detrended <- function(m, meses, ind.meses, ind.coords){    
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
        ind = 1:length(ba.serie)
        dff = data.frame(ind, ba.serie)
        ba.serie = lm(ba.serie ~ ind, data = dff)$residuals
    }        
    return (ba.serie)
}



#' @title Climate index time series calculation
#' @description Obtaining the annual climate index time series of the cpcs given. The value of each year is the average of the value of the index during the fire seasons months
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param ind.coords Vector containing the position of the points of the cluster in masked_coords dataframe
#' @param list.cpcs List containing the values of the climate indexes in the period that we are studying. Data of each index must be in a different data frame
#' @param ind.start Vector containing the positions in the masked_ba_series.log array of the starting month of the fire season
#' @param df dataframe whose new columns will be the climate indices time series
#' @param lt ...
#' @param duration ...
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



#' @title Annual correlation per cluster calculation
#' @description Obtains the correlation between the sum of the burned area of the fire season's months of each pixel of the cluster and the average of the indexes in these months for each cluster. Stores it in a dataframe and plots the results.
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param coords Array containing the x and y coordinates of each point of the ba.series dataframe
#' @param dates Array containing the dates of each burned area observation
#' @param cpc Data frame containing the value of the cpc climate index per month. First column must be the year and the other ones have to be the months.
#' @param name Name of the climate index
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param lt ...
#' @param duration. ...
#' @param mode Type of fire season whose correlation are we going to calculate. It could be 'unimodal' (by default), 'bimodal1' (main fire season of bimodal fire seasons) or 'bimodal2' (secondary fire season of bimodal fire seasons)
#' @param pvalue Threshold for considering significant each correlation
#' @return A copy of corr.df dataframe with a new column with the correlation
corr.annual.clus <- function(ba.series, fireSeasons, coords, dates, cpc, name, corr.df, lt, duration = 1, mode = 'unimodal', pvalue = 0.05){ 
    
    if (mode == 'unimodal'){
        shape = 1
    } else if (mode == 'bimodal1'){
        shape = 2
    } else if (mode == 'bimodal2'){
        shape = 2
        fireSeasons[which(fireSeasons$form == shape),]$start.1 = fireSeasons[which(fireSeasons$form == shape),]$start.2
        fireSeasons[which(fireSeasons$form == shape),]$end.1 = fireSeasons[which(fireSeasons$form == shape),]$end.2
    }
    
    corr.df$ind.cpc = NA
    corr.df$ind.cpc.pvalue = NA
    
    for (biome in 1:13){
        clusters = sort(unique(fireSeasons[which(fireSeasons$form == shape & fireSeasons$BIOME == biome),]$cl))
        n.clusters = length(clusters)
        if (n.clusters < 1){
            next
        }
        for (cl in 1:n.clusters){
            clus = clusters[cl]         
            ind.coords = which(fireSeasons$form == shape & fireSeasons$BIOME == biome & fireSeasons$cl == clus)
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
           
            # Obtaining climate indexes time series
            cpc.serie <- get.cpc.serie.detrended(fireSeasons, dates, ind.coords, list(cpc), ind.start, df = as.data.frame(ba.serie), lt, duration)
            
            # Correlation
            test = cor.test(ba.serie, cpc.serie[,-1], method = 'pearson')            
            corr.pvalue = test$p.value
            corr = test$estimate
                                    
            corr.df[ind.coords,]$ind.cpc.pvalue = corr.pvalue
            corr.df[ind.coords,]$ind.cpc = corr
        }
    }
   
    # Plot results
    arg.list <- list(col.regions = group.colors[11:1][-6],
                      at = seq(-1, 1, 0.2), main = paste(name, ' LT', lt, ' Dur', duration, sep = ''))
    v <- corr.df$ind.cpc
    v[which(is.na(corr.df$ind.cpc.pvalue))] = NA
    v[which(corr.df$ind.cpc.pvalue >= pvalue)] = NA

    df1 <- cbind.data.frame(coords, v)
    coordinates(df1) <- c(1,2)
    gridded(df1) <- TRUE
    arg.list[["sp.layout"]] <- list("sp.lines", coast.lines)
    arg.list[["obj"]] <- df1
    arg.list[["ylim"]] <- c(-90,90)
    arg.list[["xlim"]] <- c(-180,180)
    do.call("spplot", arg.list) %>% print()
    
    # Change the column names of the dataframe
    if (mode == 'bimodal2'){
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = paste(name, 'lt', lt, 'dur', duration, '.2', sep = '')
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc.pvalue'] = paste(name, 'lt', lt, 'dur', duration, '.2.pvalue', sep = '')
    } else if (mode == 'bimodal1'){
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = paste(name, 'lt', lt, 'dur', duration, '.1', sep = '')
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc.pvalue'] = paste(name, 'lt', lt, 'dur', duration, '.1.pvalue', sep = '')
    } else {
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = paste(name, 'lt', lt, 'dur', duration, sep = '.')
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc.pvalue'] = paste(name, 'lt', lt, 'dur', duration, 'pvalue', sep = '.')
    }    
    
    return (corr.df)
}
