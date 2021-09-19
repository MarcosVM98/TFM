# Functions used for obtaining the correlation between burned area during the fire season and climate indexes.

# Packages required
library(tidyr)
library(dplyr)
require(sp)
require(magrittr)
require(RColorBrewer)

# Color palette for the plots
#group.colors <- colorRampPalette(c(brewer.pal(8, "Dark2"), brewer.pal(8, "Accent")))
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
#' @description Obtaining the annual climate index time series of the cpcs given. The value of each year is the average of the value of the index during the fire seasons months
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



#' @title Annual correlation per cluster calculation
#' @description Obtains the correlation between the sum of the burned area of the fire season's months of each pixel of the cluster and the average of the indexes in these months for each cluster. Stores it in a dataframe and plots the results.
#' @param ba.series Dataframe containing the monthly burned area time series for each point
#' @param fireSeasons Data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
#' @param coords Array containing the x and y coordinates of each point of the ba.series dataframe
#' @param dates Array containing the dates of each burned area observation
#' @param cpc Data frame containing the value of the cpc climate index per month. First column must be the year and the other ones have to be the months.
#' @param name Name of the climate index
#' @param corr.df Dataframe with the same form as masked_coords containing the correlation between each cluster and each climate index
#' @param mode Type of fire season whose correlation are we going to calculate. It could be 'unimodal' (by default), 'bimodal1' (main fire season of bimodal fire seasons) or 'bimodal2' (secondary fire season of bimodal fire seasons)
#' @param pvalue Threshold for considering significant each correlation
#' @param t persistence index. Default to 0
#' @return A copy of corr.df dataframe with a new column with the correlation
corr.annual.clus <- function(ba.series, fireSeasons, coords, dates, cpc, name, corr.df, mode = 'unimodal', pvalue = 0.05, t = 0, useDeltas = F){    
    
    if (mode == 'unimodal'){
        form = 1
    } else if (mode == 'bimodal1'){
        form = 2
    } else if (mode == 'bimodal2'){
        form = 2
        fireSeasons[which(fireSeasons$form == form),]$start.1 = fireSeasons[which(fireSeasons$form == form),]$start.2
        fireSeasons[which(fireSeasons$form == form),]$end.1 = fireSeasons[which(fireSeasons$form == form),]$end.2
    }
    
    corr.df$ind.cpc = NA
    corr.df$ind.cpc.pvalue = NA
    
    for (biome in 1:13){
        clusters = sort(unique(fireSeasons[which(fireSeasons$form == form 
                                                     & fireSeasons$BIOME == biome),]$cl))
        n.clusters = length(clusters)
        if (n.clusters < 1){
            next
        }
        for (cl in 1:n.clusters){
            clus = clusters[cl]         
            ind.coords = which(fireSeasons$form == form & fireSeasons$BIOME == biome & fireSeasons$cl == clus)
            
            if (fireSeasons[ind.coords,]$start.1[1] <= fireSeasons[ind.coords,]$end.1[1]){
                meses = seq(fireSeasons[ind.coords,]$start.1[1], fireSeasons[ind.coords,]$end.1[1])
                ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)

            } else {
                meses = c(seq(1, fireSeasons[ind.coords,]$end.1[1]), seq(fireSeasons[ind.coords,]$start.1[1], 12))                
                ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)
                
                if (fireSeasons[ind.coords,]$end.1[1] <= 4){# because the last month with burned area data is April
                    ind.meses = ind.meses[(fireSeasons[ind.coords,]$end.1[1]+1):(length(ind.meses))]
                } else {
                    ind.meses = ind.meses[(fireSeasons[ind.coords,]$end.1[1]+1):(length(ind.meses)-4-(13-fireSeasons[ind.coords,]$start.1[1]))]
                }
            }
            
            ind.meses = ind.meses[1:(length(meses) * floor(length(ind.meses)/length(meses)))]
                
            m = ba.series[ind.meses, ind.coords]                

            # Obtaining burned area time series as the sum of the burned area during the fire season of all the points in the cluster
            ba.serie <- get.ba.serie(m, meses, ind.meses, ind.coords)   

            # Obtaining climate indexes time series
            cpc.serie <- get.cpc.serie(fireSeasons, dates, ind.coords, list(cpc), meses, ind.meses, df = as.data.frame(ba.serie), t)
            
            if (useDeltas == F){
                test = cor.test(ba.serie, cpc.serie[,-1], method = 'pearson')
            } else {
                test = cor.test(diff(ba.serie), diff(cpc.serie[,-1]), method = 'pearson')
            }
            corr.pvalue = test$p.value
            corr = test$estimate
                                    
            corr.df[ind.coords,]$ind.cpc.pvalue = corr.pvalue
            corr.df[ind.coords,]$ind.cpc = corr
        }
    }

    # Plot results
    arg.list <- list(col.regions = group.colors[11:1][-6],
                      at = seq(-1, 1, 0.2), main = paste(name, mode, "correlation per cluster"))
    v <- corr.df$ind.cpc
    v[which(is.na(corr.df$ind.cpc.pvalue))] = NA
    v[which(corr.df$ind.cpc.pvalue >= pvalue)] = NA

    df1 <- cbind.data.frame(coords, v)
    coordinates(df1) <- c(1,2)
    gridded(df1) <- TRUE
    arg.list[["sp.layout"]] <- list("sp.lines", coast.lines)
    arg.list[["obj"]] <- df1
    arg.list[["zcol"]] <- 1
    arg.list[["ylim"]] <- c(-90,90)
    arg.list[["xlim"]] <- c(-180,180)
    do.call("spplot", arg.list) %>% print()

    # Change the column names of the dataframe
    if (mode == 'bimodal2'){
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = paste(name, '.2', sep = '')
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc.pvalue'] = paste(name, '.2.pvalue', sep = '')
    } else if (mode == 'bimodal1'){
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = paste(name, '.1', sep = '')
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc.pvalue'] = paste(name, '.1.pvalue', sep = '')
    } else {
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = name
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc.pvalue'] = paste(name, '.pvalue', sep = '')
    }
    
    return (corr.df)
}
