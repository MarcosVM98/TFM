# Functions used for obtaining the correlation between burned area during the fire season and climate indexes.

# Packages required
library(tidyr)
library(dplyr)
require(sp)
require(magrittr)
require(RColorBrewer)

#group.colors <- colorRampPalette(c(brewer.pal(8, "Dark2"), brewer.pal(8, "Accent")))
group.colors <- brewer.pal(11, 'Spectral')

# Data required
load("scripts/worldmap.Rdata", verbose = T)
load("data/ba_mon_time_series_masked.Rdata", verbose = T)
load("fireSeasonPer75_def.Rdata", verbose = T)

masked_ba_series.log = log1p(masked_ba_series)
fireSeasonMedian_def = fireSeasonPer75_def




#' @title Monthly correlation calculation
#' @description Obtains the correlation bbetween burned arrea of each fire season's month and the value of the indexes in this month for each pixel. Stores it in a dataframe and plots the results.
#' @param cpc Data frame containing the value of the cpc climate index per month. First column must be the year and the other ones have to be the months.
#' @param name Name of the climate index
#' @param corr.df.monthly Dataframe with the same form as masked_coords
#' @param threshold Value used as a threshold for the plots. Only pixels which abs(cor) > threshold are plotted
#' @return A copy of corr.df.monthly dataframe with a new column with the correlation
corr.monthly <- function(cpc, name, corr.df.monthly){
    corr.df.monthly$ind.cpc = NA
    
    for (biome in 1:13){
        clusters = sort(unique(fireSeasonMedian_def[which(fireSeasonMedian_def$BIOME == biome),]$cl))
        n.clusters = length(clusters)
        for (cl in 1:n.clusters){
            clus = clusters[cl]
            ind.coords = which(fireSeasonMedian_def$BIOME == biome & fireSeasonMedian_def$cl == clus)
            
            if (fireSeasonMedian_def[ind.coords,]$form[1] == 1){
                if (fireSeasonMedian_def[ind.coords,]$start.1[1] <= fireSeasonMedian_def[ind.coords,]$end.1[1]){
                    meses = seq(fireSeasonMedian_def[ind.coords,]$start.1[1], fireSeasonMedian_def[ind.coords,]$end.1[1])
                } else {
                    meses = c(seq(1, fireSeasonMedian_def[ind.coords,]$end.1[1]), seq(fireSeasonMedian_def[ind.coords,]$start.1[1], 12))
                }
            } else if (fireSeasonMedian_def[ind.coords,]$form[1] == 2){
                if (fireSeasonMedian_def[ind.coords,]$start.1[1] <= fireSeasonMedian_def[ind.coords,]$end.1[1]
                    & fireSeasonMedian_def[ind.coords,]$start.2[1] <= fireSeasonMedian_def[ind.coords,]$end.2[1]){
                    meses = c(seq(fireSeasonMedian_def[ind.coords,]$start.1[1], fireSeasonMedian_def[ind.coords,]$end.1[1]),
                             seq(fireSeasonMedian_def[ind.coords,]$start.2[1], fireSeasonMedian_def[ind.coords,]$end.2[1]))
                } else if (fireSeasonMedian_def[ind.coords,]$start.1[1] > fireSeasonMedian_def[ind.coords,]$end.1[1]){
                    meses = c(seq(1, fireSeasonMedian_def[ind.coords,]$end.1[1]), seq(fireSeasonMedian_def[ind.coords,]$start.1[1], 12),
                             seq(fireSeasonMedian_def[ind.coords,]$start.2[1], fireSeasonMedian_def[ind.coords,]$end.2[1]))
                } else {
                    meses = c(seq(1, fireSeasonMedian_def[ind.coords,]$end.2[1]), seq(fireSeasonMedian_def[ind.coords,]$start.2[1], 12),
                             seq(fireSeasonMedian_def[ind.coords,]$start.1[1], fireSeasonMedian_def[ind.coords,]$end.1[1]))
                }
            } else {
                next
            }
            ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)
                
            m = masked_ba_series.log[ind.meses, ind.coords]

            ind.years = as.numeric(substr(dates[ind.meses], 1, 4))
            n = cpc[which(cpc$Year %in% ind.years), c(meses + 1)]
            t = as.vector(t(n))

            corr.vector = c()
            for (j in 1:length(ind.coords)){
                l = m[,j]
                
                if (sum(l)>0){
                    if (cor.test(l,t[1:length(ind.meses)], method = 'pearson')$p.value < 0.05){
                        corr.vector = c(corr.vector, cor.test(l,t[1:length(ind.meses)], method = 'pearson')$estimate)
                    } else {
                        corr.vector = c(corr.vector, NA)
                    }   
                } else {
                    corr.vector = c(corr.vector, NA)
                }                
            }
            corr.df.monthly[ind.coords,]$ind.cpc = corr.vector                            
        }
    }

    arg.list <- list(col.regions = group.colors[11:1][-6], at = seq(-1, 1, 0.2), main = paste(name, "monthly correlation"))
    v <- corr.df.monthly$ind.cpc

    df1 <- cbind.data.frame(masked_coords, v)
    coordinates(df1) <- c(1,2)
    gridded(df1) <- TRUE
    arg.list[["sp.layout"]] <- list("sp.lines", coast.lines)
    arg.list[["obj"]] <- df1
    arg.list[["zcol"]] <- 1
    arg.list[["ylim"]] <- c(-90,90)
    arg.list[["xlim"]] <- c(-180,180)
    do.call("spplot", arg.list) %>% print()

    colnames(corr.df.monthly)[colnames(corr.df.monthly) == 'ind.cpc'] = name
    
    return (corr.df.monthly)
}




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



#' @title Annual correlation per cluster calculation
#' @description Obtains the correlation between the sum of the burned area of the fire season's months of each pixel of the cluster and the average of the indexes in these months for each cluster. Stores it in a dataframe and plots the results.
#' @param cpc Data frame containing the value of the cpc climate index per month. First column must be the year and the other ones have to be the months.
#' @param name Name of the climate index
#' @param corr.df Dataframe with the same form as masked_coords
#' @param mode Type of fire season whose correlation are we going to calculate. It could be 'unimodal' (by default), 'bimodal1' (main fire season of bimodal fire seasons) or 'bimodal2' (secondary fire season of bimodal fire seasons)
#' @param threshold Value used as a threshold for the plots. Only pixels which abs(cor) > threshold are plotted
#' @param plot Boolean deciding if the plot is shown
#' @return A copy of corr.df dataframe with a new column with the correlation
corr.annual.clus <- function(cpc, name, corr.df, mode = 'unimodal'){    
    
    if (mode == 'unimodal'){
        form = 1
    } else if (mode == 'bimodal1'){
        form = 2
    } else if (mode == 'bimodal2'){
        form = 2
        fireSeasonMedian_def[which(fireSeasonMedian_def$form == form),]$start.1 = fireSeasonMedian_def[which(fireSeasonMedian_def$form == form),]$start.2
        fireSeasonMedian_def[which(fireSeasonMedian_def$form == form),]$end.1 = fireSeasonMedian_def[which(fireSeasonMedian_def$form == form),]$end.2
    }
    
    corr.df$ind.cpc = NA
    
    for (biome in 1:13){
        clusters = sort(unique(fireSeasonMedian_def[which(fireSeasonMedian_def$form == form 
                                                     & fireSeasonMedian_def$BIOME == biome),]$cl))
        n.clusters = length(clusters)
        if (n.clusters < 1){
            next
        }
        for (cl in 1:n.clusters){
            clus = clusters[cl]         
            ind.coords = which(fireSeasonMedian_def$form == form & fireSeasonMedian_def$BIOME == biome & fireSeasonMedian_def$cl == clus)
            if (fireSeasonMedian_def[ind.coords,]$start.1[1] == fireSeasonMedian_def[ind.coords,]$end.1[1]){
                meses = fireSeasonMedian_def[ind.coords,]$start.1[1]
                ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)
                
                if (meses > 4){
                    ind.meses = ind.meses[-length(meses)]
                }           

            } else if (fireSeasonMedian_def[ind.coords,]$start.1[1] < fireSeasonMedian_def[ind.coords,]$end.1[1]){
                meses = seq(fireSeasonMedian_def[ind.coords,]$start.1[1], fireSeasonMedian_def[ind.coords,]$end.1[1])
                ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)
                    
                ind.meses = ind.meses[1:(length(meses) * floor(length(ind.meses)/length(meses)))]

            } else {
                meses = c(seq(1, fireSeasonMedian_def[ind.coords,]$end.1[1]), seq(fireSeasonMedian_def[ind.coords,]$start.1[1], 12))
                
                ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)
                
                if (fireSeasonMedian_def[ind.coords,]$end.1[1] <= 4){
                    ind.meses = ind.meses[(fireSeasonMedian_def[ind.coords,]$end.1[1]+1):(length(ind.meses))]
                } else {
                    ind.meses = ind.meses[(fireSeasonMedian_def[ind.coords,]$end.1[1]+1):(length(ind.meses)-4-(13-fireSeasonMedian_def[ind.coords,]$start.1[1]))]
                }
            }
                
            m = masked_ba_series.log[ind.meses, ind.coords]                

            # Obtaining burned area time series as the sum of the burned area during the fire season of all the points in the cluster
            ba.serie <- get.ba.serie(m, meses, ind.meses, ind.coords)   

            # Obtaining climate indexes time series
            cpc.serie <- get.cpc.serie(fireSeasonMedian_def, ind.coords, list(cpc), meses, ind.meses, df = as.data.frame(ba.serie))

            if (cor.test(ba.serie, cpc.serie[,-1], method = 'pearson')$p.value < 0.05){
                corr = cor.test(ba.serie, cpc.serie[,-1], method = 'pearson')$estimate
            } else {
                corr = NA
            }                        

            corr.df[ind.coords,]$ind.cpc = corr
        }
    }

    # Plot results
    arg.list <- list(col.regions = group.colors[11:1][-6],
                      at = seq(-1, 1, 0.2), main = paste(name, mode, "correlation per cluster"))
    v <- corr.df$ind.cpc

    df1 <- cbind.data.frame(masked_coords, v)
    coordinates(df1) <- c(1,2)
    gridded(df1) <- TRUE
    arg.list[["sp.layout"]] <- list("sp.lines", coast.lines)
    arg.list[["obj"]] <- df1
    arg.list[["zcol"]] <- 1
    arg.list[["ylim"]] <- c(-90,90)
    arg.list[["xlim"]] <- c(-180,180)
    do.call("spplot", arg.list) %>% print()

    if (mode == 'bimodal2'){
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = paste(name, '.2', sep = '')
    } else if (mode == 'bimodal1'){
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = paste(name, '.1', sep = '')
    } else {
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = name
    }
    
    return (corr.df)
}
