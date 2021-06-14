# Functions used for obtaining the correlation between burned area during the fire season and climate indexes.

# Packages required
library(tidyr)
library(dplyr)
require(sp)
require(magrittr)
require(RColorBrewer)

group.colors <- colorRampPalette(c(brewer.pal(8, "Dark2"), brewer.pal(8, "Accent")))

# Data required
load("scripts/worldmap.Rdata", verbose = T)
load("data/ba_mon_time_series_masked.Rdata", verbose = T)
#load("fireSeasonMedian_def.Rdata", verbose = T)
load("fireSeasonPer75_def.Rdata", verbose = T)

fireSeasonMedian_def = fireSeasonPer75_def


#' @title Annual correlation calculation
#' @description Obtains the correlation between the sum of the burned area of the fire season's months and the average of the indexes in these months for each pixel. Stores it in a dataframe and plots the results.
#' @param cpc Data frame containing the value of the cpc climate index per month. First column must be the year and the other ones have to be the months.
#' @param name Name of the climate index
#' @param corr.df Dataframe with the same form as masked_coords
#' @param mode Type of fire season whose correlation are we going to calculate. It could be 'unimodal' (by default), 'bimodal1' (main fire season of bimodal fire seasons) or 'bimodal2' (secondary fire season of bimodal fire seasons)
#' @param threshold Value used as a threshold for the plots. Only pixels which abs(cor) > threshold are plotted
#' @param plot Boolean deciding if the plot is shown
#' @return A copy of corr.df dataframe with a new column with the correlation
corr.annual <- function(cpc, name, corr.df, mode = 'unimodal', threshold = 0.4, plot = T){    
    
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
        for (cl in 1:n.clusters){
            clus = clusters[cl]         
            ind.coords = which(fireSeasonMedian_def$form == form & fireSeasonMedian_def$BIOME == biome & fireSeasonMedian_def$cl == clus)
            if (fireSeasonMedian_def[ind.coords,]$start.1[1] == fireSeasonMedian_def[ind.coords,]$end.1[1] & length(ind.coords) > 0){
                meses = fireSeasonMedian_def[ind.coords,]$start.1[1]
                ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)
                
                if (meses > 4){
                    ind.meses = ind.meses[-length(meses)]
                }            
                
                m = masked_ba_series[ind.meses, ind.coords]
                
                ind.years = as.numeric(substr(dates[ind.meses], 1, 4))
                t = cpc[which(cpc$Year %in% ind.years), c(meses + 1)]

                corr.vector = c()
                for (j in 1:length(ind.coords)){                    
                    l = m[,j]
                    
                    if (sum(l)>0){
                        if (cor.test(l, t, method = 'pearson')$p.value < 0.05){
                            corr.vector = c(corr.vector, cor.test(l, t, method = 'pearson')$estimate)
                        } else {
                            corr.vector = c(corr.vector, NA)
                        }    
                    } else {
                        corr.vector = c(corr.vector, NA)
                    }                
                }
                corr.df[ind.coords,]$ind.cpc = corr.vector

            } else if (fireSeasonMedian_def[ind.coords,]$start.1[1] < fireSeasonMedian_def[ind.coords,]$end.1[1] & length(ind.coords) > 0){
                meses = seq(fireSeasonMedian_def[ind.coords,]$start.1[1], fireSeasonMedian_def[ind.coords,]$end.1[1])
                ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)
                    
                ind.meses = ind.meses[1:(length(meses) * floor(length(ind.meses)/length(meses)))]
                
                m = masked_ba_series[ind.meses, ind.coords]

                # cpc
                ind.years = as.numeric(substr(dates[ind.meses], 1, 4))
                n = cpc[which(cpc$Year %in% ind.years), c(meses + 1)]
                t = apply(n, 1, mean)
                

                corr.vector = c()
                for (j in 1:length(ind.coords)){
                    l = c()
                    for (i in 1:(length(ind.meses)/length(meses))){
                        l = c(l, sum(m[((i-1)*(length(meses))+1):(i*length(meses)),j]))
                    }
                    
                    if (sum(l)>0){
                        if (cor.test(l, t, method = 'pearson')$p.value < 0.05){
                            corr.vector = c(corr.vector, cor.test(l, t, method = 'pearson')$estimate)
                        } else {
                            corr.vector = c(corr.vector, NA)
                        }    
                    } else {
                        corr.vector = c(corr.vector, NA)
                    }                
                }
                corr.df[ind.coords,]$ind.cpc = corr.vector

            } else if (length(ind.coords) > 0) {
                meses = c(seq(1, fireSeasonMedian_def[ind.coords,]$end.1[1]), seq(fireSeasonMedian_def[ind.coords,]$start.1[1], 12))
                
                ind.meses = which(as.numeric(substr(dates, 6, 7)) %in% meses)
                
                if (fireSeasonMedian_def[ind.coords,]$end.1[1] <= 4){
                    ind.meses = ind.meses[(fireSeasonMedian_def[ind.coords,]$end.1[1]+1):(length(ind.meses))]
                } else {
                    ind.meses = ind.meses[(fireSeasonMedian_def[ind.coords,]$end.1[1]+1):(length(ind.meses)-4
                                                                                          -(13-fireSeasonMedian_def[ind.coords,]$start.1[1]))]
                }
                
                m = masked_ba_series[ind.meses, ind.coords]

                ind.years = as.numeric(substr(dates[ind.meses], 1, 4))
                n = cpc[which(cpc$Year %in% ind.years), c(1,(meses + 1))]
                
                n = as.vector(t(n))
                t = c()
                for (i in 1:(length(ind.meses)/length(meses))){
                    t = c(t, mean(n[(fireSeasonMedian_def[ind.coords,]$end.1[1] + 1 + (i-1)*length(meses)):
                                    (fireSeasonMedian_def[ind.coords,]$end.1[1] + i*length(meses))]))
                }

                corr.vector = c()
                for (j in 1:length(ind.coords)){
                    l = c()
                    for (i in 1:(length(ind.meses)/length(meses))){
                        l = c(l, sum(m[((i-1)*(length(meses)) + 1):(i*length(meses)),j]))
                    }
                    
                    if (sum(l)>0){
                        if (cor.test(l, t, method = 'pearson')$p.value < 0.05){
                            corr.vector = c(corr.vector, cor.test(l, t, method = 'pearson')$estimate)
                        } else {
                            corr.vector = c(corr.vector, NA)
                        }    
                    } else {
                        corr.vector = c(corr.vector, NA)
                    }
                }
                corr.df[ind.coords,]$ind.cpc = corr.vector
            }
        }
    }

    if (plot == T){
        arg.list <- list(col.regions = group.colors(11),
                          at = seq(threshold, 1, 0.1), main = paste(name, mode, " correlation"))
        v <- corr.df$ind.cpc
        v[which(abs(v) < threshold)] <- NA

        df1 <- cbind.data.frame(masked_coords, abs(v))
        coordinates(df1) <- c(1,2)
        gridded(df1) <- TRUE
        arg.list[["sp.layout"]] <- list("sp.lines", coast.lines)
        arg.list[["obj"]] <- df1
        arg.list[["zcol"]] <- 1
        arg.list[["ylim"]] <- c(-90,90)
        arg.list[["xlim"]] <- c(-180,180)
        do.call("spplot", arg.list) %>% print()
    }

    if (mode == 'bimodal2'){
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = paste(name, '.2', sep = '')
    } else if (mode == 'bimodal1'){
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = paste(name, '.1', sep = '')
    } else {
        colnames(corr.df)[colnames(corr.df) == 'ind.cpc'] = name
    }
    
    return (corr.df)
}


#' @title Monthly correlation calculation
#' @description Obtains the correlation bbetween burned arrea of each fire season's month and the value of the indexes in this month for each pixel. Stores it in a dataframe and plots the results.
#' @param cpc Data frame containing the value of the cpc climate index per month. First column must be the year and the other ones have to be the months.
#' @param name Name of the climate index
#' @param corr.df.monthly Dataframe with the same form as masked_coords
#' @param threshold Value used as a threshold for the plots. Only pixels which abs(cor) > threshold are plotted
#' @return A copy of corr.df.monthly dataframe with a new column with the correlation
corr.monthly <- function(cpc, name, corr.df.monthly, threshold = 0.4){
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
                
            m = masked_ba_series[ind.meses, ind.coords]

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

    arg.list <- list(col.regions = group.colors(11), at = seq(threshold, 1, 0.1), main = paste(name, "monthly correlation"))
    v <- corr.df.monthly$ind.cpc
    v[which(abs(v) < threshold)] <- NA

    df1 <- cbind.data.frame(masked_coords, abs(v))
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