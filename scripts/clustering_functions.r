## Functions used for defining and analyzing the clusters and also for calculating its fire seasons

# Required packages
library("mclust")
require(sp)
require(magrittr)
library(dplyr)
require(RColorBrewer)

# Data loading
load("scripts/worldmap.Rdata", verbose = T)
load("data/Fire/ba_mon_clim_masked_df.Rdata", verbose = T)
load("data/Fire/biome_dataframe_masked.Rdata", verbose = T)


group.colors <- colorRampPalette(c(brewer.pal(8, "Dark2"), brewer.pal(8, "Accent")))#a color palette for the plots


#' @title Best GMM model selection
#' @description Retain the highest BIC model after an arbitary number of GMM initializations
#' @param data Data frame
#' @param K maximun number of groups
#' @param n.inits Number of initializations. Default to 3
#' @return A list of mclust-class models and the index of the list of the model with the highest BIC
#' @importFrom mclust Mclust
selectBestGMM <- function(data, K, n.inits = 3) {
    i <- 1
    gmm.list <- rep(list(bquote(), n.inits))
    while (i <= n.inits) {
        set.seed(i)
        gmm.list[[i]] <- Mclust(data, G = 1:K)
        i <- i + 1
    }
    bics <- numeric(length(gmm.list))
    for (i in 1:length(gmm.list)) {
        bics[i] <- gmm.list[[i]]$bic
    }
    ind = which.max(bics)
    gmm <- gmm.list[[ind]]
    message("Highest BIC: ", round(max(bics), 2))
    message("Number of clusters:", gmm$G)
    for (i in 1:n.inits){
        cat(gmm.list[[i]]$G, gmm.list[[i]]$bic, gmm.list[[i]]$icl, '\n')
    }
    return (list("gmm" = gmm.list, "ind" = ind))
}


#' @title Clusters' spatial location plot
#' @description Plots the clusters' spatial location in a map
#' @param df Data frame containing the burned area data
#' @param coords Data frame containing the spatial coordinates
#' @param clus A mclust-class model
#' @param naind Vector containing the indexes where there are NAs values in df
#' @return Data of the plot
plotClust.gmm <- function(df, coords, clus, naind = c(), ...) {
    arg.list <- list(col.regions = group.colors(clus$G),
                      at = seq(0, clus$G, 1), main = paste(toString(clus$G), "clusters"))
    v <- numeric(nrow(df))
    if (length(naind) == 0){
        v <- clus$classification
    } else {
        v[naind] <- NA
        v[-naind] <- clus$classification
    }
    df1 <- cbind.data.frame(coords, v)
    coordinates(df1) <- c(1,2)
    gridded(df1) <- TRUE
    arg.list[["sp.layout"]] <- list("sp.lines", coast.lines)
    arg.list[["obj"]] <- df1
    arg.list[["zcol"]] <- 1
    arg.list[["ylim"]] <- c(-90,90)
    arg.list[["xlim"]] <- c(-180,180)
    do.call("spplot", arg.list) %>% print()
    return(df1@data)
}




#' @title Clusters' centroids plot
#' @description Plots the clusters' centroids and the 25 and 75 percentiles
#' @param df Data frame containing the burned area data without NAs
#' @param clus A mclust-class model
#' @param par A two-dimentional vector containing the dimensions of the grid
plotCentroids.gmm <- function(df, clus, par) {

    par(mfrow = par)

    for (i in 1:clus$G){
        med = apply(df[which(clus$classification == i),], 2, median)
        per.25 = vector("numeric", length = 12)
        per.75 = vector("numeric", length = 12)
        
        for (j in 1:12){
            per.25[j] = quantile(df[which(clus$classification == i),j], prob=0.25)
            per.75[j] = quantile(df[which(clus$classification == i),j], prob=0.75)
        }

        plot(per.75, col = 'blue', type = 'l', xlab = "Month", ylab = "Burned Area",
            main = paste('Cluster', toString(i)), sub = paste('Size', toString(sum(clus$classification == i))))
        lines(1:12, med, col = "red")
        lines(1:12, per.25, col = "green")
    }
}


#' @title Fire season calculation
#' @description Calculates the fire season of each of the clusters of the mclust model using the 75th percentile and the 80% of burned area as the threshold
#' @param df Data frame containing the no Nas burned area data
#' @param clus A mclust-class model
#' @return List containing the fire season of each of the clusters
fireSeason = function(df, clus){
    fire.season = list()
    for (i in 1:clus$G){
        centroide = apply(df[which(clus$classification == i),], 2, quantile, prob = 0.75)
        anual = as.data.frame(centroide)
        anual$mes = 1:12#if January is the first month
        anual = arrange(anual, -centroide)
        total = sum(anual$centroide)
        suma = 0
        mes = 1
        f.season = c()
        while (suma < 0.8 * total){
            suma = suma + anual$centroide[mes]
            mes = mes + 1
        }
        f.season = sort(anual[1:(mes-1),]$mes)
        cat(f.season, "\n")
        fire.season[[i]] = f.season
    }
    return (fire.season)
}



# Data is log-tranformed
df.log = log1p(df_masked)



#' @title Fire season calculation of each of the clusters of the biome
#' @description Performs the clustering of the points of the biome, plots the spatial location and the centroid for each cluster and calculates the fire season of each cluster
#' @param biome Number between 1 and 13. Each number corresponds to one particular biome as you can see in the dataframe legend.biomes
#' @param min.size Minimum size of each the clusters in percentage of the whole numbers of points of the biome
#' @return List containing the fire season of each of the clusters and the mclust clustering model used for calculating the fire seasons
biome.clustering <- function (biome, min.size = 5){
    # Select the appropiate subset of the data
    ind.coords.biome = which(biomes$BIOME == biome)
    df.log.biome = df.log[ind.coords.biome,]
    coords.biome = masked_coords[ind.coords.biome,]

    # Find the NAs points
    naind.biome <- which(is.na(df.log.biome), arr.ind = TRUE)
    
    # Find the smallest number of points that a cluster can have and the maximum number of clusters that are allowed
    max.n.clusters = floor(100 / min.size)
    min.size = min.size * (length(coords.biome[,1]) - length(naind.biome)) / 100

    K = 10 - 5
    minsize = 0

    print(toString(legend.biomes$Name[biome]))

    if (length(naind.biome) == 0){
        # Perform the clustering: the size of clusters must be greater than min.size        
        G = K
        while (G == K & K < max.n.clusters){
            K = K + 5
            gmm.biome <- selectBestGMM(data = df.log.biome, K = min(K, max.n.clusters), n.inits = 3)
            ind = gmm.biome$ind
            G = gmm.biome$gmm[[ind]]$G                
        }            
        size = vector("integer", length = gmm.biome$gmm[[ind]]$G)
        for (i in 1:gmm.biome$gmm[[ind]]$G){
            size[i] = sum(gmm.biome$gmm[[ind]]$classification == i)
        }
        minsize = min(size, na.rm = T)
        while (minsize < min.size){
            K = K - 1
            gmm.biome <- selectBestGMM(data = df.log.biome, K = K, n.inits = 3)
            ind = gmm.biome$ind
            size = vector("integer", length = gmm.biome$gmm[[ind]]$G)
            for (i in 1:gmm.biome$gmm[[ind]]$G){
                size[i] = sum(gmm.biome$gmm[[ind]]$classification == i)
            }
            minsize = min(size, na.rm = T)
        }

        # Plot the clusters' location
        gmm.plot.biome <- plotClust.gmm(df.log.biome, coords.biome, gmm.biome$gmm[[ind]])

        # Plot the clusters' centroids
        if (gmm.biome$gmm[[ind]]$G < 5) {
            plotCentroids.gmm(df.log.biome, gmm.biome$gmm[[ind]], c(2,2))
        } else {
            plotCentroids.gmm(df.log.biome, gmm.biome$gmm[[ind]], c(3,3))
        }

        # Calculate the fire season for each cluster
        fireSeason.biome = fireSeason(df.log.biome, gmm.biome$gmm[[ind]])

    } else {
        # Perform the clustering: the size of clusters must be greater than min.size
        G = K
        while (G == K & K < max.n.clusters){
            K = K + 5
            gmm.biome <- selectBestGMM(data = df.log.biome, K = min(K, max.n.clusters), n.inits = 3)
            ind = gmm.biome$ind
            G = gmm.biome$gmm[[ind]]$G                
        }            
        size = vector("integer", length = gmm.biome$gmm[[ind]]$G)
        for (i in 1:gmm.biome$gmm[[ind]]$G){
            size[i] = sum(gmm.biome$gmm[[ind]]$classification == i)
        }
        minsize = min(size, na.rm = T)
        while (minsize < min.size){
            K = K - 1
            gmm.biome <- selectBestGMM(data = df.log.biome[-naind.biome,], K = K, n.inits = 3)
            ind = gmm.biome$ind
            size = vector("integer", length = gmm.biome$gmm[[ind]]$G)
            for (i in 1:gmm.biome$gmm[[ind]]$G){
                size[i] = sum(gmm.biome$gmm[[ind]]$classification == i)
            }
            minsize = min(size, na.rm = T)
        }

        # Plot the clusters' location
        gmm.plot.biome <- plotClust.gmm(df.log.biome, coords.biome, gmm.biome$gmm[[ind]], naind.biome)

        # Plot the clusters' centroids
        if (gmm.biome$gmm[[ind]]$G < 5) {
            plotCentroids.gmm(df.log.biome[-naind.biome,], gmm.biome$gmm[[ind]], c(2,2))
        } else {
            plotCentroids.gmm(df.log.biome[-naind.biome,], gmm.biome$gmm[[ind]], c(3,3))
        }

        # Calculate the fire season for each cluster
        fireSeason.biome = fireSeason(df.log.biome[-naind.biome,], gmm.biome$gmm[[ind]])
    }

    return (list("fs" = fireSeason.biome, "gmm" = gmm.biome$gmm[[ind]]))
}



#' @title Construction of a data frame including the fire season for each pixel
#' @description Builds a data frame with the fire season information. For each point the data frame contains its coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season
#' @param fireSeasonsComplete List containing the elements that are returned by biome.clustering function
#' @return A data frame containing coordinates, biome, cluster, start and end months of the fire season, start and end months of the secondary fire season if exists and form of the fire season for each pixel
fireSeasonsDfConstruction <- function (fireSeasonsComplete){
    # Store the fire season of each pixel in the dataframe "fireSeasonMedian_def"
    cl = rep(NA ,length(biomes$x))# cluster of the pixel
    start.1 = rep(NA ,length(biomes$x))# start month of the main fire season
    end.1 = rep(NA ,length(biomes$x))# end month of the main fire season
    start.2 = rep(NA ,length(biomes$x))# start month of the secondary fire season
    end.2 = rep(NA ,length(biomes$x))# end month of the secondary fire season
    form = rep(NA ,length(biomes$x))# form of the fire season
    fireSeasonMedian_def = cbind.data.frame(biomes, cl, start.1, end.1, start.2, end.2, form)

    for (b in 1:13){
        ind.coords.biome = which(biomes$BIOME == b)
        biom = fireSeasonsComplete[[b]]
        fireSeasonMedian_def[ind.coords.biome,]$cl = as.factor(biom$gmm$classification)

        for (i in 1:length(biom$fs)){
            if (length(biom$fs[[i]]) > 1){# If the fire season has more than 1 month
                dif = diff(biom$fs[[i]])
                if (sum(dif > 1) == 0 | (sum(dif > 1) == 1 & 12 %in% biom$fs[[i]] & 1 %in% biom$fs[[i]])){# Unimodal fire season          
                    j = 1
                    if (12 %in% biom$fs[[i]] & 1 %in% biom$fs[[i]]){# if the fire season starts in one year and finishes in the following one
                        j = length(dif)
                        while (dif[j] == 1 & j > 1){
                            j = j - 1
                        }
                        k = 1
                        while (dif[k] == 1 & k < length(dif)){
                            k = k + 1
                        }
                        cat(i, biom$fs[[i]][j+1], biom$fs[[i]][k], ' UNIMODAL ', biom$fs[[i]], '\n')
                        fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$start.1 = biom$fs[[i]][j+1]
                        fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$end.1 = biom$fs[[i]][k]            
                        fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$form = 1
                    } else {# if the whole fire season is in one year
                        k = 1
                        while (dif[k] == 1 & k < length(dif)){
                            k = k + 1
                        }
                        cat(i, biom$fs[[i]][j], biom$fs[[i]][k+1], ' UNIMODAL ', biom$fs[[i]], '\n')
                        fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$start.1 = biom$fs[[i]][j]
                        fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$end.1 = biom$fs[[i]][k+1]            
                        fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$form = 1
                    }
                } else if (sum(dif > 1) == 1 | (sum(dif > 1) == 2 & 12 %in% biom$fs[[i]] & 1 %in% biom$fs[[i]])) {# bimodal fire season
                    centroide = apply(df.log[ind.coords.biome,][which(biom$gmm$classification == i),], 2, median)
                    j1 = 1
                    if (12 %in% biom$fs[[i]] & 1 %in% biom$fs[[i]]){# if one the fire seasons starts in one year and finishes in the following one
                        j1 = length(dif)
                        while (dif[j1] == 1 & j1 > 1){
                            j1 = j1 - 1
                        }
                        k1 = 1
                        while (dif[k1] == 1 & k1 < length(dif)){
                            k1 = k1 + 1
                        }

                        # Deciding which is the main fire season
                        if (sum(centroide[1:biom$fs[[i]][k1]]) + sum(centroide[biom$fs[[i]][j1+1]:12]) >= sum(centroide[biom$fs[[i]][k1+1]:biom$fs[[i]][j1]])){
                            cat(i, ' F1 ', biom$fs[[i]][j1+1], biom$fs[[i]][k1], ' F2 ', biom$fs[[i]][k1+1], biom$fs[[i]][j1], ' BIMODAL1 ', biom$fs[[i]], '\n')
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$start.1 = biom$fs[[i]][j1+1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$end.1 = biom$fs[[i]][k1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$start.2 = biom$fs[[i]][k1+1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$end.2 = biom$fs[[i]][j1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$form = 2
                        } else {
                            cat(i, ' F1 ', biom$fs[[i]][k1+1], biom$fs[[i]][j1], ' F2 ', biom$fs[[i]][j1+1], biom$fs[[i]][k1], ' BIMODAL2 ', biom$fs[[i]], '\n')
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$start.1 = biom$fs[[i]][k1+1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$end.1 = biom$fs[[i]][j1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$start.2 = biom$fs[[i]][j1+1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$end.2 = biom$fs[[i]][k1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$form = 2
                        }               
                    } else {# If both fire seasons take place in one year
                        k1 = 1
                        while (dif[k1] == 1 & k1 < length(dif)){
                            k1 = k1 + 1
                        }
                        if (k1 == 1){
                            k1 = 1
                        }

                        # Deciding which is the main fire season
                        if (sum(centroide[biom$fs[[i]][1]:biom$fs[[i]][k1]]) >= sum(centroide[biom$fs[[i]][k1+1]:biom$fs[[i]][length(biom$fs[[i]])]])){
                            cat(i, ' F1 ', biom$fs[[i]][1], biom$fs[[i]][k1], ' F2 ', biom$fs[[i]][k1+1], biom$fs[[i]][length(biom$fs[[i]])], ' BIMODAL1 ', biom$fs[[i]], '\n')   
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$start.1 = biom$fs[[i]][1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$end.1 = biom$fs[[i]][k1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$start.2 = biom$fs[[i]][k1+1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$end.2 = biom$fs[[i]][length(biom$fs[[i]])]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$form = 2
                        } else {
                            cat(i, ' F1 ', biom$fs[[i]][k1+1], biom$fs[[i]][length(biom$fs[[i]])], ' F2 ', biom$fs[[i]][1], biom$fs[[i]][k1], ' BIMODAL2 ', biom$fs[[i]], '\n')
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$start.1 = biom$fs[[i]][k1+1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$end.1 = biom$fs[[i]][length(biom$fs[[i]])]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$start.2 = biom$fs[[i]][1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$end.2 = biom$fs[[i]][k1]
                            fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$form = 2
                        }
                    }            

                } else {# if the fire season has more than 2 fire seasons
                    cat(i, 'STRANGE FIRE SEASON', biom$fs[[i]], '\n')
                    fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$form = 3
                }
            } else {# if the fire season has only 1 month
                centroide = apply(df.log[ind.coords.biome,][which(biom$gmm$classification == i),], 2, median)
                # deciding if there is enough burned area to consider the fire season or not
                if (max(centroide, na.rm = T) < 0.2){
                    cat(i, biom$fs[[i]][1], biom$fs[[i]][1], ' PLANA ', biom$fs[[i]], '\n')
                    fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$form = 0
                } else {
                    cat(i, biom$fs[[i]][1], biom$fs[[i]][1], ' UNIMENSUAL ', biom$fs[[i]], '\n')
                    fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$start.1 = biom$fs[[i]][1]
                    fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$end.1 = biom$fs[[i]][1]            
                    fireSeasonMedian_def[ind.coords.biome,][which(fireSeasonMedian_def[ind.coords.biome,]$cl == i),]$form = 1
                }
            }    
        }
    }
    return (fireSeasonMedian_def)
}
