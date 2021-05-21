library("mclust")
require(sp)
require(magrittr)
library(dplyr)
require(RColorBrewer)
group.colors <- colorRampPalette(c(brewer.pal(8, "Dark2"), brewer.pal(8, "Accent")))


## Functions used for defining and analyzing the clusters and also for calculating its fire seasons

#' @title Best GMM model selection
#' @description Retain the highest BIC model after an arbitary number of GMM initializations
#' @param data Data frame
#' @param K maximun number of groups
#' @param n.inits Number of initializations. Default to 3
#' @return A mclust-class model
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
    message("Lowest BIC: ", round(max(bics), 2))
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
    arg.list[["obj"]] <- df1
    arg.list[["zcol"]] <- 1
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
        centroide = colMeans(df[which(clus$classification == i),])

        per.25 = c()
        per.75 = c()
        for (j in 1:12){
            per.25[j] = quantile(df[which(clus$classification == i),j], prob=0.25)
            per.75[j] = quantile(df[which(clus$classification == i),j], prob=0.75)
        }

        plot(per.75, col = 'blue', type = 'l', xlab = "Month", ylab = "Burned Area",
            main = paste('Cluster', toString(i)), sub = paste('Size', toString(sum(clus$classification == i))))
        lines(1:12, centroide, col = "red")
        lines(1:12, per.25, col = "green")        
    }
}


#' @title Fire season calculation 
#' @description Calculates the fire season of each of the clusters of the mclust model
#' @param df Data frame containing the no Nas burned area data
#' @param clus A mclust-class model
#' @return List containing the fire season of each of the clusters

fireSeason = function(df, clus){
    fire.season = list()
    for (i in 1:clus$G){
        centroide = colMeans(df[which(clus$classification == i),])
        anual = as.data.frame(centroide)
        anual$mes = 1:12
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


## Function for performing the clustering and calculating the fire seasons. It uses the four previous functions and the data stored in "biome_dataframe.Rdata" and in "ba_dataframe.Rdata"

load("biome_dataframe.Rdata")
load("ba_dataframe.Rdata")

df.log = log1p(df)

#' @title Fire season calculation of each of the clusters of the biome 
#' @description Performs the clustering of the points with the biome, plots the spatial location and the centroid for each cluster and calculates the fire season of each cluster
#' @param biome Number between 1 and 13. Each number corresponds to one particular biome as you can see in the dataframe legend.biomes
#' @param min.size Minimum size of each the clusters 
#' @return List containing the fire season of each of the clusters

biome.clustering <- function (biome, min.size = 5){
    # Select the appropiate subset of the data
    ind.coords.biome = which(biomes$BIOME == biome)
    df.log.biome = df.log[ind.coords.biome,]
    coords.biome = coords[ind.coords.biome,]
    
    # Find the NAs points
    naind.biome <- which(is.na(df.log.biome), arr.ind = TRUE)
    
    K = 9
    minsize = 0
    
    print(toString(legend.biomes$Name[biome]))
    
    if (length(naind.biome) == 0){
        # Perform the clustering: clusters' size must be greater than min.size
        while (minsize < min.size){
            gmm.biome <- selectBestGMM(data = df.log.biome, K, n.inits = 3)
            size = c()
            for (i in 1:gmm.biome$gmm[[1]]$G){
                size = c(size, sum(gmm.biome$gmm[[1]]$classification == i))
            }
            minsize = min(size)
            K = K - 1
        }

        # Plot the clusters' location
        gmm.plot.biome <- plotClust.gmm(df.log.biome, coords.biome, gmm.biome$gmm[[1]])

        # Plot the clusters' centroids
        if (gmm.biome$gmm[[1]]$G < 5) {
            plotCentroids.gmm(df.log.biome, gmm.biome$gmm[[1]], c(2,2))
        } else {
            plotCentroids.gmm(df.log.biome, gmm.biome$gmm[[1]], c(3,3))
        }

        # Calculate the fire season for each cluster
        fireSeason.biome = fireSeason(df.log.biome, gmm.biome$gmm[[1]])
        
    } else {
        # Perform the clustering: clusters' size must be greater than min.size
        while (minsize < min.size){
            gmm.biome <- selectBestGMM(data = df.log.biome[-naind.biome,], K, n.inits = 3)
            size = c()
            for (i in 1:gmm.biome$gmm[[1]]$G){
                size = c(size, sum(gmm.biome$gmm[[1]]$classification == i))
            }
            minsize = min(size)
            K = K - 1
        }

        # Plot the clusters' location
        gmm.plot.biome <- plotClust.gmm(df.log.biome, coords.biome, gmm.biome$gmm[[1]], naind.biome)

        # Plot the clusters' centroids
        if (gmm.biome$gmm[[1]]$G < 5) {
            plotCentroids.gmm(df.log.biome[-naind.biome,], gmm.biome$gmm[[1]], c(2,2))
        } else {
            plotCentroids.gmm(df.log.biome[-naind.biome,], gmm.biome$gmm[[1]], c(3,3))
        }

        # Calculate the fire season for each cluster
        fireSeason.biome = fireSeason(df.log.biome[-naind.biome,], gmm.biome$gmm[[1]])
    }
        
    return (fireSeason.biome)
}
