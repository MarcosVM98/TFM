library(rgdal)
library(raster)
library(sp)
library(magrittr)
library(RColorBrewer)

## 1. SPATIAL OVERLAY ----------------------------------------------------------
# Load reference BA data
load("ba_dataframe.Rdata", verbose = TRUE)
# Geo-reference of BA grid
spoints <- SpatialPoints(coords)
projection(spoints) <- CRS("+proj=longlat +ellps=WGS84")

# Read biome shapefile <https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world>
a <- readOGR("ignore/official_teow/official/wwf_terr_ecos.shp")
projection(a) <- CRS("+proj=longlat +ellps=WGS84")

# Spatial Overlay
ov <- over(spoints, a)

## 2. MAP OF BIOMES ON BA GRID -------------------------------------------------

# Create distinct colors
n <- 16
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                            rownames(qual_col_pals)))
set.seed(6) # Seed for random color choice
model.colors <- sample(col_vector, n)

# Plot biome map as SpatialGridDataFrame
df <- cbind.data.frame(coords, ov)
coordinates(df) <- c(1,2)
df$BIOME <- as.factor(df$BIOME)
gridded(df) <- TRUE
spplot(df, zcol = "BIOME", col.regions = model.colors)

# Legend of biomes
legend.biomes <- matrix(data = c("1", "Tropical and Subtropical Moist Broadleaf Forests", "TrMoistBrFor",
                                 "2", "Tropical and Subtropical Dry Broadleaf Forests", "TrDryBrFor",
                                 "3", "Tropical Conifer Forests", "TrConFor",
                                 "4", "Temperate Broadleaf and Mixed Forests", "TemBrFor",
                                 "5", "Temperate Conifer Forests", "TemConFor",
                                 "6", "Boreal Forests/Taiga", "Taiga",
                                 "7", "Tropical and Subtropical Grasslands, Savannas and Shrublands", "TrGrass",
                                 "8", "Temperate Grasslands, Savannas and Shrublands", "TemGrass",
                                 "9", "Flooded Grasslands and Savannas","FlGrass",
                                 "10","Montane Grasslands and Shrublands","MnGrass",
                                 "11", "Tundra", "Tundra",
                                 "12", "Meditarranean Forests, Woodlands and Scrub", "Med",
                                 "13", "Desert and Xeric Shrublands", "Desert",
                                 "14", "Mangroves", "Mangroves"),
                        ncol = 3, byrow = TRUE) %>% data.frame()

names(legend.biomes) <- c("ID", "Name", "Label")

## 3. SAVING FINAL DATASET -----------------------------------------------------

# biomes <- cbind.data.frame(coords, "BIOME" = df@data$BIOME)
# save(legend.biomes, biomes, file = "biome_dataframe.Rdata")





