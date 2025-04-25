#Sweeney, Harrison and Vander Linden 2024. 
#European forest cover during the Holocene reconstructed from pollen records
# ---------------------------------------------------------

#There are six seperate scripts associated with this research:
# 1. Map_setup.R (this script)
# 2. Modern_pollen.R 
# 3. Model_import_data.R
# 4. Tree_model.R
# 5. Tree_recon.R
# 6. Figures.R

# ---------------------------------------------------------

############ Map setup ############ 

# This script sets boundary areas, sets up empty raster files, imports landcover 
# data, and generates a masked map of modern tree cover


# ---------------------------------------------------------
# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Set map extents
# 3. Load and crop European base map and Hengl pnv map 
# 4. Generate masked modern tree cover map


# ---------------------------------------------------------


# 1. Packages
# ---------------------------------------------------------
# ---------------------------------------------------------

# Install
# ---------------------------------------------------------
# install.packages("analogue")
# install.packages("betareg")
# install.packages("broom")
# install.packages("car")
# install.packages("exactextractr")
# install.packages("ggplotify")
# install.packages("ggpubr")
# install.packages("ggthemes")
# install.packages("grDevices")
# install.packages("gridExtra")
# install.packages("lmtest")
# install.packages("locfit")
# install.packages("parallel")
# install.packages("qmap")
# install.packages("raster")
# install.packages("remotes")
# install.packages("rnaturalearth")
# install.packages("raster")
# install.packages("rio")
# install.packages("sf")
# remotes::install_github("special-uor/smpds")
# install.packages("terra")
# install.packages("tidytable")
# install.packages("tidyverse")
# install.packages("vegan")
# install.packages("viridis")
# install.packages("WorldFlora")



# Library
# ---------------------------------------------------------
library(tidyverse)



# Paths
# ---------------------------------------------------------
# The script(s) currently run(s) with the following sub-folders:
# /data: any data that is used for the analysis
#   /input                          : unmodified data inputs
#       /copernicus_frac_cover       
#         /2015                     : landcover data from https://zenodo.org/records/3939038
#         /2016                     : landcover data from https://zenodo.org/records/3518026
#         /2017                     : landcover data from https://zenodo.org/records/3518036
#         /2018                     : landcover data from https://zenodo.org/records/3518038
#         /2019                     : landcover data from https://zenodo.org/records/3939050
#       /EEA_bioregions             : biogeographical regions map from https://www.eea.europa.eu/data-and-maps/figures/biogeographical-regions-in-europe-2
#       /Hengl                      : potential natural vegetation map from https://zenodo.org/records/10513520 
#       /Serge                      
#         /TERRA_RVresults_RPPs.st1 
#           /RV_mean_RPPs.st1       : Serge et al. (2023) mean vegetation data from https://data.indores.fr/dataset.xhtml?persistentId=doi:10.48579/PRO/J5GZUO
#       /SMPDS                      : updated information regarding SMPDS metadata (see SMPDSv2_updated_meta_info.csv)
#       /SPECIAL-EPD                : SPECIAL.EPD data from https://researchdata.reading.ac.uk/1295/
#       /ZANON
#        /forest_cover              : Zanon et al. (2017) tree cover data from https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2018.00253/full
#   /intermediate_output  
#     /vegetation                   : saved intermediate data
#       /copernicus                 : amalgamated maps based on tree cover data
#       /euro_map                   : maps of europe
#         /3035                     : crs 3035 shapefile
#         /4258                     : crs 4258 shapefile
#         /4326                     : crs 4326 shapefile
#       /hengl                      : european pnv, crs 3035
#       /SMPDS                      : adjusted modern pollen data
#       /raster_tree                : adjusted modern pollen data
# /figs                             : figures/tables outputted from the analysis
#   /supplement                     : supplementary figures/tables
#   /tree_cover_rasters             : rasters of tree cover reconstructions



# If the folders have not yet been created, the following code
# will generate these from the working directory, with the exception of
# the code file, where this and the other scripts should be saved. Input data
# should be saved into the relevant input filepath:

# dir.create("data")
# dir.create("data/input")
# dir.create("data/input/copernicus_frac_cover")
# dir.create("data/input/copernicus_frac_cover/2015")
# dir.create("data/input/copernicus_frac_cover/2016")
# dir.create("data/input/copernicus_frac_cover/2017")
# dir.create("data/input/copernicus_frac_cover/2018")
# dir.create("data/input/copernicus_frac_cover/2019")
# dir.create("data/input/EEA_bioregions")
# dir.create("data/input/Hengl")
# dir.create("data/input/Serge")
# dir.create("data/input/Serge/TERRA_RVresults_RPPs.st1")
# dir.create("data/input/Serge/TERRA_RVresults_RPPs.st1/RV_mean_RPPs.st1")
# dir.create("data/input//SMPDS")
# dir.create("data/input/SPECIAL-EPD")
# dir.create("data/input/ZANON")
# dir.create("data/input/ZANON/forest_cover")
# dir.create("data/intermediate_output")
# dir.create("data/intermediate_output/vegetation")
# dir.create("data/intermediate_output/vegetation/copernicus")
# dir.create("data/intermediate_output/vegetation/euro_map")
# dir.create("data/intermediate_output/vegetation/hengl")
# dir.create("data/intermediate_output/vegetation/SMPDS")
# dir.create("data/intermediate_output/vegetation/raster_tree")
# dir.create("figs")
# dir.create("figs/supplement")
# dir.create("figs/tree_cover_rasters")


# ---------------------------------------------------------





# 2. Set map extents
# ---------------------------------------------------------
# ---------------------------------------------------------
euro_extent_4258_ymin <- 34 
euro_extent_4258_ymax <- 73
euro_extent_4258_xmin <- -12
euro_extent_4258_xmax <- 45

euro_extent_4326_ymin <- 34 
euro_extent_4326_ymax <- 73
euro_extent_4326_xmin <- -12 
euro_extent_4326_xmax <- 45

boundary_area_4258 <- sf::st_bbox(c(xmin = euro_extent_4258_xmin, xmax = euro_extent_4258_xmax, ymin = euro_extent_4258_ymin, ymax = euro_extent_4258_ymax), crs = 4258)
boundary_area_4326 <- sf::st_bbox(c(xmin = euro_extent_4326_xmin, xmax = euro_extent_4326_xmax, ymin = euro_extent_4326_ymin, ymax = euro_extent_4326_ymax), crs = 4326)
boundary_area_3035 <- boundary_area_4258 %>% 
  sf::st_as_sfc() %>% 
  sf::st_transform(crs = 3035) %>% 
  sf::st_bbox()
boundary_area_4326 <- sf::st_bbox(c(xmin = euro_extent_4326_xmin, xmax = euro_extent_4326_xmax, ymin = euro_extent_4326_ymin, ymax = euro_extent_4326_ymax), crs = 4326)

# ---------------------------------------------------------


# 3. Load and crop European base map and Hengl pnv map
# ---------------------------------------------------------
# ---------------------------------------------------------
global_map_4326 <- rnaturalearth::ne_download(scale = 50, type = "land", category = "physical", returnclass = "sf")
euro_map_4326 <- terra::vect(global_map_4326) %>% 
  terra::crop(., boundary_area_4326) %>% #crop to boundary area
  terra::as.data.frame(., geom = "WKT") %>% 
  sf::st_as_sf(., wkt = "geometry")
sf::st_write(euro_map_4326, "data/intermediate_output/vegetation/euro_map/4326/euro_map_4326.shp", append = FALSE) #save map

global_map_4258 <- sf::st_transform(global_map_4326, crs = 4258) #transform crs
euro_map_4258 <- terra::vect(global_map_4258) %>% 
  terra::crop(., boundary_area_4258) %>% #crop to boundary area
  terra::as.data.frame(., geom = "WKT") %>% 
  sf::st_as_sf(., wkt = "geometry")
sf::st_write(euro_map_4258, "data/intermediate_output/vegetation/euro_map/4258/euro_map_4258.shp", append = FALSE) #save map

global_map_3035 <- sf::st_transform(global_map_4326, crs = 3035) #transform crs
euro_map_3035 <- sf::st_crop(global_map_3035, boundary_area_3035) #crop to boundary area
sf::st_write(euro_map_3035, "data/intermediate_output/vegetation/euro_map/3035/euro_map_3035.shp", append = FALSE) #save map

euro_map_3035_rast_50 <- raster::raster(euro_map_3035, res = c(50000, 50000)) #save raster at 50km resolution
euro_map_3035_rast_7 <- raster::raster(euro_map_3035, res = c(7000, 7000)) #save raster at 10km resolution
euro_map_3035_rast_80 <- raster::raster(euro_map_3035, res = c(80000, 80000)) #save raster at 80km resolution

raster::writeRaster(euro_map_3035_rast_50, "data/intermediate_output/vegetation/euro_map/euro_map_3035_rast_50.tif", overwrite = TRUE)
raster::writeRaster(euro_map_3035_rast_7, "data/intermediate_output/vegetation/euro_map/euro_map_3035_rast_7.tif", overwrite = TRUE)
raster::writeRaster(euro_map_3035_rast_80, "data/intermediate_output/vegetation/euro_map/euro_map_3035_rast_80.tif", overwrite = TRUE)


europe_biomes_3035 <- terra::rast("data/input/Hengl/hengl_pnv.tif") %>% 
  terra::crop(., boundary_area_4326) %>% 
  terra::project(., "EPSG:3035", method = "near")
terra::writeRaster(europe_biomes_3035, "data/intermediate_output/vegetation/hengl/europe_biomes_3035.tif", overwrite = TRUE)

# ---------------------------------------------------------


# 4. Generate masked modern tree cover map
# ---------------------------------------------------------
# ---------------------------------------------------------
#Copernicus fractional cover
cop_file_list <- list.files("data/input/copernicus_frac_cover/", pattern = "\\.tif$", recursive = TRUE, full.names = TRUE) #filenames
f_cop_mask <- function(a,b,c){
  cop_file_list <- cop_file_list[stringr::str_detect(basename(cop_file_list), a)] #select cover files
  cop_cover <- terra::rast(cop_file_list) %>% #import and crop to extent
    terra::crop(., boundary_area_4326)
  cop_cover_max <- max(cop_cover, na.rm = TRUE) #identify the maximum value across layers
  cop_cover_mask <- terra::clamp(cop_cover_max, lower = c, values = FALSE) %>% #set all values less than 50% to NA then convert to 0
    terra::subst(., NA, 0)
  assign(paste0("cop_cover_", b, "_mask"), cop_cover_mask, envir = parent.frame()) #store masked files
}
  
#Mask where 1. Crop is 50% or greater; 2. Snow and ice cover is 50%; 3. Bare is 50%; 4. Built-up is 50%; 5. Permanent water is 50%; 6. Moss is 50%
#Load and crop to boundary
f_cop_mask("Crops","crop",49.5)
f_cop_mask("Crops","crop0.75",74.9)
f_cop_mask("Crops","crop0.65",64.9)
f_cop_mask("Crops","crop0.85",84.9)
f_cop_mask("Crops","crop0.25",24.9)
f_cop_mask("Crops","crop0.35",34.9)
f_cop_mask("Snow","ice",49.5)
f_cop_mask("Bare","bare",49.5)
f_cop_mask("BuiltUp","built",49.5)
f_cop_mask("PermanentWater","water",49.5)
f_cop_mask("MossLichen","moss",49.5)


#Combined mask layers
cop_mask_comb <- cop_cover_water_mask + cop_cover_built_mask + cop_cover_bare_mask + cop_cover_ice_mask + cop_cover_crop_mask + cop_cover_moss_mask    #add layers than set values below 50% to NA
cop_mask <- terra::clamp(cop_mask_comb, lower = 49.5, values = FALSE) #any layer above 50%
terra::writeRaster(cop_mask, "data/intermediate_output/vegetation/copernicus/cop_mask_100m.tif", overwrite = TRUE)

cop_mask_comb_not_crop <- cop_cover_water_mask + cop_cover_built_mask + cop_cover_bare_mask + cop_cover_ice_mask + cop_cover_moss_mask    #add layers than set values below 50% to NA
cop_mask_not_crop <- terra::clamp(cop_mask_comb_not_crop, lower = 49.5, values = FALSE) #any layer above 50%
terra::writeRaster(cop_mask_not_crop, "data/intermediate_output/vegetation/copernicus/cop_mask_not_crop_100m.tif", overwrite = TRUE)

cop_mask_comb_crop0.75 <- cop_cover_water_mask + cop_cover_built_mask + cop_cover_bare_mask + cop_cover_ice_mask + cop_cover_moss_mask +  cop_cover_crop0.75_mask  #add layers than set values below 50% to NA
cop_mask_crop0.75_mask <- terra::clamp(cop_mask_comb_crop0.75, lower = 49.5, values = FALSE) #any layer above 50%
terra::writeRaster(cop_mask_crop0.75_mask, "data/intermediate_output/vegetation/copernicus/cop_mask_crop0.75_100m.tif", overwrite = TRUE)

cop_mask_comb_crop0.65 <- cop_cover_water_mask + cop_cover_built_mask + cop_cover_bare_mask + cop_cover_ice_mask + cop_cover_moss_mask +  cop_cover_crop0.65_mask  #add layers than set values below 50% to NA
cop_mask_crop0.65_mask <- terra::clamp(cop_mask_comb_crop0.65, lower = 49.5, values = FALSE) #any layer above 50%
terra::writeRaster(cop_mask_crop0.65_mask, "data/intermediate_output/vegetation/copernicus/cop_mask_crop0.65_100m.tif", overwrite = TRUE)

cop_mask_comb_crop0.85 <- cop_cover_water_mask + cop_cover_built_mask + cop_cover_bare_mask + cop_cover_ice_mask + cop_cover_moss_mask +  cop_cover_crop0.85_mask  #add layers than set values below 50% to NA
cop_mask_crop0.85_mask <- terra::clamp(cop_mask_comb_crop0.85, lower = 49.5, values = FALSE) #any layer above 50%
terra::writeRaster(cop_mask_crop0.85_mask, "data/intermediate_output/vegetation/copernicus/cop_mask_crop0.85_100m.tif", overwrite = TRUE)

cop_mask_comb_crop0.35 <- cop_cover_water_mask + cop_cover_built_mask + cop_cover_bare_mask + cop_cover_ice_mask + cop_cover_moss_mask +  cop_cover_crop0.35_mask  #add layers than set values below 50% to NA
cop_mask_crop0.35_mask <- terra::clamp(cop_mask_comb_crop0.35, lower = 34.5, values = FALSE) #any layer above 35%
terra::writeRaster(cop_mask_crop0.35_mask, "data/intermediate_output/vegetation/copernicus/cop_mask_crop0.35_100m.tif", overwrite = TRUE)

cop_mask_comb_crop0.25 <- cop_cover_water_mask + cop_cover_built_mask + cop_cover_bare_mask + cop_cover_ice_mask + cop_cover_moss_mask +  cop_cover_crop0.25_mask  #add layers than set values below 50% to NA
cop_mask_crop0.25_mask <- terra::clamp(cop_mask_comb_crop0.25, lower = 24.5, values = FALSE) #any layer above 25%
terra::writeRaster(cop_mask_crop0.25_mask, "data/intermediate_output/vegetation/copernicus/cop_mask_crop0.25_100m.tif", overwrite = TRUE)

#Tree cover copfract
cop_file_list_tree <- cop_file_list[stringr::str_detect(basename(cop_file_list), "Tree")] #select cover files
cop_cover_tree <- terra::rast(cop_file_list_tree) %>% #import and crop to extent
  terra::crop(., boundary_area_4326)
cop_cover_tree_mean <- terra::mean(cop_cover_tree, na.rm = TRUE) #take average across years
terra::writeRaster(cop_cover_tree_mean, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_100m.tif", overwrite = TRUE)
cop_cover_tree_mean_3035 <- terra::project(cop_cover_tree_mean, "EPSG:3035")
terra::writeRaster(cop_cover_tree_mean_3035, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_3035_100m.tif", overwrite = TRUE)

cop_tree_cover_masked <- terra::mask(cop_cover_tree_mean, cop_mask, inverse = TRUE)
terra::writeRaster(cop_tree_cover_masked, "data/intermediate_output/vegetation/copernicus/cop_masked_tree_cover_100m.tif", overwrite = TRUE)
cop_tree_cover_masked_3035 <- terra::project(cop_tree_cover_masked, "EPSG:3035")
terra::writeRaster(cop_tree_cover_masked_3035, "data/intermediate_output/vegetation/copernicus/cop_masked_tree_cover_3035_100m.tif", overwrite = TRUE)

cop_tree_cover_masked_not_crop <- terra::mask(cop_cover_tree_mean, cop_mask_not_crop, inverse = TRUE)
terra::writeRaster(cop_tree_cover_masked_not_crop, "data/intermediate_output/vegetation/copernicus/cop_masked_not_crop_tree_cover_100m.tif", overwrite = TRUE)
cop_tree_cover_masked_not_crop_3035 <- terra::project(cop_tree_cover_masked_not_crop, "EPSG:3035")
terra::writeRaster(cop_tree_cover_masked_not_crop_3035, "data/intermediate_output/vegetation/copernicus/cop_masked_not_crop_tree_cover_3035_100m.tif", overwrite = TRUE)

cop_tree_cover_masked_crop0.75 <- terra::mask(cop_cover_tree_mean, cop_mask_crop0.75_mask, inverse = TRUE)
terra::writeRaster(cop_tree_cover_masked_crop0.75, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_masked_crop0.75_tree_cover_100m.tif", overwrite = TRUE)
cop_tree_cover_masked_crop0.75_3035 <- terra::project(cop_tree_cover_masked_crop0.75, "EPSG:3035")
terra::writeRaster(cop_tree_cover_masked_crop0.75_3035, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_masked_crop0.75_3035_100m.tif", overwrite = TRUE)

cop_tree_cover_masked_crop0.65 <- terra::mask(cop_cover_tree_mean, cop_mask_crop0.65_mask, inverse = TRUE)
terra::writeRaster(cop_tree_cover_masked_crop0.65, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_masked_crop0.65_tree_cover_100m.tif", overwrite = TRUE)
cop_tree_cover_masked_crop0.65_3035 <- terra::project(cop_tree_cover_masked_crop0.65, "EPSG:3035")
terra::writeRaster(cop_tree_cover_masked_crop0.65_3035, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_masked_crop0.65_3035_100m.tif", overwrite = TRUE)

cop_tree_cover_masked_crop0.85 <- terra::mask(cop_cover_tree_mean, cop_mask_crop0.85_mask, inverse = TRUE)
terra::writeRaster(cop_tree_cover_masked_crop0.85, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_masked_crop0.85_tree_cover_100m.tif", overwrite = TRUE)
cop_tree_cover_masked_crop0.85_3035 <- terra::project(cop_tree_cover_masked_crop0.85, "EPSG:3035")
terra::writeRaster(cop_tree_cover_masked_crop0.85_3035, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_masked_crop0.85_3035_100m.tif", overwrite = TRUE)

cop_tree_cover_masked_crop0.35 <- terra::mask(cop_cover_tree_mean, cop_mask_crop0.35_mask, inverse = TRUE)
terra::writeRaster(cop_tree_cover_masked_crop0.35, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_masked_crop0.35_tree_cover_100m.tif", overwrite = TRUE)
cop_tree_cover_masked_crop0.35_3035 <- terra::project(cop_tree_cover_masked_crop0.35, "EPSG:3035")
terra::writeRaster(cop_tree_cover_masked_crop0.35_3035, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_masked_crop0.35_3035_100m.tif", overwrite = TRUE)

cop_tree_cover_masked_crop0.25 <- terra::mask(cop_cover_tree_mean, cop_mask_crop0.25_mask, inverse = TRUE)
terra::writeRaster(cop_tree_cover_masked_crop0.25, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_masked_crop0.25_tree_cover_100m.tif", overwrite = TRUE)
cop_tree_cover_masked_crop0.25_3035 <- terra::project(cop_tree_cover_masked_crop0.25, "EPSG:3035")
terra::writeRaster(cop_tree_cover_masked_crop0.25_3035, "data/intermediate_output/vegetation/copernicus/cop_tree_cover_masked_crop0.25_3035_100m.tif", overwrite = TRUE)

#Crop cover
crop_cover_max <- terra::rast("data/input/copernicus_frac_cover//2019/PROBAV_LC100_global_v3.0.1_2019-nrt_Crops-CoverFraction-layer_EPSG-4326.tif" ) %>% 
  terra::crop(., boundary_area_4326) %>% 
  max(., na.rm = TRUE)
terra::writeRaster(crop_cover_max, "data/intermediate_output/vegetation/copernicus/crop_cover_max.tif", overwrite = TRUE)
crop_cover_max_3035 <- terra::project(crop_cover_max, "EPSG:3035")
terra::writeRaster(crop_cover_max_3035, "data/intermediate_output/vegetation/copernicus/crop_cover_max_3035.tif", overwrite = TRUE)



