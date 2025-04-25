#Sweeney, Harrison and Vander Linden 2024. 
#European forest cover during the Holocene reconstructed from pollen records
# ---------------------------------------------------------

#There are six seperate scripts associated with this research:
# 1. Map_setup.R 
# 2. Modern_pollen.R 
# 3. Model_import_data.R
# 4. Tree_model.R
# 5. Tree_recon.R (this script)
# 6. Figures.R

# ---------------------------------------------------------

############ Tree recon ############ 

# Using the model of modern tree cover and fossil pollen data, this script 
# reconstructs tree cover throught time. Comparison with other reconstructions
# is also performed


# ---------------------------------------------------------
# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Load in data and quantile mapping adjustment model
# 3. Import fossil pollen data and reshape
# 4. Reconstructions
# 5. (as 2,3,and partially 4 above, but where no shannon index included)

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
# install.packages("grDevices")
# install.packages("gridExtra")
# install.packages("ggthemes")
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





# 2. Load in data and quantile mapping adjustment model

#Input files
ap_tree_cover <- rio::import("data/intermediate_output/vegetation/ap_tree_cover.csv")
tree_model <- readRDS("data/intermediate_output/vegetation/tree_beta.rda") #modern tree cover model
boot_tree_model <- readRDS("data/intermediate_output/vegetation/boot_tree_beta.rda")
taxa_cat <- rio::import("data/input/taxa_cat.csv") %>% #import information regarding the taxa
  dplyr::arrange(taxon_name) 
taxa_cat_single <- taxa_cat %>% #for single clean_taxon_name
  dplyr::select(-taxon_name) %>% 
  dplyr::distinct()
loocv_predict <- rio::import( "data/intermediate_output/vegetation/loocv.csv")
boot_predict <- readRDS( "data/intermediate_output/vegetation/boot_dataset_pred_df.rda")

##Set map extents (as 1.Map_setup)
euro_extent_4258_ymin <- 34 
euro_extent_4258_ymax <- 73
euro_extent_4258_xmin <- -12
euro_extent_4258_xmax <- 45

#Load in necessary maps and rasters
euro_map_3035 <- sf::st_read("data/intermediate_output/vegetation/euro_map/3035/euro_map_3035.shp") #load map (1.Map_setup)
euro_map_3035_rast_50 <- raster::raster("data/intermediate_output/vegetation/euro_map/euro_map_3035_rast_50.tif") #(1.Map_setup)
euro_map_3035_rast_7 <- raster::raster("data/intermediate_output/vegetation/euro_map/euro_map_3035_rast_7.tif") #(1.Map_setup) - for use with Zanon data
euro_map_3035_rast_80 <- raster::raster("data/intermediate_output/vegetation/euro_map/euro_map_3035_rast_80.tif") #(1.Map_setup) - for use with Serge data

#Analyse relationships observed tree cover against LOOCV predictions 
fitted_tree_plot <- loocv_predict %>% 
  dplyr::select(ID_ENTITY, tree, prediction, E, latitude, longitude) %>% 
  dplyr::rename(observed = tree, predicted = prediction, difference = E) %>% #rename variables
  dplyr::mutate(difference = - difference) %>% #observation minus predictions
  dplyr::mutate(observed100 = observed*100, predicted100 = predicted*100) %>% #convert to precentage
  dplyr::mutate(observed_group = ggplot2::cut_width(observed100, width = 10, center = 5)) %>% 
  dplyr::arrange(observed_group)

rio::export(fitted_tree_plot, "data/intermediate_output/vegetation/fitted_tree_plot.csv")

fitted_tree_plot_nonaobs <- dplyr::filter(fitted_tree_plot, !is.na(observed)) #exclude records with  no observed cover
round(cor(fitted_tree_plot_nonaobs$observed, fitted_tree_plot_nonaobs$predicted),3) #correlation between obersevations and predictions
lm(fitted_tree_plot_nonaobs$predicted~fitted_tree_plot_nonaobs$observed)
round(max(fitted_tree_plot$predicted, na.rm = TRUE),2) #investigate maximum predicted values
round(max(fitted_tree_plot$observed, na.rm = TRUE),2) #investigate maximum observed values



#Empirical adjustment quantile mapping
tree_qmap_model_ssplin <- qmap::fitQmapSSPLIN(fitted_tree_plot_nonaobs$observed, fitted_tree_plot_nonaobs$predicted, wet.day = TRUE, qstep = 0.001) #quantile fit using smoothing spline
fitted_tree_plot_qmap <- fitted_tree_plot_nonaobs %>% 
  dplyr::mutate(refitted = qmap::doQmapSSPLIN(fitted_tree_plot_nonaobs$predicted, tree_qmap_model_ssplin)) %>% #re-fit values
  dplyr::mutate(difference2 = refitted-observed) %>%   #refitted minus observed
  dplyr::mutate(refitted100 = refitted*100) %>% 
  dplyr::arrange(observed_group)

rio::export(fitted_tree_plot_qmap, "data/intermediate_output/vegetation/fitted_tree_plot_qmap.csv")

refitted_tree_plot_nonaobs <- dplyr::filter(fitted_tree_plot_qmap, !is.na(observed))
round(cor(refitted_tree_plot_nonaobs$observed, refitted_tree_plot_nonaobs$refitted),3) #correlation between observations and adjusted predictions

refitted_model_3035 <- fitted_tree_plot_qmap %>%
  dplyr::select(longitude, latitude, difference2) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>%
  sf::st_transform(crs = 3035)

refitted_model_rast <- terra::rasterize(refitted_model_3035, euro_map_3035_rast_50, fun = mean)
refitted_model_df <- terra::as.data.frame(refitted_model_rast, xy = TRUE, na.rm = TRUE) %>% 
  dplyr::mutate(difference3 = ggplot2::cut_width(difference2, width = 0.4)) %>% 
  dplyr::mutate(Percentage_difference = dplyr::if_else(difference3 == "(-0.2,0.2]", "-20 to 20", difference3)) %>% 
  dplyr::mutate(Percentage_difference = dplyr::if_else(difference3 == "(-0.6,-0.2]", "-60 to -20", Percentage_difference)) %>%
  dplyr::mutate(Percentage_difference = dplyr::if_else(difference3 == "[-1,-0.6]", "-100 to -60", Percentage_difference)) %>%
  dplyr::mutate(Percentage_difference = dplyr::if_else(difference3 == "(0.2,0.6]", "20 to 60", Percentage_difference)) %>%
  dplyr::mutate(Percentage_difference = dplyr::if_else(difference3 == "(0.6,1]", "60 to 100", Percentage_difference)) 

refitted_model_df$Percentage_difference <- factor(refitted_model_df$Percentage_difference, levels = c("-100 to -60",
                                                                                                      "-60 to -20",
                                                                                                      "-20 to 20",
                                                                                                      "20 to 60",
                                                                                                      "60 to 100"))
rio::export(refitted_model_df, "data/intermediate_output/vegetation/refitted_model_df.csv")



#Analyse relationships observed tree cover against AP values
ap_tree_plot <- loocv_predict %>% 
  dplyr::select(ID_ENTITY, tree, ap_cover, latitude, longitude) %>% 
  dplyr::rename(observed = tree) %>% #rename variables
  dplyr::mutate(difference = observed - ap_cover) %>% 
  dplyr::mutate(difference = - difference) %>% # predictions minus observations
  dplyr::mutate(observed100 = observed*100, ap100 = ap_cover*100) %>% #convert to precentage
  dplyr::mutate(observed_group = ggplot2::cut_width(observed100, width = 10, center = 5)) %>% 
  dplyr::arrange(observed_group)

rio::export(ap_tree_plot, "data/intermediate_output/vegetation/ap_tree_plot.csv")

ap_tree_plot_nonaobs <- dplyr::filter(ap_tree_plot, !is.na(observed)) #exclude records with  no observed cover
round(cor(ap_tree_plot_nonaobs$observed, ap_tree_plot_nonaobs$ap_cover),3) #correlation between obersevations and predictions
round(max(ap_tree_plot$ap_cover, na.rm = TRUE),2) #investigate maximum ap values
round(max(ap_tree_plot$observed, na.rm = TRUE),2) #investigate maximum observed values


#Bootstraps
#Analyse relationships observed tree cover against LOOCV predictions 
boot_fitted_tree_plot <- boot_predict %>% 
  dplyr::mutate(E = tree - prediction) %>% 
  dplyr::select(ID_ENTITY, tree, prediction, E, latitude, longitude, boot) %>% 
  dplyr::rename(observed = tree, predicted = prediction, difference = E) %>% #rename variables
  dplyr::mutate(difference = - difference) %>% #observation minus predictions
  dplyr::mutate(observed100 = observed*100, predicted100 = predicted*100) %>% #convert to precentage
  dplyr::mutate(observed_group = ggplot2::cut_width(observed100, width = 10, center = 5)) %>% 
  dplyr::arrange(boot,observed_group)

rio::export(boot_fitted_tree_plot, "data/intermediate_output/vegetation/boot_fitted_tree_plot.csv")
boot_fitted_tree_plot_nonaobs <- dplyr::filter(boot_fitted_tree_plot, !is.na(observed)) #exclude records with  no observed cover
boot_fitted_tree_plot_nonaobs_ls <- boot_fitted_tree_plot_nonaobs %>% 
  dplyr::group_by(boot) %>% 
  dplyr::group_split()

#Empirical adjustment quantile mapping
boot_tree_qmap_model_ssplin <- lapply(boot_fitted_tree_plot_nonaobs_ls, function(x){
  qmap::fitQmapSSPLIN(x$observed, x$predicted, wet.day = TRUE, qstep = 0.001) #quantile fit using smoothing spline
})



# ---------------------------------------------------------




# 3. Import fossil pollen data and reshape
# ---------------------------------------------------------
# ---------------------------------------------------------

#Load special EPD data 
age_model <- rio::import("data/input/SPECIAL-EPD/age_model.csv") %>% 
  dplyr::filter(median != "unknown") %>% 
  dplyr::mutate(median = as.numeric(median))
dates <- rio::import("data/input/SPECIAL-EPD/dates.csv") %>% 
  dplyr::rename(depth = `depth (cm)`)
entity <- rio::import("data/input/SPECIAL-EPD/metadata.csv")
pollen_count <- rio::import("data/input/SPECIAL-EPD/pollen_counts_amalgamated.csv")
sample <- rio::import("data/input/SPECIAL-EPD/samples.csv") %>% 
  dplyr::rename(depth = `depth (cm)`)

#Rename and reshape data
multiple_entity_site_check <- sample %>% #limit entities to one per site
  dplyr::left_join(dplyr::select(age_model, ID_SAMPLE, median), by = "ID_SAMPLE") %>%
  dplyr::filter(median <= 12000) %>%
  dplyr::group_by(ID_ENTITY) %>%
  dplyr::summarise(number_per_entity = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(entity, ID_ENTITY, ID_SITE), by = "ID_ENTITY")  %>%
  dplyr::group_by(ID_SITE)  %>%
  dplyr::slice_max(number_per_entity, with_ties = FALSE) %>% #if same number then just one allowed
  dplyr::ungroup()

lq_multiple_entity_site_check <- sample %>% #limit entities to one per site
  dplyr::left_join(dplyr::select(age_model, ID_SAMPLE, median, UNCERT_25), by = "ID_SAMPLE")  %>%
  dplyr::filter(UNCERT_25 != "unknown") %>% #for quartile analysis of age model
  dplyr::mutate(lowerq = median + as.numeric(UNCERT_25)) %>%
  dplyr::filter(lowerq >= -70) %>%  #To prevent dates into the future
  dplyr::filter(lowerq <= 12000) %>%
  dplyr::group_by(ID_ENTITY) %>%
  dplyr::summarise(number_per_entity = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(entity, ID_ENTITY, ID_SITE), by = "ID_ENTITY")  %>%
  dplyr::group_by(ID_SITE) %>%
  dplyr::slice_max(number_per_entity, with_ties = FALSE) %>% #if same number then just one allowed
  dplyr::ungroup()

uq_multiple_entity_site_check <- sample %>% #limit entities to one per site
  dplyr::left_join(dplyr::select(age_model, ID_SAMPLE, median, UNCERT_75), by = "ID_SAMPLE")  %>%
  dplyr::filter(UNCERT_75 != "unknown") %>% #for quartile analysis of age model
  dplyr::mutate(upperq = median - as.numeric(UNCERT_75)) %>%
  dplyr::filter(upperq >= -70) %>%  #To prevent dates into the future
  dplyr::filter(upperq <= 12000) %>%
  dplyr::group_by(ID_ENTITY) %>%
  dplyr::summarise(number_per_entity = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(entity, ID_ENTITY, ID_SITE), by = "ID_ENTITY")  %>%
  dplyr::group_by(ID_SITE) %>%
  dplyr::slice_max(number_per_entity, with_ties = FALSE) %>% #if same number then just one allowed
  dplyr::ungroup()

entity_epd <- entity %>% 
  dplyr::mutate(handle = entity_name) %>% 
  dplyr::filter(ID_ENTITY %in% c(multiple_entity_site_check$ID_ENTITY)) 

lq_entity_epd <- entity %>% 
  dplyr::mutate(handle = entity_name) %>% 
  dplyr::filter(ID_ENTITY %in% c(lq_multiple_entity_site_check$ID_ENTITY)) 

uq_entity_epd <- entity %>% 
  dplyr::mutate(handle = entity_name) %>% 
  dplyr::filter(ID_ENTITY %in% c(uq_multiple_entity_site_check$ID_ENTITY)) 

sample_epd <- sample
age_model_epd <- age_model
pollen_count_epd <- pollen_count

pollen_sample_ages <- sample_epd %>% #Depth and date info regarding samples
  dplyr::left_join(age_model_epd, by = "ID_SAMPLE") %>% 
  dplyr::left_join(entity_epd, by = "ID_ENTITY") %>% 
  dplyr::filter(!site_type %in% c("archaeological site", 
                                  "cave",
                                  "coastal, estuarine",
                                  "coastal, lagoon",
                                  "glacial",
                                  "fluvial",
                                  "marine"))  

lq_pollen_sample_ages <- sample_epd %>% #Depth and date info regarding samples
  dplyr::left_join(age_model_epd, by = "ID_SAMPLE")  %>% 
  dplyr::filter(UNCERT_25 != "unknown") %>% 
  dplyr::mutate(lowerq = median + as.numeric(UNCERT_25)) %>% 
  dplyr::filter(lowerq >= -70) %>%  #To prevent dates into the future
  dplyr::left_join(lq_entity_epd, by = "ID_ENTITY") %>% 
  dplyr::filter(!site_type %in% c("archaeological site", 
                                  "cave",
                                  "coastal, estuarine",
                                  "coastal, lagoon",
                                  "glacial",
                                  "fluvial",
                                  "marine"))   

uq_pollen_sample_ages <- sample_epd %>% #Depth and date info regarding samples
  dplyr::left_join(age_model_epd, by = "ID_SAMPLE")  %>% 
  dplyr::filter(UNCERT_75 != "unknown") %>% 
  dplyr::mutate(upperq = median - as.numeric(UNCERT_75)) %>% 
  dplyr::filter(upperq >= -70) %>%  #To prevent dates into the future
  dplyr::left_join(uq_entity_epd, by = "ID_ENTITY") %>% 
  dplyr::filter(!site_type %in% c("archaeological site", 
                                  "cave",
                                  "coastal, estuarine",
                                  "coastal, lagoon",
                                  "glacial",
                                  "fluvial",
                                  "marine"))  

pollen_sample_age_nona <- pollen_sample_ages %>% #filter to records with a median date
  dplyr::filter(!is.na(median))
total_epd_records <- length(unique(pollen_sample_age_nona$handle))

lq_pollen_sample_age_nona <- lq_pollen_sample_ages %>% 
  dplyr::filter(!is.na(lowerq))
lq_total_epd_records <- length(unique(lq_pollen_sample_age_nona$handle))

uq_pollen_sample_age_nona <- uq_pollen_sample_ages %>% 
  dplyr::filter(!is.na(upperq))
uq_total_epd_records <- length(unique(uq_pollen_sample_age_nona$handle))


#Select data
epd_counts_taxon_am2 <- pollen_count_epd %>% 
  tidyr::pivot_longer(!ID_SAMPLE, names_to = "taxon_name", values_to = "count") %>% 
  tidyr::drop_na() %>% 
  dplyr::left_join(dplyr::select(pollen_sample_ages, ID_SAMPLE, depth, median, handle), by = "ID_SAMPLE") %>% #Ensure required details added
  dplyr::rename(quantity = count) %>% #Rename for later
  dplyr::mutate(Notes = 0) %>% #Add for later combining
  dplyr::select(ID_SAMPLE, handle, depth, median, taxon_name, quantity, Notes) %>% 
  dplyr::filter(!is.na(median)) #To ensure all count info has a date

epd_pollen_counts <- epd_counts_taxon_am2 %>% #All of the pollen data, with dates
  dplyr::filter(median <= 14000) %>% #drop older data
  dplyr::filter(!is.na(quantity)) %>%  #drop na quantities
  dplyr::filter(quantity >= 0) %>% #get rid of negative quantities
  dplyr::arrange(handle, depth, taxon_name) %>% 
  dplyr::left_join(dplyr::select(entity_epd, handle, latitude, longitude, elevation), by = "handle") %>%  #add in coordinate and elevation info
  dplyr::filter(dplyr::between(latitude, euro_extent_4258_ymin, euro_extent_4258_ymax)) %>% #Limit to European range
  dplyr::filter(dplyr::between(longitude, euro_extent_4258_xmin, euro_extent_4258_xmax))

lq_epd_counts_taxon_am2 <- pollen_count_epd %>% 
  tidyr::pivot_longer(!ID_SAMPLE, names_to = "taxon_name", values_to = "count") %>%  #CSVs
  tidyr::drop_na() %>% #CSVs
  dplyr::left_join(dplyr::select(lq_pollen_sample_ages, ID_SAMPLE, depth, lowerq, handle), by = "ID_SAMPLE") %>% #Ensure required details added
  dplyr::rename(quantity = count) %>% #Rename for later
  dplyr::mutate(Notes = 0) %>% #Add for later combining
  dplyr::select(ID_SAMPLE, handle, depth, lowerq, taxon_name, quantity, Notes) %>% 
  dplyr::filter(!is.na(lowerq)) #To ensure all count info has a date

lq_epd_pollen_counts <- lq_epd_counts_taxon_am2 %>% #All of the pollen data, with dates
  dplyr::filter(lowerq <= 14000) %>% #drop older data
  dplyr::filter(!is.na(quantity)) %>%  #drop na quantities
  dplyr::filter(quantity >= 0) %>% #get rid of negative quantities
  dplyr::arrange(handle, depth, taxon_name) %>% 
  dplyr::left_join(dplyr::select(lq_entity_epd, handle, latitude, longitude, elevation), by = "handle") %>%  #add in coordinate and elevation info
  dplyr::filter(dplyr::between(latitude, euro_extent_4258_ymin, euro_extent_4258_ymax)) %>% #Limit to European range
  dplyr::filter(dplyr::between(longitude, euro_extent_4258_xmin, euro_extent_4258_xmax))

uq_epd_counts_taxon_am2 <- pollen_count_epd %>% 
  tidyr::pivot_longer(!ID_SAMPLE, names_to = "taxon_name", values_to = "count") %>%  #CSVs
  tidyr::drop_na() %>% #CSVs
  dplyr::left_join(dplyr::select(uq_pollen_sample_ages, ID_SAMPLE, depth, upperq, handle), by = "ID_SAMPLE") %>% #Ensure required details added
  dplyr::rename(quantity = count) %>% #Rename for later
  dplyr::mutate(Notes = 0) %>% #Add for later combining
  dplyr::select(ID_SAMPLE, handle, depth, upperq, taxon_name, quantity, Notes) %>% 
  dplyr::filter(!is.na(upperq))   #To ensure all count info has a date

uq_epd_pollen_counts <- uq_epd_counts_taxon_am2 %>% #All of the pollen data, with dates
  dplyr::filter(upperq <= 14000) %>% #drop older data
  dplyr::filter(!is.na(quantity)) %>%  #drop na quantities
  dplyr::filter(quantity >= 0) %>% #get rid of negative quantities
  dplyr::arrange(handle, depth, taxon_name) %>% 
  dplyr::left_join(dplyr::select(uq_entity_epd, handle, latitude, longitude, elevation), by = "handle") %>%  #add in coordinate and elevation info
  dplyr::filter(dplyr::between(latitude, euro_extent_4258_ymin, euro_extent_4258_ymax)) %>% #Limit to European range
  dplyr::filter(dplyr::between(longitude, euro_extent_4258_xmin, euro_extent_4258_xmax))

#Add taxon cleaner data and filter
epd_pollen_counts_cleaner <- epd_pollen_counts %>% 
  dplyr::left_join(taxa_cat, by = c("taxon_name")) %>% 
  dplyr::rename(date_value = median)
lq_epd_pollen_counts_cleaner <- lq_epd_pollen_counts %>% 
  dplyr::left_join(taxa_cat, by = c("taxon_name")) %>%
  dplyr::rename(date_value = lowerq)
uq_epd_pollen_counts_cleaner <- uq_epd_pollen_counts %>% 
  dplyr::left_join(taxa_cat, by = c("taxon_name")) %>%
  dplyr::rename(date_value = upperq)

epd_pollen_counts_cleaner_ls <- list(epd_pollen_counts_cleaner,lq_epd_pollen_counts_cleaner,uq_epd_pollen_counts_cleaner)
names(epd_pollen_counts_cleaner_ls) <- c("median", "lq", "uq")

epd_pollen_counts_tps_ls <- parallel::mclapply(epd_pollen_counts_cleaner_ls, function(x){
  x %>% 
    dplyr::filter(terrestrial_pollen_sum == "yes") %>% #just those within TPS
    dplyr::select(-taxon_name) %>% 
    dplyr::distinct() %>% #remove potential duplicate rows caused by clean_taxa_name being the same for some taxa_name
    dplyr::mutate(europe = dplyr::if_else(is.na(europe), "none", europe)) %>% 
    dplyr::filter(europe != "not native to Europe") %>%  #remove non native species
    dplyr::select(-europe) %>% 
    dplyr::select(ID_SAMPLE, handle, depth, date_value, latitude, longitude, elevation, clean_taxon_name, quantity) %>%
    dplyr::mutate(clean_taxon_name = stringr::str_replace_all(clean_taxon_name, " type", "")) %>% #to stop type being a seperate species
    dplyr::group_by(ID_SAMPLE, handle, depth, date_value, latitude, longitude, elevation, clean_taxon_name) %>% #to sum quantities over same species
    dplyr::summarise(quantity = sum(quantity)) %>% 
    dplyr::ungroup() 
}, mc.cores = 6)


##Calculate indices
epd_pollen_wider_ls <- parallel::mclapply(epd_pollen_counts_tps_ls, function(x){
  x %>% 
    dplyr::arrange(ID_SAMPLE, clean_taxon_name) %>% 
    dplyr::select(ID_SAMPLE, clean_taxon_name, quantity) %>% 
    tidyr::pivot_wider(names_from = clean_taxon_name, values_from = quantity, values_fill = 0) 
},mc.cores = 6)

#Pollen percentages
epd_pollen_wider_per_ls <- parallel::mclapply(epd_pollen_wider_ls, function(x){
  x %>% 
    smpds::normalise_taxa(cols = 1)  #calculate percentages, ignoring first column
},mc.cores = 6)


#Shannon index
epd_shannon_ls <- parallel::mclapply(epd_pollen_wider_per_ls, function(x){
  x %>% 
    dplyr::select(-ID_SAMPLE) %>% 
    as.matrix() %>% 
    vegan::diversity(., index = "shannon") %>% 
    dplyr::as_tibble() %>%   
    dplyr::mutate(ID_SAMPLE = x$ID_SAMPLE, .before = 1) %>% 
    dplyr::rename(shannon = value)
}, mc.cores = 6)
epd_tree_shannon_ls <- parallel::mclapply(epd_pollen_wider_per_ls, function(x){
  x %>% 
    tidyr::pivot_longer(!c(ID_SAMPLE), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% #re-shape 
      dplyr::left_join(dplyr::select(taxa_cat_single, clean_taxon_name, ap_sp_hp), by = "clean_taxon_name") %>%
      dplyr::filter(ap_sp_hp == "AP") %>%
      dplyr::select(ID_SAMPLE, clean_taxon_name, percentage_cover) %>%
      tidyr::pivot_wider(names_from = clean_taxon_name, values_from = percentage_cover) %>%
      dplyr::select(-ID_SAMPLE) %>%
      as.matrix() %>%
      vegan::diversity(., index = "shannon") %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(ID_SAMPLE = x$ID_SAMPLE, .before = 1) %>%
      dplyr::rename(tree_shannon = value)
}, mc.cores = 6)

epd_tree_shannon_bin200 <- epd_tree_shannon_ls$median %>% 
  dplyr::left_join(dplyr::select(pollen_sample_ages, ID_SAMPLE, entity_name, latitude, longitude, median, site_type, elevation), by = "ID_SAMPLE")  %>% 
  dplyr::filter(site_type %in% c("lake", "terrestrial bog/mire/fen", "terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, bog/mire/fen", "terrestrial, marsh")) %>% 
  dplyr::filter(elevation < 1000) %>%  #Limit to records below 1000m)) 
  dplyr::select(entity_name, tree_shannon, median) %>% 
  dplyr::mutate(bin = ggplot2::cut_width(median, width = 200, center = 0, labels = F)) %>%    #place in 200 year bins
  dplyr::mutate(bin_age = (bin * 200) - 200) %>%
  dplyr::filter(bin_age <=14000)  %>% 
  dplyr::group_by(entity_name, bin_age) %>% 
  dplyr::summarise(tree_shannon = mean(tree_shannon)) %>%
  dplyr::ungroup() 

rio::export(epd_tree_shannon_bin200,"data/intermediate_output/vegetation/epd_tree_shannon_bin200.csv")

#Hills N2 filter
epd_hillsN2_ls <- parallel::mclapply(epd_pollen_wider_ls, function(x){
  x %>% 
    dplyr::select(-ID_SAMPLE) %>% 
    analogue::n2(., "sites") %>%
    dplyr::as_tibble() %>%
    dplyr::rename(HillsN2 = value) %>% 
    dplyr::mutate(ID_SAMPLE = x$ID_SAMPLE, .before = 1)
},mc.cores = 6)


#AP percentages
epd_pollen_counts_APinfo_ls <- parallel::mclapply(epd_pollen_wider_per_ls, function(x){
  x %>% 
    tidyr::pivot_longer(!c(ID_SAMPLE), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% 
    dplyr::left_join(taxa_cat_single, by = "clean_taxon_name") 
},mc.cores = 6)

epd_pollen_ap_per_ls <- parallel::mclapply(epd_pollen_counts_APinfo_ls, function(x){
  x %>% 
    dplyr::filter(ap_sp_hp == "AP") %>% #Just AP
    dplyr::group_by(ID_SAMPLE) %>% 
    dplyr::summarise(ap_cover = sum(percentage_cover)) %>% 
    dplyr::ungroup()
},mc.cores = 6)


#Needleshare 
epd_pollen_ap_needleshare_per_ls1 <- lapply(epd_pollen_counts_APinfo_ls,function(x){
  x %>%
    dplyr::filter(ap_sp_hp == "AP") %>%
    dplyr::mutate(needle = dplyr::if_else(is.na(ap_needle_broad), "unknown", ap_needle_broad)) %>%   #need to remove NA
    dplyr::filter(needle == "needle") %>%
    dplyr::group_by(ID_SAMPLE) %>%
    dplyr::summarise(ap_cover_needle = sum(percentage_cover)) %>%
    dplyr::ungroup()
})
epd_pollen_ap_needleshare_per_ls2 <- purrr::map2(epd_pollen_ap_needleshare_per_ls1, epd_pollen_ap_per_ls, dplyr::left_join)
epd_pollen_ap_needleshare_per_ls3 <- lapply(epd_pollen_ap_needleshare_per_ls2, function(x){
  x %>% 
    dplyr::mutate(needle_share = ap_cover_needle/ap_cover)  %>%
    dplyr::mutate(needle_share = tidyr::replace_na(needle_share, 0))
})

#Wind pollination
epd_pollen_ap_windshare_per_ls1 <- lapply(epd_pollen_counts_APinfo_ls,function(x){
  x %>%
    dplyr::filter(ap_sp_hp == "AP") %>%
    dplyr::mutate(wind = dplyr::if_else(is.na(pollination), "unknown", pollination)) %>%   #need to remove NA
    dplyr::filter(wind == "wind") %>%
    dplyr::group_by(ID_SAMPLE) %>%
    dplyr::summarise(ap_cover_wind = sum(percentage_cover)) %>%
    dplyr::ungroup()
})
epd_pollen_ap_windshare_per_ls2 <- purrr::map2(epd_pollen_ap_windshare_per_ls1, epd_pollen_ap_per_ls, dplyr::left_join)
epd_pollen_ap_wind_per_ls3 <- lapply(epd_pollen_ap_wind_per_ls2, function(x){
  x %>% 
    dplyr::mutate(wind_share = ap_cover_wind/ap_cover)  %>%
    dplyr::mutate(wind_share = tidyr::replace_na(wind_share, 0))
})

#Pinus
epd_pollen_ap_pinusshare_per_ls1 <- lapply(epd_pollen_counts_APinfo_ls,function(x){
  x %>%
    dplyr::filter(ap_sp_hp == "AP") %>%
    dplyr::filter(clean_taxon_name %in% c("Pinus","Pinus (diploxylon)","Pinus (haploxylon)")) %>%
    dplyr::group_by(ID_SAMPLE) %>%
    dplyr::summarise(ap_cover_pinus = sum(percentage_cover)) %>%
    dplyr::ungroup()
})
epd_pollen_ap_pinusshare_per_ls2 <- purrr::map2(epd_pollen_ap_pinusshare_per_ls1, epd_pollen_ap_per_ls, dplyr::left_join)
epd_pollen_ap_pinusshare_per_ls3 <- lapply(epd_pollen_ap_pinus_per_ls2, function(x){
  x %>% 
    dplyr::mutate(pinus_share = ap_cover_pinus/ap_cover)  %>%
    dplyr::mutate(pinus_share = tidyr::replace_na(pinus_share, 0))
})

#SP cover
epd_pollen_counts_SPinfo_ls <- parallel::mclapply(epd_pollen_wider_per_ls, function(x){
  x %>% 
    tidyr::pivot_longer(!c(ID_SAMPLE), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% 
    dplyr::left_join(taxa_cat_single, by = "clean_taxon_name") 
}, mc.cores = 6)

epd_pollen_sp_per_ls <- parallel::mclapply(epd_pollen_counts_SPinfo_ls, function(x){
  x %>% 
    dplyr::filter(ap_sp_hp == "SP") %>% #Just SP
    dplyr::group_by(ID_SAMPLE) %>% 
    dplyr::summarise(sp_cover = sum(percentage_cover)) %>% 
    dplyr::ungroup()
}, mc.cores = 6)


##Input for downcore
epd_pollen_counts_tps <- epd_pollen_counts_tps_ls$median
epd_shannon <- epd_shannon_ls$median
epd_tree_shannon <- epd_tree_shannon_ls$median
epd_hillsN2 <- epd_hillsN2_ls$median
epd_pollen_ap_per <- epd_pollen_ap_per_ls$median
epd_pollen_sp_per <- epd_pollen_sp_per_ls$median
epd_pollen_ap_needleshare_per <- epd_pollen_ap_needleshare_per_ls3$median
# epd_pollen_ap_windshare_per <- epd_pollen_ap_windshare_per_ls3$median
# epd_pollen_ap_pinusshare_per <- epd_pollen_ap_pinusshare_per_ls3$median

lq_epd_pollen_counts_tps <- epd_pollen_counts_tps_ls$lq
lq_epd_shannon <- epd_shannon_ls$lq
lq_epd_tree_shannon <- epd_tree_shannon_ls$lq
lq_epd_hillsN2 <- epd_hillsN2_ls$lq
lq_epd_pollen_ap_per <- epd_pollen_ap_per_ls$lq
lq_epd_pollen_sp_per <- epd_pollen_sp_per_ls$lq
lq_epd_pollen_ap_needleshare_per <- epd_pollen_ap_needleshare_per_ls3$lq
# lq_epd_pollen_ap_windshare_per <- epd_pollen_ap_windshare_per_ls3$lq
# lq_epd_pollen_ap_pinusshare_per <- epd_pollen_ap_pinusshare_per_ls3$lq

uq_epd_pollen_counts_tps <- epd_pollen_counts_tps_ls$uq
uq_epd_shannon <- epd_shannon_ls$uq
uq_epd_tree_shannon <- epd_tree_shannon_ls$uq
uq_epd_hillsN2 <- epd_hillsN2_ls$uq
uq_epd_pollen_ap_per <- epd_pollen_ap_per_ls$uq
uq_epd_pollen_sp_per <- epd_pollen_sp_per_ls$uq
uq_epd_pollen_ap_needleshare_per <- epd_pollen_ap_needleshare_per_ls3$uq
# uq_epd_pollen_ap_windshare_per <- epd_pollen_ap_windshare_per_ls3$uq
# uq_epd_pollen_ap_pinusshare_per <- epd_pollen_ap_pinusshare_per_ls3$uq

epd_downcore_input <- epd_pollen_counts_tps %>% 
  dplyr::select(ID_SAMPLE, elevation) %>%
  dplyr::distinct() %>% 
  dplyr::left_join(epd_shannon, by = "ID_SAMPLE") %>%
  dplyr::left_join(epd_tree_shannon, by = "ID_SAMPLE") %>%
  dplyr::left_join(epd_hillsN2, by = "ID_SAMPLE") %>% 
  dplyr::left_join(epd_pollen_ap_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(ap_cover = ap_cover / 100) %>% 
  dplyr::left_join(epd_pollen_sp_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(sp_cover = sp_cover / 100) %>% 
  dplyr::left_join(dplyr::select(epd_pollen_ap_needleshare_per, ID_SAMPLE, needle_share), by = "ID_SAMPLE") %>%
  # dplyr::left_join(dplyr::select(epd_pollen_ap_windshare_per, ID_SAMPLE, wind_share), by = "ID_SAMPLE") %>%
  # dplyr::left_join(dplyr::select(epd_pollen_ap_pinusshare_per, ID_SAMPLE, pinus_share), by = "ID_SAMPLE") %>%
  dplyr::filter(HillsN2 >= 2) %>% #Hills filter
  dplyr::select(-HillsN2) %>% 
  dplyr::left_join(dplyr::select(pollen_sample_ages, ID_SAMPLE, site_type), by = "ID_SAMPLE") %>% 
  dplyr::mutate(site_model = dplyr::if_else(site_type == "lake", 1, 0)) %>% 
  dplyr::select(-site_type)

rio::export(epd_downcore_input, "data/intermediate_output/vegetation/epd_downcore_input.csv")

lq_epd_downcore_input <- lq_epd_pollen_counts_tps %>% 
  dplyr::select(ID_SAMPLE, elevation) %>%
  dplyr::distinct() %>% 
  dplyr::left_join(lq_epd_shannon, by = "ID_SAMPLE") %>%
  dplyr::left_join(lq_epd_tree_shannon, by = "ID_SAMPLE") %>%
  dplyr::left_join(lq_epd_hillsN2, by = "ID_SAMPLE") %>% 
  dplyr::left_join(lq_epd_pollen_ap_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(ap_cover = ap_cover / 100) %>% 
  dplyr::left_join(lq_epd_pollen_sp_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(sp_cover = sp_cover / 100) %>% 
  dplyr::left_join(dplyr::select(lq_epd_pollen_ap_needleshare_per, ID_SAMPLE, needle_share), by = "ID_SAMPLE") %>% 
  # dplyr::left_join(dplyr::select(lq_epd_pollen_ap_windshare_per, ID_SAMPLE, wind_share), by = "ID_SAMPLE") %>%
  # dplyr::left_join(dplyr::select(lq_epd_pollen_ap_pinusshare_per, ID_SAMPLE, pinus_share), by = "ID_SAMPLE") %>%
  dplyr::filter(HillsN2 >= 2) %>% #Hills filter
  dplyr::select(-HillsN2) %>% 
  dplyr::left_join(dplyr::select(lq_pollen_sample_ages, ID_SAMPLE, site_type), by = "ID_SAMPLE") %>% 
  dplyr::mutate(site_model = dplyr::if_else(site_type == "lake", 1, 0)) %>% 
  dplyr::select(-site_type)

rio::export(lq_epd_downcore_input, "data/intermediate_output/vegetation/lq_epd_downcore_input.csv")

uq_epd_downcore_input <- uq_epd_pollen_counts_tps %>% 
  dplyr::select(ID_SAMPLE, elevation) %>%
  dplyr::distinct() %>% 
  dplyr::left_join(uq_epd_shannon, by = "ID_SAMPLE") %>%
  dplyr::left_join(uq_epd_tree_shannon, by = "ID_SAMPLE") %>%
  dplyr::left_join(uq_epd_hillsN2, by = "ID_SAMPLE") %>% 
  dplyr::left_join(uq_epd_pollen_ap_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(ap_cover = ap_cover / 100) %>% 
  dplyr::left_join(uq_epd_pollen_sp_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(sp_cover = sp_cover / 100) %>% 
  dplyr::left_join(dplyr::select(uq_epd_pollen_ap_needleshare_per, ID_SAMPLE, needle_share), by = "ID_SAMPLE") %>% 
  # dplyr::left_join(dplyr::select(uq_epd_pollen_ap_windshare_per, ID_SAMPLE, wind_share), by = "ID_SAMPLE") %>%
  # dplyr::left_join(dplyr::select(uq_epd_pollen_ap_pinusshare_per, ID_SAMPLE, pinus_share), by = "ID_SAMPLE") %>%
  dplyr::filter(HillsN2 >= 2) %>% #Hills filter
  dplyr::select(-HillsN2) %>% 
  dplyr::left_join(dplyr::select(uq_pollen_sample_ages, ID_SAMPLE, site_type), by = "ID_SAMPLE") %>% 
  dplyr::mutate(site_model = dplyr::if_else(site_type == "lake", 1, 0)) %>% 
  dplyr::select(-site_type)

rio::export(uq_epd_downcore_input, "data/intermediate_output/vegetation/uq_epd_downcore_input.csv")

epd_downcore_input_records <- epd_downcore_input %>% 
  dplyr::left_join(dplyr::select(pollen_sample_ages, ID_SAMPLE, ID_ENTITY, latitude, longitude, site_type), by = "ID_SAMPLE") %>% 
  dplyr::filter(site_type %in% c("lake", "terrestrial bog/mire/fen", "terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, bog/mire/fen", "terrestrial, marsh")) %>% 
  dplyr::filter(elevation < 1000) %>% #Limit to records below 1000m
  dplyr::select(ID_ENTITY, latitude, longitude) %>% 
  dplyr::distinct()

pollen_records2 <- rio::import("data/intermediate_output/vegetation/pollen_records2.csv")
pollen_records2_filt <- pollen_records2 %>% 
  dplyr::filter(ID_ENTITY %in% c(ap_tree_cover$ID_ENTITY))

epd_downcore_input_records_vect <- epd_downcore_input_records %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>% 
  sf::st_transform(crs = 3035) %>% 
  terra::vect()
epd_downcore_input_records_buf5k <- epd_downcore_input_records %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>% 
  sf::st_transform(crs = 3035) %>% 
  sf::st_buffer(5000)
epd_downcore_input_records_buf_sourcemedian <- epd_downcore_input_records %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>% 
  sf::st_transform(crs = 3035) %>% 
  sf::st_buffer(median(pollen_records2_filt$i0.75))

cop_tree_cover_masked_3035 <-  terra::rast("data/intermediate_output/vegetation/copernicus/cop_masked_tree_cover_3035_100m.tif")
epd_downcore_input_observedtree <- terra::extract(cop_tree_cover_masked_3035, epd_downcore_input_records_vect)
epd_downcore_input_observedtree_5k <- exactextractr::exact_extract(cop_tree_cover_masked_3035, epd_downcore_input_records_buf5k, 'mean', force_df = TRUE)
epd_downcore_input_observedtree_sourcemedian <- exactextractr::exact_extract(cop_tree_cover_masked_3035, epd_downcore_input_records_buf_sourcemedian, 'mean', force_df = TRUE)
median(epd_downcore_input_observedtree$mean, na.rm = TRUE)
median(epd_downcore_input_observedtree_5k$mean, na.rm = TRUE)
median(epd_downcore_input_observedtree_sourcemedian$mean, na.rm = TRUE)

# ---------------------------------------------------------




# 4. Reconstructions
# ---------------------------------------------------------
# ---------------------------------------------------------
## Modelled fit
epd_downcore_input_ls <- list(rio::import("data/intermediate_output/vegetation/epd_downcore_input.csv"),
                              rio::import("data/intermediate_output/vegetation/lq_epd_downcore_input.csv"),
                              rio::import("data/intermediate_output/vegetation/uq_epd_downcore_input.csv"))
names(epd_downcore_input_ls) <- c("median", "lq", "uq")
pred_EPD_tree_ls <- lapply(epd_downcore_input_ls, function(x){
  betareg::predict(tree_model, newdata = dplyr::select(x, ap_cover, elevation, needle_share, tree_shannon, site_model, sp_cover),type="response") #modelled
})

pred_EPD_tree_refit_ls <- lapply(pred_EPD_tree_ls, function(x){
  qmap::doQmapSSPLIN(as.numeric(x), tree_qmap_model_ssplin)
})

recon_tree <- cbind(epd_downcore_input_ls$median, pred_EPD_tree_refit_ls$median)  %>%
  dplyr::rename(tree_cover = dplyr::last_col()) %>%
  dplyr::mutate(tree_cover = dplyr::if_else(tree_cover>1, 1, tree_cover)) %>% #Because otherwise we may have re-fitted values greater than 100%
  dplyr::mutate(tree_cover = 100*tree_cover) %>% 
  dplyr::left_join(dplyr::select(pollen_sample_ages, ID_SAMPLE, entity_name, latitude, longitude, median, site_type), by = "ID_SAMPLE")  %>% 
  dplyr::filter(site_type %in% c("lake", "terrestrial bog/mire/fen", "terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, bog/mire/fen", "terrestrial, marsh")) %>% 
  dplyr::filter(elevation < 1000) #Limit to records below 1000m
rio::export(recon_tree, "data/intermediate_output/vegetation/recon_tree.csv")

lq_recon_tree <- cbind(epd_downcore_input_ls$lq, pred_EPD_tree_refit_ls$lq)  %>%
  dplyr::rename(tree_cover = dplyr::last_col()) %>%
  dplyr::mutate(tree_cover = dplyr::if_else(tree_cover>1, 1, tree_cover)) %>% #Because otherwise we may have re-fitted values greater than 100%
  dplyr::mutate(tree_cover = 100*tree_cover) %>% 
  dplyr::left_join(dplyr::select(lq_pollen_sample_ages, ID_SAMPLE, entity_name, latitude, longitude, lowerq, site_type), by = "ID_SAMPLE")  %>% 
  dplyr::filter(site_type %in% c("lake", "terrestrial bog/mire/fen", "terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, bog/mire/fen", "terrestrial, marsh")) %>% 
  dplyr::filter(elevation < 1000) #Limit to records below 1000m
rio::export(lq_recon_tree, "data/intermediate_output/vegetation/lq_recon_tree.csv")

uq_recon_tree <- cbind(epd_downcore_input_ls$uq, pred_EPD_tree_refit_ls$uq)  %>%
  dplyr::rename(tree_cover = dplyr::last_col()) %>%
  dplyr::mutate(tree_cover = dplyr::if_else(tree_cover>1, 1, tree_cover)) %>% #Because otherwise we may have re-fitted values greater than 100%
  dplyr::mutate(tree_cover = 100*tree_cover) %>% 
  dplyr::left_join(dplyr::select(uq_pollen_sample_ages, ID_SAMPLE, entity_name, latitude, longitude, upperq, site_type), by = "ID_SAMPLE")  %>% 
  dplyr::filter(site_type %in% c("lake", "terrestrial bog/mire/fen", "terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, bog/mire/fen", "terrestrial, marsh")) %>% 
  dplyr::filter(elevation < 1000) #Limit to records below 1000m
rio::export(uq_recon_tree, "data/intermediate_output/vegetation/uq_recon_tree.csv")

ap_cover_tree <- recon_tree %>% 
  dplyr::mutate(tree_cover = ap_cover*100) #convert to ap cover subsequent analysis
rio::export(ap_cover_tree, "data/intermediate_output/vegetation/ap_cover_tree.csv")

#Bin data into 200yrs bins
recon_tree <- rio::import("data/intermediate_output/vegetation/recon_tree.csv")
lq_recon_tree <- rio::import("data/intermediate_output/vegetation/lq_recon_tree.csv")
uq_recon_tree <- rio::import("data/intermediate_output/vegetation/uq_recon_tree.csv")
ap_cover_tree <- rio::import("data/intermediate_output/vegetation/ap_cover_tree.csv")

recon_tree_ls <- list(recon_tree, lq_recon_tree, uq_recon_tree, ap_cover_tree)
names(recon_tree_ls) <- c("median", "lq", "uq", "ap")

recon_tree_binned_200_ls <- lapply(recon_tree_ls, function(x){
  x %>% 
    dplyr::rename(date_value = 13) %>% 
    dplyr::mutate(bin = ggplot2::cut_width(date_value, width = 200, center = 0, labels = F)) %>%    #place in 200 year bins
    dplyr::mutate(bin_age = (bin * 200) - 200) %>%
    dplyr::filter(bin_age <=14000)  %>% 
    dplyr::group_by(entity_name, bin_age) %>% 
    dplyr::summarise(tree_cover = mean(tree_cover), number = n()) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(dplyr::select(pollen_sample_ages, entity_name, latitude, longitude), by = "entity_name") %>% 
    dplyr::distinct()
})

bin_200_mean_recon_tree_ls <- lapply(recon_tree_binned_200_ls, function(x){ #Calculate mean through time
    x %>% 
      dplyr::group_by(bin_age) %>% 
      dplyr::summarise(tree_cover_mean = mean(tree_cover), tree_cover_sd = sd(tree_cover)) %>% 
      dplyr::ungroup()
}) 

bin_200_median_recon_tree_ls <- lapply(recon_tree_binned_200_ls, function(x){ #Calculate median through time
  x %>% 
    dplyr::group_by(bin_age) %>% 
    dplyr::summarise(tree_cover_median = median(tree_cover), tree_cover_lower = quantile(tree_cover,probs = c(0.25)),tree_cover_higher = quantile(tree_cover,probs = c(0.75))) %>% 
    dplyr::ungroup()
}) 

bin_200_max_recon_tree_ls <- lapply(recon_tree_binned_200_ls, function(x){ #Calculate max through time
  x %>% 
    dplyr::group_by(bin_age) %>% 
    dplyr::summarise(tree_cover_max = max(tree_cover)) %>% 
    dplyr::ungroup()
}) 

bin_200_n_recon_tree_ls <- lapply(recon_tree_binned_200_ls, function(x){ #Calculate number through time
  x %>% 
    dplyr::select(-tree_cover, -number) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(bin_age) %>% 
    dplyr::summarise(number_records = n()) %>% 
    dplyr::ungroup()
}) 

rio::export(recon_tree_binned_200_ls$median, "data/intermediate_output/vegetation/recon_tree_binned_200.csv")
rio::export(bin_200_mean_recon_tree_ls$median, "data/intermediate_output/vegetation/bin_200_mean_recon_tree.csv")
rio::export(bin_200_median_recon_tree_ls$median, "data/intermediate_output/vegetation/bin_200_median_recon_tree.csv")
rio::export(bin_200_max_recon_tree_ls$median, "data/intermediate_output/vegetation/bin_200_max_recon_tree.csv")
rio::export(bin_200_n_recon_tree_ls$median, "data/intermediate_output/vegetation/bin_200_n_recon_tree.csv")

rio::export(recon_tree_binned_200_ls$lq, "data/intermediate_output/vegetation/lq_recon_tree_binned_200.csv")
rio::export(bin_200_mean_recon_tree_ls$lq, "data/intermediate_output/vegetation/lq_bin_200_mean_recon_tree.csv")
rio::export(bin_200_median_recon_tree_ls$lq, "data/intermediate_output/vegetation/lq_bin_200_median_recon_tree.csv")
rio::export(bin_200_max_recon_tree_ls$lq, "data/intermediate_output/vegetation/lq_bin_200_max_recon_tree.csv")
rio::export(bin_200_n_recon_tree_ls$lq, "data/intermediate_output/vegetation/lq_bin_200_n_recon_tree.csv")

rio::export(recon_tree_binned_200_ls$uq, "data/intermediate_output/vegetation/uq_recon_tree_binned_200.csv")
rio::export(bin_200_mean_recon_tree_ls$uq, "data/intermediate_output/vegetation/uq_bin_200_mean_recon_tree.csv")
rio::export(bin_200_median_recon_tree_ls$uq, "data/intermediate_output/vegetation/uq_bin_200_median_recon_tree.csv")
rio::export(bin_200_max_recon_tree_ls$uq, "data/intermediate_output/vegetation/uq_bin_200_max_recon_tree.csv")
rio::export(bin_200_n_recon_tree_ls$uq, "data/intermediate_output/vegetation/uq_bin_200_n_recon_tree.csv")

rio::export(recon_tree_binned_200_ls$ap, "data/intermediate_output/vegetation/ap_tree_binned_200.csv")
rio::export(bin_200_mean_recon_tree_ls$ap, "data/intermediate_output/vegetation/bin_200_mean_ap_tree.csv")
rio::export(bin_200_median_recon_tree_ls$ap, "data/intermediate_output/vegetation/bin_200_median_ap_tree.csv")
rio::export(bin_200_max_recon_tree_ls$ap, "data/intermediate_output/vegetation/bin_200_max_ap_tree.csv")
rio::export(bin_200_n_recon_tree_ls$ap, "data/intermediate_output/vegetation/bin_200_n_ap_tree.csv")

recon_shannon_binned_200_ls <- lapply(recon_tree_ls, function(x){
  x %>% 
    dplyr::rename(date_value = 13) %>% 
    dplyr::mutate(bin = ggplot2::cut_width(date_value, width = 200, center = 0, labels = F)) %>%    #place in 200 year bins
    dplyr::mutate(bin_age = (bin * 200) - 200) %>%
    dplyr::filter(bin_age <=14000)  %>% 
    dplyr::group_by(entity_name, bin_age) %>% 
    dplyr::summarise(shannon = mean(shannon)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(dplyr::select(pollen_sample_ages, entity_name, latitude, longitude), by = "entity_name") %>% 
    dplyr::distinct()
})

ggplot2::ggplot(data = recon_shannon_binned_200_ls$median, mapping = aes(x = bin_age, y = shannon))+
  geom_point()+
  geom_smooth()+
  theme_bw()+
  scale_x_reverse(limit=c(12000,0))

bin_200_median_recon_shannon_ls <- lapply(recon_shannon_binned_200_ls, function(x){ #Calculate median through time
  x %>% 
    dplyr::group_by(bin_age) %>% 
    dplyr::summarise(shannon_median = median(shannon), shannon_lower = quantile(shannon,probs = c(0.25)),shannon_higher = quantile(shannon,probs = c(0.75))) %>% 
    dplyr::ungroup()
}) 

ggplot2::ggplot(data = bin_200_median_recon_shannon_ls$median, mapping = aes(x = bin_age, y = shannon_median))+
  geom_point()+
  geom_smooth()+
  theme_bw()+
  scale_x_reverse(limit=c(12000,0))

##Reconstruction analysis
recon_tree <- rio::import("data/intermediate_output/vegetation/recon_tree.csv")
recon_tree_binned_200 <- rio::import("data/intermediate_output/vegetation/recon_tree_binned_200.csv")
bin_200_mean_recon_tree <- rio::import("data/intermediate_output/vegetation/bin_200_mean_recon_tree.csv")
bin_200_median_recon_tree <- rio::import("data/intermediate_output/vegetation/bin_200_median_recon_tree.csv")
bin_200_max_recon_tree <- rio::import("data/intermediate_output/vegetation/bin_200_max_recon_tree.csv")
bin_200_n_recon_tree <- rio::import("data/intermediate_output/vegetation/bin_200_n_recon_tree.csv")

binned_200yr_reconstructions <- recon_tree_binned_200 %>% 
  dplyr::rename(bin_centre = bin_age, number_samples = number) %>% 
  dplyr::select(entity_name, longitude, latitude, bin_centre, number_samples, tree_cover)
rio::export(binned_200yr_reconstructions, "figs/binned_200yr_reconstructions.csv")

#Bootstraps medians
nrecords <- length(unique(recon_tree_binned_200$entity_name)) #Number of records
recon_tree_binned_200_records <- recon_tree_binned_200 %>% #Add row number of each entity
  dplyr::select(entity_name) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(entity_row = dplyr::row_number())
recon_tree_binned_200_num <- recon_tree_binned_200 %>% 
  dplyr::left_join(recon_tree_binned_200_records, by = "entity_name")

boot_ls <- list()
set.seed(42)
nreps <- 1000 #number of reps
for (i in 1:nreps){ #generate list of boostrapped resampled dfs
  print(i)
  record_boot <- sample(seq(1:nrecords), nrecords, replace = TRUE) #select records, with replacement
  bin_200_recon_tree_boot <- dplyr::tibble(entity_row = record_boot) %>% 
    dplyr::left_join(recon_tree_binned_200_num, by = "entity_row", relationship = "many-to-many") %>% #add info by entity
    dplyr::group_by(bin_age) %>% 
    dplyr::summarise(tree_cover_mean = mean(tree_cover), tree_cover_median = median(tree_cover),tree_cover_max = max(tree_cover) ) %>% #calculate bin mean, median and max
    dplyr::ungroup() 
  boot_ls[[i]] <- bin_200_recon_tree_boot #add to list
  
}
bin_200_recon_tree_boots <- dplyr::bind_rows(boot_ls, .id = "column_label") #reduce to single df
rio::export(bin_200_recon_tree_boots, "data/intermediate_output/vegetation/bin_200_recon_tree_boots.csv")

bin_200_recon_tree_boots_5_95 <- bin_200_recon_tree_boots %>% #Calculate 95% CIs
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(tree_cover_lower_mean = quantile(tree_cover_mean,probs = c(0.025)),
                   tree_cover_higher_mean = quantile(tree_cover_mean,probs = c(0.975)),
                   tree_cover_lower_median = quantile(tree_cover_median,probs = c(0.025)),
                   tree_cover_higher_median = quantile(tree_cover_median,probs = c(0.975)),
                   tree_cover_lower_max = quantile(tree_cover_max,probs = c(0.025)),
                   tree_cover_higher_max = quantile(tree_cover_max,probs = c(0.975)),
                   ) %>% 
  dplyr::ungroup()  %>% 
  tidyr::pivot_longer(!bin_age, names_to = "lower_upper", values_to = "tree_cover_5_95") %>%  
  dplyr::mutate(column_label = dplyr::if_else(stringr::str_detect(lower_upper, "tree_cover_lower"), nreps+2, nreps+3)) #need to add for plotting
rio::export(bin_200_recon_tree_boots_5_95, "data/intermediate_output/vegetation/bin_200_recon_tree_boots_5_95.csv")

bin_200_mean_recon_tree_plot <- bin_200_mean_recon_tree %>% 
  dplyr::mutate(column_label = nreps+1)
rio::export(bin_200_mean_recon_tree_plot, "data/intermediate_output/vegetation/bin_200_mean_recon_tree_plot.csv")

bin_200_median_recon_tree_plot <- bin_200_median_recon_tree %>% 
  dplyr::mutate(column_label = nreps+1) #need to include within plot
rio::export(bin_200_median_recon_tree_plot, "data/intermediate_output/vegetation/bin_200_median_recon_tree_plot.csv")

bin_200_max_recon_tree_plot <- bin_200_max_recon_tree %>% 
  dplyr::mutate(column_label = nreps+1)
rio::export(bin_200_max_recon_tree_plot, "data/intermediate_output/vegetation/bin_200_max_recon_tree_plot.csv")

bin_200_recon_tree_boots_5_95_mean <- bin_200_recon_tree_boots_5_95 %>% 
  dplyr::filter(stringr::str_detect(lower_upper, "mean"))
rio::export(bin_200_recon_tree_boots_5_95_mean, "data/intermediate_output/vegetation/bin_200_recon_tree_boots_5_95_mean.csv")

bin_200_recon_tree_boots_5_95_median <- bin_200_recon_tree_boots_5_95 %>% 
  dplyr::filter(stringr::str_detect(lower_upper, "median"))
rio::export(bin_200_recon_tree_boots_5_95_median, "data/intermediate_output/vegetation/bin_200_recon_tree_boots_5_95_median.csv")

bin_200_recon_tree_boots_5_95_max <- bin_200_recon_tree_boots_5_95 %>% 
  dplyr::filter(stringr::str_detect(lower_upper, "max"))
rio::export(bin_200_recon_tree_boots_5_95_max, "data/intermediate_output/vegetation/bin_200_recon_tree_boots_5_95_max.csv")



#Locfit smoothing of tree cover through time
f_locfit_tree_average <- function(d){
  tree_loc <- d %>% 
  dplyr::select(bin_age, 2) %>% 
  dplyr::rename(tree_cover = 2)
  x <- as.vector(tree_loc$bin_age)
  y <- as.vector(tree_loc$tree_cover)
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
  loc01 <- locfit::locfit(y ~ locfit::lp(x, deg=1, h=500), maxk=800, family="qrgauss")
  pred01 <- predict(loc01, newdata=tree_loc$bin_age, se.fit=TRUE)
  loc02 <- locfit::locfit(y ~ locfit::lp(x, deg=1, h=1000), maxk=800, family="qrgauss")
  pred02 <- predict(loc02, newdata=tree_loc$bin_age, se.fit=TRUE)
  loc03 <- locfit::locfit(y ~ locfit::lp(x, deg=1, h=2000), maxk=800, family="qrgauss")
  pred03 <- predict(loc03, newdata=tree_loc$bin_age, se.fit=TRUE)
  locfit <- data.frame(tree_loc$bin_age, tree_loc$tree_cover, pred01$fit, pred02$fit, pred03$fit)
  name_average <- stringr::str_extract(deparse(substitute(d)), "(?<=_)[^_]+(?=_[^_]*_[^_]*$)")
  colnames(locfit) <- c("bin_age", "None", "500-year", "1000-year", "2000-year")
  locfit_long <- locfit %>% 
    tidyr::pivot_longer(!bin_age, names_to = "half_width", values_to = paste0(name_average)) %>% 
    dplyr::mutate(half_width = as.factor(half_width))
  assign(paste0(deparse(substitute(d)),"_loc"), locfit_long, envir = parent.frame())
}
f_locfit_tree_average(bin_200_median_recon_tree)
f_locfit_tree_average(bin_200_mean_recon_tree)

rio::export(bin_200_median_recon_tree_loc, "data/intermediate_output/vegetation/bin_200_median_recon_tree_loc.csv")
rio::export(bin_200_mean_recon_tree_loc, "data/intermediate_output/vegetation/bin_200_mean_recon_tree_loc.csv")




#Bioregions through time
europe_bioregions <- terra::rast("data/input/EEA_bioregions/EEA_bioregions.tif")
recon_buffer_5000 <- recon_tree %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>% 
  sf::st_transform(crs = 3035) %>% 
  sf::st_buffer(5000)

epd_bioregions_extracted <- exactextractr::exact_extract(europe_bioregions, recon_buffer_5000, "mode", force_df = TRUE) %>%   
  dplyr::mutate(entity_name = recon_tree$entity_name)  %>% 
  dplyr::mutate(latitude = recon_tree$latitude) %>% 
  dplyr::mutate(longitude = recon_tree$longitude) %>% 
  dplyr::rename(bioregion = mode) %>%  #buffer
  dplyr::mutate(bioregion = dplyr::if_else(is.na(bioregion), 9, bioregion))

epd_bioregion_names <- tibble(bioregion = sort(unique(epd_bioregions_extracted$bioregion)), #Order is important!
                          bioregion_name = c("Alpine", #0
                                         "Anatolian", #1
                                         "Arctic", #2
                                         "Atlantic", #3
                                         "BlackSea", #4
                                         "Boreal", #5
                                         "Continental", #6
                                         # "Macaronesia", #7
                                         "Mediterranean", #8
                                         "Outside", #9
                                         "Pannonian", #10
                                         "Steppic"))#11

epd_bioregion_map <- epd_bioregions_extracted %>% 
  dplyr::distinct() %>% 
  dplyr::left_join(epd_bioregion_names, by = "bioregion") %>% 
  tidyr::drop_na() %>% 
  dplyr::filter(bioregion_name != "Outside") %>% #For map, leave these records out. Due to projection they are off the map
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>% 
  sf::st_transform(crs = 3035)


#200 year bins
bioregion_recon_tree_binned_200 <- recon_tree_binned_200 %>% 
  dplyr::left_join(dplyr::select(epd_bioregions_extracted, entity_name, bioregion), by = "entity_name") %>%   
  dplyr::left_join(dplyr::select(epd_bioregion_names, bioregion, bioregion_name), by = "bioregion") %>% 
  dplyr::mutate(bioregion_name = dplyr::if_else(bioregion_name == "Outside" | is.na(bioregion_name), "Mediterranean", bioregion_name)) %>%  #sites that were in Northern Africa and Malta, Cyprus, converted to Med
  dplyr::distinct()
rio::export(bioregion_recon_tree_binned_200, "data/intermediate_output/vegetation/bioregion_recon_tree_binned_200.csv")
  
entity_bioregions <-  bioregion_recon_tree_binned_200 %>% 
  dplyr::select(entity_name, bioregion_name) %>% 
  dplyr::distinct()
rio::export(entity_bioregions, "data/intermediate_output/vegetation/entity_bioregions.csv")

bioregion_bin_200_mean_recon_tree <- bioregion_recon_tree_binned_200 %>% 
  dplyr::group_by(bin_age, bioregion_name) %>% 
  dplyr::summarise(tree_cover_mean = mean(tree_cover), tree_cover_sd = sd(tree_cover)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!is.na(bioregion_name))

bioregion_bin_200_median_recon_tree <- bioregion_recon_tree_binned_200 %>% 
  dplyr::group_by(bin_age, bioregion_name) %>% 
  dplyr::summarise(tree_cover_median = median(tree_cover), tree_cover_low = quantile(tree_cover, probs = 0.25), tree_cover_high = quantile(tree_cover, probs = 0.75)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!is.na(bioregion_name))

bioregion_bin_200_n_recon_tree <- bioregion_recon_tree_binned_200 %>% 
  dplyr::select(-tree_cover) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(bin_age, bioregion_name) %>% 
  dplyr::summarise(number_records = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!is.na(bioregion_name))

rio::export(entity_bioregions, "data/intermediate_output/vegetation/entity_bioregions.csv")
rio::export(bioregion_bin_200_mean_recon_tree, "data/intermediate_output/vegetation/bioregion_bin_200_mean_recon_tree.csv")
rio::export(bioregion_bin_200_median_recon_tree, "data/intermediate_output/vegetation/bioregion_bin_200_median_recon_tree.csv")
rio::export(bioregion_bin_200_n_recon_tree, "data/intermediate_output/vegetation/bioregion_bin_200_n_recon_tree.csv")

#Locfit smoothing of bioregion averages
bioregion_bin_200_median_recon_tree_ls <- bioregion_bin_200_median_recon_tree %>% 
  tidytable::group_split(bioregion_name, .keep = TRUE, .named = TRUE)
bioregion_bin_200_mean_recon_tree_ls <- bioregion_bin_200_mean_recon_tree %>% 
  tidytable::group_split(bioregion_name, .keep = TRUE, .named = TRUE)

f_locfit_tree_average_bioregion <- function(d){
  tree_loc <- d %>% 
    dplyr::select(bin_age, 2, 3) %>% 
    dplyr::rename(tree_cover = 3)
  x <- as.vector(tree_loc$bin_age)
  y <- as.vector(tree_loc$tree_cover)
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
  loc01 <- locfit::locfit(y ~ locfit::lp(x, deg=1, h=1000), maxk=800, family="qrgauss")
  pred01 <- predict(loc01, newdata=tree_loc$bin_age, se.fit=TRUE)
  locfit <- data.frame(tree_loc$bin_age, tree_loc$tree_cover, pred01$fit)
  colnames(locfit) <- c("bin_age", "None", "1000-year")
  locfit_long <- locfit %>% 
    tidyr::pivot_longer(!bin_age, names_to = "half_width", values_to = "average") %>% 
    dplyr::mutate(half_width = as.factor(half_width))
  return(locfit_long)
}
bioregion_bin_200_median_recon_tree_loc <- lapply(bioregion_bin_200_median_recon_tree_ls, f_locfit_tree_average_bioregion) %>% 
  dplyr::bind_rows(., .id = "bioregion_name")
bioregion_bin_200_mean_recon_tree_loc <- lapply(bioregion_bin_200_mean_recon_tree_ls, f_locfit_tree_average_bioregion) %>% 
  dplyr::bind_rows(., .id = "bioregion_name")




#Reduced grouping
bioregion_bin_200_mean_recon_tree_gr4 <- bioregion_bin_200_mean_recon_tree %>% dplyr::filter(bioregion_name %in% c("Atlantic", "Boreal", "Continental", "Mediterranean"))
bioregion_bin_200_median_recon_tree_gr4 <- bioregion_bin_200_median_recon_tree %>% dplyr::filter(bioregion_name %in% c("Atlantic", "Boreal", "Continental", "Mediterranean"))
bioregion_bin_200_median_recon_tree_gr4_all <- bioregion_recon_tree_binned_200 %>% dplyr::filter(bioregion_name %in% c("Atlantic", "Boreal", "Continental", "Mediterranean")) %>% 
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(tree_cover_median = median(tree_cover), tree_cover_low = quantile(tree_cover, probs = 0.25), tree_cover_high = quantile(tree_cover, probs = 0.75)) %>% 
  dplyr::ungroup()
bioregion_bin_200_n_recon_tree_gr4 <- bioregion_bin_200_n_recon_tree %>% dplyr::filter(bioregion_name %in% c("Atlantic", "Boreal", "Continental", "Mediterranean"))

epd_bioregion_map_gr4 <- bioregion_recon_tree_binned_200 %>% 
  dplyr::select(entity_name, latitude, longitude, bioregion_name) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(bioregion_name %in% c("Atlantic", "Boreal", "Continental", "Mediterranean"))  

rio::export(epd_bioregion_map_gr4, "data/intermediate_output/vegetation/epd_bioregion_map_gr4.csv")

bioregion_bin_200_median_recon_tree_loc_gr4 <- bioregion_bin_200_median_recon_tree_loc %>% dplyr::filter(bioregion_name %in% c("Atlantic", "Boreal", "Continental", "Mediterranean"))
rio::export(bioregion_bin_200_median_recon_tree_loc_gr4, "data/intermediate_output/vegetation/bioregion_bin_200_median_recon_tree_loc_gr4.csv")
bioregion_bin_200_mean_recon_tree_loc_gr4 <- bioregion_bin_200_mean_recon_tree_loc %>% dplyr::filter(bioregion_name %in% c("Atlantic", "Boreal", "Continental", "Mediterranean"))
rio::export(bioregion_bin_200_mean_recon_tree_loc_gr4, "data/intermediate_output/vegetation/bioregion_bin_200_mean_recon_tree_loc_gr4.csv")



#By bioregion, observed tree cover values
epd_downcore_input_observedtree_5k_bioregion <- epd_downcore_input_observedtree_5k %>% 
  dplyr::mutate(ID_ENTITY = epd_downcore_input_records_buf5k$ID_ENTITY) %>% 
  dplyr::left_join(dplyr::select(pollen_sample_ages, ID_ENTITY, handle), by = "ID_ENTITY") %>% 
  dplyr::rename(entity_name = handle) %>% 
  dplyr::left_join(dplyr::select(epd_bioregion_map, entity_name, bioregion_name), by = "entity_name") %>% 
  dplyr::distinct()
epd_downcore_input_observedtree_5k_bioregion_median <- epd_downcore_input_observedtree_5k_bioregion %>% 
  tidyr::drop_na() %>% 
  dplyr::group_by(bioregion_name) %>% 
  dplyr::summarise(median_tree = median(mean, na.rm = TRUE))
  

##AP values
ap_tree_binned_200 <- rio::import("data/intermediate_output/vegetation/ap_tree_binned_200.csv")
bin_200_mean_ap_tree <- rio::import("data/intermediate_output/vegetation/bin_200_mean_ap_tree.csv")
bin_200_median_ap_tree <- rio::import("data/intermediate_output/vegetation/bin_200_median_ap_tree.csv")
bin_200_max_ap_tree <- rio::import("data/intermediate_output/vegetation/bin_200_max_ap_tree.csv")
bin_200_n_ap_tree <- rio::import("data/intermediate_output/vegetation/bin_200_n_ap_tree.csv")


binned_200yr_ap <- ap_tree_binned_200 %>% 
  dplyr::rename(bin_centre = bin_age, number_samples = number) %>% 
  dplyr::select(entity_name, longitude, latitude, bin_centre, number_samples, tree_cover)
rio::export(binned_200yr_ap, "figs/binned_200yr_ap.csv")

#Bootstraps medians
nrecords <- length(unique(ap_tree_binned_200$entity_name)) #Number of records
ap_tree_binned_200_records <- ap_tree_binned_200 %>% #Add row number of each entity
  dplyr::select(entity_name) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(entity_row = dplyr::row_number())
ap_tree_binned_200_num <- ap_tree_binned_200 %>% 
  dplyr::left_join(ap_tree_binned_200_records, by = "entity_name")

boot_ls <- list()
set.seed(42)
nreps <- 1000 #number of reps
for (i in 1:nreps){ #generate list of boostrapped resampled dfs
  print(i)
  record_boot <- sample(seq(1:nrecords), nrecords, replace = TRUE) #select records, with replacement
  bin_200_ap_tree_boot <- dplyr::tibble(entity_row = record_boot) %>% 
    dplyr::left_join(ap_tree_binned_200_num, by = "entity_row", relationship = "many-to-many") %>% #add info by entity
    dplyr::group_by(bin_age) %>% 
    dplyr::summarise(tree_cover_mean = mean(tree_cover), tree_cover_median = median(tree_cover),tree_cover_max = max(tree_cover) ) %>% #calculate bin mean, median and max
    dplyr::ungroup() 
  boot_ls[[i]] <- bin_200_ap_tree_boot #add to list
  
}
bin_200_ap_tree_boots <- dplyr::bind_rows(boot_ls, .id = "column_label") #reduce to single df
rio::export(bin_200_ap_tree_boots, "data/intermediate_output/vegetation/bin_200_ap_tree_boots.csv")

bin_200_ap_tree_boots_5_95 <- bin_200_ap_tree_boots %>% #Calculate 95% CIs
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(tree_cover_lower_mean = quantile(tree_cover_mean,probs = c(0.025)),
                   tree_cover_higher_mean = quantile(tree_cover_mean,probs = c(0.975)),
                   tree_cover_lower_median = quantile(tree_cover_median,probs = c(0.025)),
                   tree_cover_higher_median = quantile(tree_cover_median,probs = c(0.975)),
                   tree_cover_lower_max = quantile(tree_cover_max,probs = c(0.025)),
                   tree_cover_higher_max = quantile(tree_cover_max,probs = c(0.975))) %>% 
  dplyr::ungroup()  %>% 
  tidyr::pivot_longer(!bin_age, names_to = "lower_upper", values_to = "tree_cover_5_95") %>%  
  dplyr::mutate(column_label = dplyr::if_else(stringr::str_detect(lower_upper, "tree_cover_lower"), nreps+2, nreps+3)) #need to add for plotting
rio::export(bin_200_ap_tree_boots_5_95, "data/intermediate_output/vegetation/bin_200_ap_tree_boots_5_95.csv")

bin_200_mean_ap_tree_plot <- bin_200_mean_ap_tree %>% 
  dplyr::mutate(column_label = nreps+1)
rio::export(bin_200_mean_ap_tree_plot, "data/intermediate_output/vegetation/bin_200_mean_ap_tree_plot.csv")

bin_200_median_ap_tree_plot <- bin_200_median_ap_tree %>% 
  dplyr::mutate(column_label = nreps+1) #need to include within plot
rio::export(bin_200_median_ap_tree_plot, "data/intermediate_output/vegetation/bin_200_median_ap_tree_plot.csv")

bin_200_max_ap_tree_plot <- bin_200_max_ap_tree %>% 
  dplyr::mutate(column_label = nreps+1)
rio::export(bin_200_max_ap_tree_plot, "data/intermediate_output/vegetation/bin_200_max_ap_tree_plot.csv")

bin_200_ap_tree_boots_5_95_mean <- bin_200_ap_tree_boots_5_95 %>% 
  dplyr::filter(stringr::str_detect(lower_upper, "mean"))
rio::export(bin_200_ap_tree_boots_5_95_mean, "data/intermediate_output/vegetation/bin_200_ap_tree_boots_5_95_mean.csv")

bin_200_ap_tree_boots_5_95_median <- bin_200_ap_tree_boots_5_95 %>% 
  dplyr::filter(stringr::str_detect(lower_upper, "median"))
rio::export(bin_200_ap_tree_boots_5_95_median, "data/intermediate_output/vegetation/bin_200_ap_tree_boots_5_95_median.csv")

bin_200_ap_tree_boots_5_95_max <- bin_200_ap_tree_boots_5_95 %>% 
  dplyr::filter(stringr::str_detect(lower_upper, "max"))
rio::export(bin_200_ap_tree_boots_5_95_max, "data/intermediate_output/vegetation/bin_200_ap_tree_boots_5_95_max.csv")



#Locfit smoothing of tree cover through time
f_locfit_tree_average(bin_200_median_ap_tree)
rio::export(bin_200_median_ap_tree_loc, "data/intermediate_output/vegetation/bin_200_median_ap_tree_loc.csv")



##Model bootstraps for prediction interval
## Modelled fit
epd_downcore_input <- rio::import( "data/intermediate_output/vegetation/epd_downcore_input.csv")

boot_pred_EPD_tree_ls <- lapply(boot_tree_model, function(x){ #fit reconstruction tree cover values based on bootstrapped models
  as.numeric(betareg::predict(x, newdata = dplyr::select(epd_downcore_input, ap_cover, elevation, needle_share, tree_shannon, site_model, sp_cover),type="response")) #modelled
})
boot_res_EPD_tree_ls <- lapply(boot_tree_model, function(x){ #fit reconstruction tree cover values based on bootstrapped models
 set.seed(42)
  as.numeric(sample(x$residuals, nrow(epd_downcore_input), replace = TRUE))
})

boot_predres_EPD_tree_ls <- purrr::map2(boot_pred_EPD_tree_ls, boot_res_EPD_tree_ls, `+`)

boot_pred_EPD_tree_refit_ls <- purrr::map2(boot_predres_EPD_tree_ls, boot_tree_qmap_model_ssplin, qmap::doQmapSSPLIN) #Apply respective quantile adjustment to reconstructions

boot_recon_tree_ls <- lapply(boot_pred_EPD_tree_refit_ls, function(x){ #add input data
  cbind(epd_downcore_input, x)  %>% 
    dplyr::filter(elevation < 1000) %>%  #Limit to records below 1000m
    dplyr::rename(tree_cover = dplyr::last_col()) %>%
    dplyr::mutate(tree_cover = dplyr::if_else(tree_cover>1, 1, tree_cover)) %>% #Because otherwise we may have re-fitted values greater than 100%
    dplyr::mutate(tree_cover = 100*tree_cover) %>%
    dplyr::left_join(dplyr::select(pollen_sample_ages, ID_SAMPLE, entity_name, latitude, longitude, median, site_type), by = "ID_SAMPLE")  %>%
    dplyr::filter(site_type %in% c("lake", "terrestrial bog/mire/fen", "terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, bog/mire/fen", "terrestrial, marsh")) 
  
})
boot_recon_tree <- dplyr::bind_rows(boot_recon_tree_ls, .id = "boot")  #reduce to single df
saveRDS(boot_recon_tree, "data/intermediate_output/vegetation/boot_recon_tree.rda")

#Bin data into 200yrs bins
boot_recon_tree <- readRDS("data/intermediate_output/vegetation/boot_recon_tree.rda")
boot_recon_tree_binned_200 <- boot_recon_tree %>%  #bin data by age, entity and boot
  dplyr::rename(date_value = median) %>% 
  dplyr::mutate(bin = ggplot2::cut_width(date_value, width = 200, center = 0, labels = F)) %>%    #place in 200 year bins
  dplyr::mutate(bin_age = (bin * 200) - 200) %>%
  dplyr::filter(bin_age <=14000)  %>% 
  dplyr::group_by(entity_name, bin_age, boot) %>% 
  dplyr::summarise(tree_cover = mean(tree_cover)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(pollen_sample_ages, entity_name, latitude, longitude), by = "entity_name", multiple = "first") %>% 
  dplyr::distinct()

boot_bin_200_median_recon_tree <- boot_recon_tree_binned_200 %>% 
  dplyr::group_by(boot, bin_age) %>% 
  dplyr::summarise(median_tree_cover = median(tree_cover)) %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(column_label = boot)
saveRDS(boot_bin_200_median_recon_tree, "data/intermediate_output/vegetation/boot_bin_200_median_recon_tree.rda")

boot_bin_200_median_recon_tree_5_95 <- boot_bin_200_median_recon_tree %>% #Calculate 95% CIs
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(tree_cover_lower_median = quantile(median_tree_cover,probs = c(0.025)),
                   tree_cover_higher_median = quantile(median_tree_cover,probs = c(0.975))) %>% 
  dplyr::ungroup()  %>% 
  tidyr::pivot_longer(!bin_age, names_to = "lower_upper", values_to = "tree_cover_5_95") %>%  
  dplyr::mutate(column_label = dplyr::if_else(stringr::str_detect(lower_upper, "tree_cover_lower"), 1002, 1003)) #need to add for plotting
rio::export(boot_bin_200_median_recon_tree_5_95, "data/intermediate_output/vegetation/boot_bin_200_median_recon_tree_5_95.csv")

boot_recon_tree_SE <- boot_recon_tree %>% #calculate SE based on bootstraps
  dplyr::group_by(ID_SAMPLE) %>%
  dplyr::summarise(mean_tree_cover = mean(tree_cover, na.rm = TRUE),
                   sd_tree_cover = sd(tree_cover, na.rm = TRUE),
                   n = n()) %>%
  dplyr::mutate(se_tree_cover = sd_tree_cover / sqrt(n)) %>%
  dplyr::ungroup()

boot_recon_tree_SE_binned_200 <- boot_recon_tree_SE %>% #SE
  dplyr::left_join(dplyr::select(boot_recon_tree, ID_SAMPLE:site_model,entity_name:site_type),  by = "ID_SAMPLE", multiple = "first") %>% 
  dplyr::distinct() %>% 
  dplyr::rename(date_value = median) %>% 
  dplyr::mutate(bin = ggplot2::cut_width(date_value, width = 200, center = 0, labels = F)) %>%    #place in 200 year bins
  dplyr::mutate(bin_age = (bin * 200) - 200) %>%
  dplyr::filter(bin_age <=14000)  %>% 
  dplyr::group_by(entity_name, bin_age) %>% 
  dplyr::summarise(se_tree_cover = mean(se_tree_cover)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(pollen_sample_ages, entity_name, latitude, longitude), by = "entity_name", multiple = "first") %>% 
  dplyr::distinct()

saveRDS(boot_recon_tree_SE_binned_200, "data/intermediate_output/vegetation/boot_recon_tree_SE_binned_200.rda")



########################
###Rasters
##Non-adjusted
#200 year bins
recon_tree_sf_200 <- recon_tree_binned_200 %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>%
  sf::st_transform(crs = 3035) 

recon_tree_rast_200 <- recon_tree_sf_200 %>% 
  sf::st_coordinates() %>% 
  raster::rasterize(., euro_map_3035_rast_50, fun = mean, field = recon_tree_sf_200$tree_cover) #50km2 grid
raster::plot(recon_tree_rast_200)
terra::global(terra::rast(recon_tree_rast_200), fun = "notNA")

recon_tree_200 <- ls()
recon_tree_200 <- parallel::mclapply(unique(recon_tree_sf_200$bin_age), function(x){
  dplyr::filter(recon_tree_sf_200, bin_age == x)
}, mc.cores = 8)
names(recon_tree_200) <- unique(recon_tree_sf_200$bin_age)
recon_tree_rast_200 <- ls()
recon_tree_rast_200 <- lapply(recon_tree_200, function(x){
  sf::st_coordinates(x) %>%    
  raster::rasterize(., euro_map_3035_rast_50, fun = mean, field = x$tree_cover)
})

terra::plot(recon_tree_rast_200[[1]])
recon_tree_spatrast_200 <- lapply(recon_tree_rast_200, terra::rast)
recon_tree_spatrast_200_stack <- terra::rast(recon_tree_spatrast_200)

terra::writeRaster(recon_tree_spatrast_200_stack, "data/intermediate_output/vegetation/raster_tree/recon_tree_rast_200.tif", overwrite = TRUE)

classification_matrix <- dplyr::tibble(from = c(0,10,20,30,40,50,60,70,80,90), 
                                       to = c(10,20,30,40,50,60,70,80,90,100),
                                       becomes = c(5, 15, 25, 35, 45, 55, 65, 75, 85, 95)) %>% 
  as.matrix()

recon_tree_rast_200_class <- ls() #classified version of above
recon_tree_rast_200_class <- lapply(recon_tree_rast_200, function(x){
  terra::rast(x) %>% 
    terra::classify(., classification_matrix, include.lowest = TRUE, right = TRUE)
})

recon_tree_rast_200_plotsetup <- ls()
recon_tree_rast_200_plotsetup <- lapply(recon_tree_rast_200, function(q){
  raster::as.data.frame(q, xy = TRUE) %>% 
    na.omit()
})

recon_tree_rast_200_plotsetup_4258 <- ls()
recon_tree_rast_200_plotsetup_4258 <- lapply(recon_tree_rast_200_plotsetup, function(q){
  q %>% 
    sf::st_as_sf(coords = c("x","y"), agr = "constant",crs = 3035) %>% 
    sf::st_transform(crs = 4258) %>% 
    dplyr::mutate(longitude = sf::st_coordinates(.)[,1], latitude = sf::st_coordinates(.)[,2]) %>% 
    sf::st_drop_geometry()
})

recon_tree_rast_200_plot <- ls()
recon_tree_rast_200_plot <- lapply(recon_tree_rast_200_plotsetup, function(q){
  ggplot2::ggplot(data = euro_map_3035)+
    geom_sf(fill = "seashell") +
    theme(panel.background = element_rect(fill = "aliceblue")) +
    geom_tile(data = q, aes(x = x, y = y, fill = layer))+
    scale_fill_viridis_c()
})

recon_tree_rast_200_plot$`1000`

recon_tree_rast_200_class_plotsetup <- ls()
recon_tree_rast_200_class_plotsetup <- lapply(recon_tree_rast_200_class, function(q){
  terra::as.data.frame(q, xy = TRUE) %>% 
    na.omit() %>% 
    dplyr::mutate(layer = as.character(layer)) %>% 
    dplyr::mutate(layer = dplyr::if_else(layer == "5", "0-10", layer)) %>% 
    dplyr::mutate(layer = dplyr::if_else(layer == "15", "10-20", layer)) %>%
    dplyr::mutate(layer = dplyr::if_else(layer == "25", "20-30", layer)) %>%
    dplyr::mutate(layer = dplyr::if_else(layer == "35", "30-40", layer)) %>%
    dplyr::mutate(layer = dplyr::if_else(layer == "45", "40-50", layer)) %>%
    dplyr::mutate(layer = dplyr::if_else(layer == "55", "50-60", layer)) %>%
    dplyr::mutate(layer = dplyr::if_else(layer == "65", "60-70", layer)) %>%
    dplyr::mutate(layer = dplyr::if_else(layer == "75", "70-80", layer)) %>%
    dplyr::mutate(layer = dplyr::if_else(layer == "85", "80-90", layer)) %>%
    dplyr::mutate(layer = dplyr::if_else(layer == "95", "90-100", layer)) %>% 
    dplyr::mutate(layer = factor(layer, levels = c("0-10",
                                                   "10-20",
                                                   "20-30",
                                                   "30-40",
                                                   "40-50",
                                                   "50-60",
                                                   "60-70",
                                                   "70-80",
                                                   "80-90",
                                                   "90-100")))
    
})


recon_tree_rast_200_class_plotsetup <- recon_tree_rast_200_class_plotsetup[c("12600","12400","12200","12000", #clip and reorder
                                                                             "11800","11600","11400","11200",
                                                                             "11000","10800","10600","10400",
                                                                             "10200","10000","9800","9600",
                                                                             "9400","9200","9000","8800",
                                                                             "8600","8400","8200","8000",
                                                                             "7800","7600","7400","7200",
                                                                             "7000","6800","6600","6400",
                                                                             "6200","6000","5800","5600",
                                                                             "5400","5200","5000","4800",
                                                                             "4600","4400","4200","4000",
                                                                             "3800","3600","3400","3200",
                                                                             "3000","2800","2600","2400",
                                                                             "2200","2000","1800","1600",
                                                                             "1400","1200","1000","800",
                                                                             "600","400","200",'0')]
raster_limits <- c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100")

recon_tree_rast_200_class_plot <- ls()
recon_tree_rast_200_class_plot <- lapply(seq_along(recon_tree_rast_200_class_plotsetup), function(q){
  ggplot2::ggplot(data = euro_map_3035)+
    geom_sf(fill = "seashell") +
    theme(panel.background = element_rect(fill = "aliceblue")) +
    geom_tile(data = recon_tree_rast_200_class_plotsetup[[q]], aes(x = x, y = y, fill = factor(layer))) +
    viridis::scale_fill_viridis(discrete = TRUE, name = "Tree cover %", drop = FALSE, option = "C", limits = raster_limits)+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          axis.ticks.y=element_blank(),
          legend.key.size = unit(0.5, "cm"), legend.key = element_rect(color = NA, fill = NA),
          legend.title=element_text(size=9), 
          legend.text=element_text(size=8))+
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))+
    ggtitle(paste0("Bin age ", names(recon_tree_rast_200_class_plotsetup[q]), " cal. BP"))+
    theme(plot.title = element_text(size=9))
})
names(recon_tree_rast_200_class_plot) <- names(recon_tree_rast_200_class_plotsetup)
saveRDS(recon_tree_rast_200_class_plot, "data/intermediate_output/vegetation/raster_tree/recon_tree_rast_200_class_plot.rda")


recon_tree_rast_200_class_plot$`12000`

recon_tree_rast_200_class_plot_regroup <- ls()
recon_tree_rast_200_class_plot_regroup <- split(recon_tree_rast_200_class_plot, ceiling(seq_along(recon_tree_rast_200_class_plot)/4))
recon_tree_rast_200_class_plot_regroup_im <- ls()
recon_tree_rast_200_class_plot_regroup_im <- lapply(recon_tree_rast_200_class_plot_regroup, function(x){
  gridExtra::grid.arrange(grobs = x, ncol = 2) %>% 
    ggsave(paste0("figs/tree_cover_rasters/recon_tree_rast_200_", names(x)[1], "to", tail(names(x), n=1), ".pdf"),., height = 12, width = 16, unit = "cm")
})





#####################
#Matching to Zanon - 5arc minutes, 250 year
#Non-adjusted
#250 year bins
recon_tree_binned_250 <- recon_tree %>% #readjusted bins to match Zanon
  dplyr::mutate(bin = ggplot2::cut_width(median, width = 250, center = 0, labels = F)) %>%   #place in 250 year bins
  dplyr::mutate(bin_age = (bin - 1)*250) %>% 
  dplyr::filter(bin_age <=14000)  %>% 
  dplyr::group_by(entity_name, bin_age) %>% 
  dplyr::summarise(tree_cover = mean(tree_cover), number = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(pollen_sample_ages, entity_name, latitude, longitude), by = "entity_name") %>% 
  dplyr::distinct()

bin_250_mean_recon_tree <- recon_tree_binned_250 %>% #Calculate mean through time
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(tree_cover_mean = mean(tree_cover), tree_cover_sd = sd(tree_cover)) %>% 
  dplyr::ungroup()

bin_250_median_recon_tree <- recon_tree_binned_250 %>% #Calculate median through time
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(tree_cover_median = median(tree_cover), tree_cover_lower = quantile(tree_cover, probs = 0.25),tree_cover_higher = quantile(tree_cover, probs = 0.75) ) %>% 
  dplyr::ungroup()

bin_250_n_recon_tree <- recon_tree_binned_250 %>% #Number of records through time per bin
  dplyr::select(-tree_cover, -number) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(number_records = n()) %>% 
  dplyr::ungroup()

#Extract values from Zanon files
recon_tree_binned_250_coords <- recon_tree_binned_250 %>% #Get record locations
  dplyr::select(longitude, latitude) %>% 
  dplyr::distinct() %>% 
  as.matrix()

recon_tree_binned_250_coords_buf <- recon_tree_binned_250 %>% #Add buffer
  dplyr::select(longitude, latitude, entity_name) %>% 
  dplyr::distinct() %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>% 
  sf::st_buffer(., 5000)

f_ZANON_extract <- function(b){ #this function generates a raster based on saved files then extracts the values associated with the sites from the stack
  raster_zan <- terra::rast(b)
  terra::crs(raster_zan) <- "EPSG:4258"
  ZANON_extract <- exactextractr::exact_extract(raster_zan, recon_tree_binned_250_coords_buf,'mean', force_df = TRUE)
  return(ZANON_extract)
  
}

test <- f_ZANON_extract("./data/input/ZANON/forest_cover/forest_cover_0.grd")

f_ZANON_extract_ID <- function(b){ #this function generates a raster based on saved files then extracts the cell id for sites
  raster_zan <- terra::rast(b)
  terra::crs(raster_zan) <- "EPSG:4258"
  ZANON_extract_ID <- terra::cellFromXY(raster_zan, recon_tree_binned_250_coords) %>% 
    dplyr::as_tibble()
  return(ZANON_extract_ID)
  
}

test1 <- f_ZANON_extract_ID("./data/input/ZANON/forest_cover/forest_cover_0.grd")


ZANON_forest <- data.frame()
for (i in list.files(path = "./data/input/ZANON/forest_cover/", full.names = T)){
  print(i)
  extracted <- f_ZANON_extract(i) %>%  #Buffer or not?
    dplyr::as_tibble() %>% 
    dplyr::mutate(year = stringr::str_extract(i, "(\\d+)")) %>% #Get year
    dplyr::rename(forest_cover = 1) %>% # Rename columns
    dplyr::mutate(entity_name = recon_tree_binned_250_coords_buf$entity_name, .before = 1)  #Add entity, lat and lon data
  ZANON_forest <- rbind(ZANON_forest, extracted)
  
}

ZANON_forest_ID <- data.frame()
for (i in list.files(path = "./data/input/ZANON/forest_cover/", full.names = T)){
  print(i)
  extracted <- f_ZANON_extract_ID(i) %>%  #Get cell IDs
    dplyr::as_tibble() %>% 
    dplyr::mutate(year = stringr::str_extract(i, "(\\d+)")) %>% #Get year
    dplyr::rename(cell_ID = 1) %>% # Rename columns
    dplyr::mutate(entity_name = recon_tree_binned_250_coords_buf$entity_name, .before = 1)  #Add entity, lat and lon data
  ZANON_forest_ID <- rbind(ZANON_forest_ID, extracted)
  
}

zan_forest <- ZANON_forest %>% 
  dplyr::mutate(bin_age = as.numeric(year)) %>% 
  dplyr::arrange(bin_age, entity_name) %>% 
  dplyr::rename(zan_tree = forest_cover)

zan_cell_id <- ZANON_forest_ID %>% 
  dplyr::mutate(bin_age = as.numeric(year))

#Comparison at just locations when we have data
comp_zanon <- recon_tree_binned_250 %>% 
  dplyr::left_join(dplyr::select(zan_forest, bin_age, entity_name, zan_tree), by = c("bin_age", "entity_name")) %>% 
  tidyr::drop_na() %>% 
  dplyr::left_join(zan_cell_id, by = c("entity_name", "bin_age")) %>% #Add cell ID
  dplyr::left_join(dplyr::select(epd_bioregion_map, entity_name, bioregion_name), by = "entity_name") #Add bioregion name

comp_zanon_cell <- comp_zanon %>% 
  dplyr::group_by(bin_age, cell_ID) %>%
  dplyr::summarise(zan_tree = mean(zan_tree), tree_cover = mean(tree_cover), bioregion_name = statip::mfv1(bioregion_name)) %>% #First get average value per cell
  dplyr::ungroup()
  
comp_zanon_median <- comp_zanon %>% #Calculate median values for all record locations (at 250 year intervals)
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(zan = median(zan_tree), tree_recon = median(tree_cover)) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_longer(!bin_age, names_to = "source", values_to = "tree_cover")

comp_zanon_median_cell <- comp_zanon %>% #for a single value per cell
  dplyr::group_by(bin_age, cell_ID) %>%
  dplyr::summarise(zan_tree = mean(zan_tree), tree_cover = mean(tree_cover)) %>% #First get average value per cell
  dplyr::ungroup() %>% 
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(zan = median(zan_tree), tree_recon = median(tree_cover)) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_longer(!bin_age, names_to = "source", values_to = "tree_cover")


#Scatter plots
cor(comp_zanon$tree_cover, comp_zanon$zan_tree)
comp_zanon_atlantic <- comp_zanon %>% 
  dplyr::filter(bioregion_name == "Atlantic")
cor(comp_zanon_atlantic$tree_cover, comp_zanon_atlantic$zan_tree)
comp_zanon_boreal <- comp_zanon %>% 
  dplyr::filter(bioregion_name == "Boreal")
cor(comp_zanon_boreal$tree_cover, comp_zanon_boreal$zan_tree)
comp_zanon_continental <- comp_zanon %>% 
  dplyr::filter(bioregion_name == "Continental")
cor(comp_zanon_continental$tree_cover, comp_zanon_continental$zan_tree)
comp_zanon_mediterranean <- comp_zanon %>% 
  dplyr::filter(bioregion_name == "Mediterranean")
cor(comp_zanon_mediterranean$tree_cover, comp_zanon_mediterranean$zan_tree)

data("entity", package = "smpds")

SMPDS_entity <- entity
SMPDS_entity_lookup <- SMPDS_entity %>% 
  dplyr::select(entity_name, ID_ENTITY) %>% 
  dplyr::distinct()



#Comparison with all locations
zan_forest_all <- zan_forest %>% 
  dplyr::rename(tree_cover = zan_tree) %>% 
  dplyr::mutate(source = "zanon") %>% 
  dplyr::select(entity_name, bin_age, tree_cover, source)

comp_zanon_all_locations <- recon_tree_binned_250 %>%
  dplyr::mutate(source = "recon") %>% 
  dplyr::select(entity_name, bin_age, tree_cover, source) %>% 
  rbind(zan_forest_all) %>% 
  tidyr::drop_na() %>%
  dplyr::left_join(zan_cell_id, by = c("entity_name", "bin_age")) %>% #Add cell ID
  dplyr::left_join(dplyr::select(epd_bioregion_map, entity_name, bioregion_name), by = "entity_name") #Add bioregion name

comp_zanon_cell_all_locations <- comp_zanon_all_locations %>% 
  dplyr::group_by(bin_age, cell_ID, source) %>%
  dplyr::summarise(tree_cover = mean(tree_cover), bioregion_name = statip::mfv1(bioregion_name)) %>% #First get average value per cell
  dplyr::ungroup()

comp_zanon_median_all_locations <- comp_zanon_all_locations %>% #Calculate median values for record locations where we have data at those time periods (at 250 year intervals)
  dplyr::group_by(bin_age,source) %>% 
  dplyr::summarise(tree_cover = median(tree_cover)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(column_label = 3001)
rio::export(comp_zanon_median_all_locations, "data/intermediate_output/vegetation/comp_zanon_median_all_locations.csv") 
  

comp_zanon_median_cell_all_locations <- comp_zanon_all_locations %>% #for a single value per cell
  dplyr::group_by(bin_age, cell_ID, source) %>%
  dplyr::summarise(tree_cover = median(tree_cover)) %>% #First get average value per cell
  dplyr::ungroup() %>% 
  dplyr::group_by(bin_age, source) %>% 
  dplyr::summarise(tree_cover = median(tree_cover)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(column_label = 3001)
rio::export(comp_zanon_median_cell_all_locations, "data/intermediate_output/vegetation/comp_zanon_median_cell_all_locations.csv") 




#####################
#Matching Serge - 1degree (approx 100km2), different years
#Non-adjusted
#Custom bins
recon_tree_binned_serge <- recon_tree %>% 
  dplyr::mutate(bin = cut(median, c(0,100,350,700,1200,1700,2200,2700,3200,3700,4200,4700,5200,5700,6200,
                                    6700,7200,7700,8200,8700,9200,9700,10200,10700,11200,11700,Inf),
                          labels = c("0-100", "100-350", "350-700", "700-1200","1200-1700",
                                     "1700-2200","2200-2700","2700-3200","3200-3700","3700-4200",
                                     "4200-4700","4700-5200","5200-5700","5700-6200","6200-6700",
                                     "6700-7200","7200-7700","7700-8200","8200-8700","8700-9200",
                                     "9200-9700","9700-10200","10200-10700","10700-11200","11200-11700", "11700-last"))) %>% 
  dplyr::filter(bin !="11700-last")  %>% 
  dplyr::group_by(entity_name, bin) %>% 
  dplyr::summarise(tree_cover = median(tree_cover), number = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(pollen_sample_ages, entity_name, latitude, longitude), by = "entity_name") %>% 
  dplyr::distinct()

bin_serge_mean_recon_tree <- recon_tree_binned_serge %>% #Calculate mean through time
  dplyr::group_by(bin) %>% 
  dplyr::summarise(tree_cover_mean = mean(tree_cover), tree_cover_sd = sd(tree_cover)) %>% 
  dplyr::ungroup()

bin_serge_n_recon_tree <- recon_tree_binned_serge %>% 
  dplyr::select(-tree_cover, -number) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(bin) %>% 
  dplyr::summarise(number_records = n()) %>% 
  dplyr::ungroup()

#Extract values from Serge files
recon_tree_binned_serge_coords <- recon_tree_binned_serge %>% 
  dplyr::select(longitude, latitude) %>% 
  dplyr::distinct() %>% 
  as.matrix()

recon_tree_binned_serge_coords_buf <- recon_tree_binned_serge %>% 
  dplyr::select(longitude, latitude, entity_name) %>% 
  dplyr::distinct() %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>% 
  sf::st_buffer(., 5000)

serge_ls <- rio::import_list(Sys.glob("./data/input/Serge/TERRA_RVresults_RPPs.st1/RV_mean_RPPs.st1/*.csv"))
names(serge_ls) <- c("0-100", "3700-4200", "4200-4700","4700-5200","5200-5700","5700-6200",
                        "6200-6700","6700-7200","7200-7700","7700-8200","8200-8700","100-350",
                        "8700-9200","9200-9700","9700-10200","10200-10700","10700-11200",
                        "11200-11700","350-700", "700-1200","1200-1700","1700-2200",
                        "2200-2700","2700-3200","3200-3700")
serge_tree <- lapply(serge_ls, function(x){
  x %>% 
    dplyr::mutate(serge_tree = NBTT+EMT+BSBT+BSTCT+BSTWT) %>% 
    dplyr::mutate(temp_sum = rowSums(across(Abies.alba:serge_tree))) %>% #exclude no vegetation
    dplyr::filter(temp_sum != 0) %>% 
    dplyr::select(LonDD, LatDD, serge_tree)
})

serge_tree_rast <- lapply(serge_tree, function(x){
  terra::rast(x, type = "xyz", crs = "EPSG:4258")
})

serge_tree_extract_ls <- lapply(serge_tree_rast, function(x){
  exactextractr::exact_extract(x, recon_tree_binned_serge_coords_buf, 'mean', force_df = TRUE) %>% 
    dplyr::mutate(entity_name = recon_tree_binned_serge_coords_buf$entity_name, .before = 1)
})

serge_tree_extract <- bind_rows(serge_tree_extract_ls, .id = "bin") %>% 
  dplyr::rename(serge_tree = mean)

serge_cell_id_extract_ls <- lapply(serge_tree_rast, function(x){
  terra::cellFromXY(x, recon_tree_binned_serge_coords) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(entity_name = recon_tree_binned_serge_coords_buf$entity_name, .before = 1)
})

serge_cell_id_extract <- bind_rows(serge_cell_id_extract_ls, .id = "bin") %>%   #toidentify where the same cell is being used
  dplyr::rename(cell_ID = value)

# Comparison at just locations when we have data
comp_serge <- recon_tree_binned_serge %>% 
  dplyr::left_join(dplyr::select(serge_tree_extract, bin, entity_name, serge_tree), by = c("bin", "entity_name")) %>%
  dplyr::left_join(serge_cell_id_extract, by = c("bin", "entity_name")) %>% 
  tidyr::drop_na() %>% 
  dplyr::mutate(bin1 = as.numeric(stringr::str_extract(bin, "[^-]+"))) %>% 
  dplyr::mutate(bin2 = as.numeric(stringr::str_extract(bin, '(?<=-).*'))) %>% 
  dplyr::mutate(bin_age = (bin1 + bin2)/2) %>% 
  dplyr::filter(serge_tree != 0)  %>% 
  dplyr::left_join(dplyr::select(epd_bioregion_map, entity_name, bioregion_name), by = "entity_name") #Add bioregion name

comp_serge_cell <- comp_serge %>% 
  dplyr::group_by(bin_age, cell_ID) %>% 
  dplyr::summarise(serge_tree = mean(serge_tree), tree_cover = mean(tree_cover), bioregion_name = statip::mfv1(bioregion_name)) %>% 
  dplyr::ungroup()

comp_serge_median <- comp_serge %>% 
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(serge = median(serge_tree), tree_recon = median(tree_cover)) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_longer(!bin_age, names_to = "source", values_to = "tree_cover")

comp_serge_median_cell <- comp_serge_cell %>% 
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(serge = median(serge_tree), tree_recon = median(tree_cover)) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_longer(!bin_age, names_to = "source", values_to = "tree_cover")

comp_serge_n <- comp_serge %>% 
  dplyr::select(bin_age, entity_name) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(number_records = n()) %>% 
  dplyr::ungroup()

cor(comp_serge$tree_cover, comp_serge$serge_tree)
cor(comp_serge_cell$tree_cover, comp_serge_cell$serge_tree)


#Comparison at all locations
serge_tree_extract_all <- serge_tree_extract %>% 
  dplyr::rename(tree_cover = serge_tree) %>% 
  dplyr::mutate(source = "serge") %>% 
  dplyr::mutate(bin1 = as.numeric(stringr::str_extract(bin, "[^-]+"))) %>% 
  dplyr::mutate(bin2 = as.numeric(stringr::str_extract(bin, '(?<=-).*'))) %>% 
  dplyr::mutate(bin_age = (bin1 + bin2)/2) %>% 
  dplyr::select(entity_name, bin_age, tree_cover, source) %>% 
  tidyr::drop_na() %>% 
  dplyr::filter(tree_cover !=0)

serge_cell_id_extract_all <- serge_cell_id_extract %>% 
  dplyr::mutate(bin1 = as.numeric(stringr::str_extract(bin, "[^-]+"))) %>% 
  dplyr::mutate(bin2 = as.numeric(stringr::str_extract(bin, '(?<=-).*'))) %>% 
  dplyr::mutate(bin_age = (bin1 + bin2)/2) %>% 
  dplyr::select(entity_name, bin_age, cell_ID)

comp_serge_all_locations <- recon_tree_binned_serge %>% 
  dplyr::mutate(source = "recon") %>% 
  dplyr::mutate(bin1 = as.numeric(stringr::str_extract(bin, "[^-]+"))) %>% 
  dplyr::mutate(bin2 = as.numeric(stringr::str_extract(bin, '(?<=-).*'))) %>% 
  dplyr::mutate(bin_age = (bin1 + bin2)/2) %>% 
  dplyr::select(entity_name, bin_age, tree_cover, source) %>%  
  rbind(serge_tree_extract_all) %>% 
  tidyr::drop_na() %>%
  dplyr::left_join(serge_cell_id_extract_all, by = c("bin_age", "entity_name")) %>% 
  tidyr::drop_na() %>% 
  dplyr::left_join(dplyr::select(epd_bioregion_map, entity_name, bioregion_name), by = "entity_name") #Add bioregion name

comp_serge_cell_all_locations <- comp_serge_all_locations %>% 
  dplyr::group_by(bin_age, cell_ID, source) %>% 
  dplyr::summarise(tree_cover = mean(tree_cover), bioregion_name = statip::mfv1(bioregion_name)) %>% 
  dplyr::ungroup()

comp_serge_median_all_locations <- comp_serge_all_locations %>% 
  dplyr::group_by(bin_age,source) %>% 
  dplyr::summarise(tree_cover = median(tree_cover)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(column_label = 2001)
rio::export(comp_serge_median_all_locations, "data/intermediate_output/vegetation/comp_serge_median_all_locations.csv")


comp_serge_median_cell_all_locations <- comp_serge_cell_all_locations %>% 
  dplyr::group_by(bin_age,source) %>% 
  dplyr::summarise(tree_cover = median(tree_cover)) %>% 
  dplyr::ungroup()  %>% 
  dplyr::mutate(column_label = 2001)
rio::export(comp_serge_median_cell_all_locations, "data/intermediate_output/vegetation/comp_serge_median_all_locations.csv")


#Comparison our reconstruction with other reconstructions, no changes
comp_all <- recon_tree_binned_200 %>% 
  dplyr::select(entity_name, bin_age, tree_cover) %>% 
  dplyr::mutate(source = "recon") %>% 
  dplyr::bind_rows(dplyr::select(dplyr::rename(zan_forest, tree_cover = zan_tree), entity_name, bin_age, tree_cover)) %>%  #add zanon values to our reconstructions, no adjustment bins or merging cells
  dplyr::mutate(source = dplyr::if_else(is.na(source), "zanon", source)) %>% 
  dplyr::bind_rows(dplyr::rename(dplyr::select(comp_serge, entity_name, bin_age, serge_tree), tree_cover = serge_tree)) %>%  #add serge values to our reconstructions, no adjustment bins or merging cells
  dplyr::mutate(source = dplyr::if_else(is.na(source), "serge", source)) %>% 
  tidyr::drop_na() %>% 
  dplyr::left_join(entity_bioregions, by = "entity_name") %>%
  dplyr::filter(bin_age <= 12000)
  
comp_all_median <- comp_all %>% 
  dplyr::group_by(bin_age, source) %>% 
  dplyr::summarise(tree_cover = median(tree_cover)) %>% 
  dplyr::ungroup()

comp_all_mean <- comp_all %>% 
  dplyr::group_by(bin_age, source) %>% 
  dplyr::summarise(tree_cover = mean(tree_cover)) %>% 
  dplyr::ungroup()

comp_all_n <- comp_all %>% 
  dplyr::group_by(bin_age, source) %>% 
  dplyr::summarise(n()) %>% 
  dplyr::ungroup()

comp_all_median_ls <- comp_all_median %>% 
  tidytable::group_split(source, .keep = FALSE, .named = TRUE)
comp_all_mean_ls <- comp_all_mean %>% 
  tidytable::group_split(source, .keep = FALSE, .named = TRUE)

f_locfit_tree_average_comp <- function(d){
  tree_loc <- d %>% 
    dplyr::select(bin_age, 2) %>% 
    dplyr::rename(tree_cover = 2)
  x <- as.vector(tree_loc$bin_age)
  y <- as.vector(tree_loc$tree_cover)
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
  loc01 <- locfit::locfit(y ~ locfit::lp(x, deg=1, h=1000), maxk=800, family="qrgauss")
  pred01 <- predict(loc01, newdata=tree_loc$bin_age, se.fit=TRUE)
  locfit <- data.frame(tree_loc$bin_age, tree_loc$tree_cover, pred01$fit)
  colnames(locfit) <- c("bin_age", "None", "1000-year")
  locfit_long <- locfit %>% 
    tidyr::pivot_longer(!bin_age, names_to = "half_width", values_to = "average") %>% 
    dplyr::mutate(half_width = as.factor(half_width))
}

comp_all_median_loc <- lapply(comp_all_median_ls, f_locfit_tree_average_comp) %>% 
  dplyr::bind_rows(., .id = "source")
rio::export(comp_all_median_loc, "data/intermediate_output/vegetation/comp_all_median_loc.csv")
comp_all_mean_loc <- lapply(comp_all_mean_ls, f_locfit_tree_average_comp) %>% 
  dplyr::bind_rows(., .id = "source")
rio::export(comp_all_mean_loc, "data/intermediate_output/vegetation/comp_all_mean_loc.csv")




#Group 4 bioregion
comp_all_bioregiongr4 <- comp_all %>% 
  dplyr::filter(bioregion_name %in% c("Atlantic", "Boreal", "Continental", "Mediterranean"))

comp_all_bioregiongr4_median <- comp_all_bioregiongr4 %>% 
  dplyr::group_by(bin_age, source, bioregion_name) %>% 
  dplyr::summarise(tree_cover = median(tree_cover)) %>% 
  dplyr::ungroup()

comp_all_bioregiongr4_mean <- comp_all_bioregiongr4 %>% 
  dplyr::group_by(bin_age, source, bioregion_name) %>% 
  dplyr::summarise(tree_cover = mean(tree_cover)) %>% 
  dplyr::ungroup()

comp_all_bioregiongr4_n <- comp_all_bioregiongr4 %>% 
  dplyr::group_by(bin_age, source, bioregion_name) %>% 
  dplyr::summarise(number_records = n()) %>% 
  dplyr::ungroup()

comp_all_bioregiongr4_median_EU <- comp_all_bioregiongr4 %>% 
  dplyr::group_by(bin_age, source) %>% 
  dplyr::summarise(tree_cover = median(tree_cover)) %>% 
  dplyr::ungroup()

comp_all_bioregiongr4_eqw_median_EU <- comp_all_bioregiongr4_median %>%
  dplyr::group_by(bin_age, source) %>% 
  dplyr::summarise(tree_cover = median(tree_cover)) %>% 
  dplyr::ungroup()

comp_all_bioregiongr4_median_ls <- comp_all_bioregiongr4_median %>% 
  dplyr::select(bin_age, bioregion_name, tree_cover, source) %>% 
  dplyr::mutate(grouping = paste0(source, "_", bioregion_name)) %>% 
  tidytable::group_split(grouping, .keep = TRUE, .named = TRUE)

comp_all_bioregiongr4_mean_ls <- comp_all_bioregiongr4_mean %>% 
  dplyr::select(bin_age, bioregion_name,tree_cover, source) %>% 
  dplyr::mutate(grouping = paste0(source, "_", bioregion_name)) %>% 
  tidytable::group_split(grouping, .keep = TRUE, .named = TRUE)

comp_all_bioregiongr4_n_ls <- comp_all_bioregiongr4_n %>% 
  dplyr::select(bin_age, bioregion_name,number_records, source) %>% 
  dplyr::mutate(grouping = paste0(source, "_", bioregion_name)) %>% 
  tidytable::group_split(grouping, .keep = TRUE, .named = TRUE)

comp_all_bioregiongr4_median_EU_ls <- comp_all_bioregiongr4_median_EU %>% 
  dplyr::select(bin_age, tree_cover, source) %>% 
  tidytable::group_split(source, .keep = FALSE, .named = TRUE)

comp_all_bioregiongr4_eqw_median_EU_ls <- comp_all_bioregiongr4_median_EU %>% 
  dplyr::select(bin_age, tree_cover, source) %>% 
  tidytable::group_split(source, .keep = FALSE, .named = TRUE)



f_locfit_tree_average_bioregion <- function(d){
  tree_loc <- d %>% 
    dplyr::select(bin_age,2, 3) %>% 
    dplyr::rename(tree_cover = 3)
  x <- as.vector(tree_loc$bin_age)
  y <- as.vector(tree_loc$tree_cover)
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
  loc01 <- locfit::locfit(y ~ locfit::lp(x, deg=1, h=1000), maxk=800, family="qrgauss")
  pred01 <- predict(loc01, newdata=tree_loc$bin_age, se.fit=TRUE)
  locfit <- data.frame(tree_loc$bin_age, tree_loc$tree_cover, pred01$fit)
  colnames(locfit) <- c("bin_age", "None", "1000-year")
  locfit_long <- locfit %>% 
    tidyr::pivot_longer(!bin_age, names_to = "half_width", values_to = "average") %>% 
    dplyr::mutate(half_width = as.factor(half_width))
  return(locfit_long)
}
comp_all_bioregiongr4_median_loc <- lapply(comp_all_bioregiongr4_median_ls, f_locfit_tree_average_bioregion) %>% 
  dplyr::bind_rows(., .id = "grouping") %>% 
  dplyr::mutate(bioregion_name = stringr::str_extract(grouping, "(?<=_).*"), source = stringr::str_extract(grouping, ".*(?=_)"))
rio::export(comp_all_bioregiongr4_median_loc, "data/intermediate_output/vegetation/comp_all_bioregiongr4_median_loc.csv")

comp_all_bioregiongr4_mean_loc <- lapply(comp_all_bioregiongr4_mean_ls, f_locfit_tree_average_bioregion) %>% 
  dplyr::bind_rows(., .id = "grouping") %>% 
  dplyr::mutate(bioregion_name = stringr::str_extract(grouping, "(?<=_).*"), source = stringr::str_extract(grouping, ".*(?=_)"))
rio::export(comp_all_bioregiongr4_mean_loc, "data/intermediate_output/vegetation/comp_all_bioregiongr4_mean_loc.csv")

comp_all_bioregiongr4_n_loc <- lapply(comp_all_bioregiongr4_n_ls, f_locfit_tree_average_bioregion) %>% 
  dplyr::bind_rows(., .id = "grouping") %>% 
  dplyr::mutate(bioregion_name = stringr::str_extract(grouping, "(?<=_).*"), source = stringr::str_extract(grouping, ".*(?=_)")) %>% 
  dplyr::rename(number_records = average)
rio::export(comp_all_bioregiongr4_n_loc, "data/intermediate_output/vegetation/comp_all_bioregiongr4_n_loc.csv")

comp_all_bioregiongr4_median_EU_loc <- lapply(comp_all_bioregiongr4_median_EU_ls,f_locfit_tree_average_comp) %>% 
  dplyr::bind_rows(., .id = "source")
rio::export(comp_all_bioregiongr4_median_EU_loc, "data/intermediate_output/vegetation/comp_all_bioregiongr4_median_EU_loc.csv")

comp_all_bioregiongr4_eqw_median_EU_loc <- lapply(comp_all_bioregiongr4_eqw_median_EU_ls,f_locfit_tree_average_comp) %>% 
  dplyr::bind_rows(., .id = "source")
rio::export(comp_all_bioregiongr4_eqw_median_EU_loc, "data/intermediate_output/vegetation/comp_all_bioregiongr4_eqw_median_EU_loc.csv")



#Peak analysis
peak_comp_all_bioregiongr4_median_loc <- comp_all_bioregiongr4_median_loc %>% 
  dplyr::filter(half_width == "1000-year") %>% 
  dplyr::group_by(grouping) %>% 
  dplyr::filter(average == max(average)) %>% 
  dplyr::rename(peak = average) %>% 
  dplyr::ungroup()

peak_range_comp_all_bioregiongr4_median_loc <- comp_all_bioregiongr4_median_loc %>%
  dplyr::filter(half_width == "1000-year") %>% 
  dplyr::left_join(dplyr::select(peak_comp_all_bioregiongr4_median_loc, grouping, peak), by = "grouping") %>% 
  dplyr::mutate(range_calc = peak - average) %>% 
  dplyr::filter(dplyr::between(range_calc,-5,5)) %>% 
  dplyr::mutate(range_calc1 = 1-(abs(range_calc)/5)) %>% 
  dplyr::mutate(range_calc2 = dplyr::if_else(range_calc1 == 1, 1, range_calc1/2)) %>% 
  dplyr::mutate(grouping = as.factor(grouping)) %>% 
  dplyr::mutate(grouping = forcats::fct_relevel(grouping,
                                                "zanon_Mediterranean","zanon_Continental","zanon_Boreal","zanon_Atlantic",
                                                "serge_Mediterranean", "serge_Continental", "serge_Boreal","serge_Atlantic",    
                                                "recon_Mediterranean", "recon_Continental", "recon_Boreal","recon_Atlantic"))
rio::export(peak_range_comp_all_bioregiongr4_median_loc, "data/intermediate_output/vegetation/peak_range_comp_all_bioregiongr4_median_loc.csv")
  


#Comparison our reconstruction with other cell reconstructions, no changes
comp_all_cell <- recon_tree_binned_200 %>% 
  dplyr::select(entity_name, bin_age, tree_cover) %>% 
  dplyr::left_join(entity_bioregions, by = "entity_name") %>% 
  dplyr::mutate(source = "recon") %>%
  dplyr::rename(tree_cover_temp = tree_cover) %>% 
  dplyr::bind_rows(dplyr::select(dplyr::rename(comp_zanon_cell, tree_cover_temp = zan_tree), bin_age, tree_cover_temp, bioregion_name)) %>%  #add zanon values to our reconstructions, no adjustment bins or merging cells
  dplyr::mutate(source = dplyr::if_else(is.na(source), "zanon", source)) %>% 
  dplyr::bind_rows(dplyr::rename(dplyr::select(comp_serge_cell, bin_age, serge_tree, bioregion_name), tree_cover_temp = serge_tree)) %>%  #add serge values to our reconstructions, no adjustment bins or merging cells
  dplyr::mutate(source = dplyr::if_else(is.na(source), "serge", source)) %>% 
  dplyr::select(-entity_name) %>%
  dplyr::rename(tree_cover = tree_cover_temp) %>% 
  tidyr::drop_na()

comp_all_cell_median <- comp_all_cell %>% 
  dplyr::group_by(bin_age, source) %>% 
  dplyr::summarise(tree_cover = median(tree_cover)) %>% 
  dplyr::ungroup()

comp_all_cell_mean <- comp_all_cell %>% 
  dplyr::group_by(bin_age, source) %>% 
  dplyr::summarise(tree_cover = mean(tree_cover)) %>% 
  dplyr::ungroup()

comp_all_cell_n <- comp_all_cell %>% 
  dplyr::group_by(bin_age, source) %>% 
  dplyr::summarise(n()) %>% 
  dplyr::ungroup()


comp_all_cell_median_ls <- comp_all_cell_median %>% 
  tidytable::group_split(source, .keep = FALSE, .named = TRUE)
comp_all_cell_mean_ls <- comp_all_cell_mean %>% 
  tidytable::group_split(source, .keep = FALSE, .named = TRUE)

comp_all_cell_median_loc <- lapply(comp_all_cell_median_ls, f_locfit_tree_average_comp) %>% 
  dplyr::bind_rows(., .id = "source")
rio::export(comp_all_cell_median_loc,"data/intermediate_output/vegetation/comp_all_cell_median_loc.csv")
comp_all_cell_mean_loc <- lapply(comp_all_cell_mean_ls, f_locfit_tree_average_comp) %>% 
  dplyr::bind_rows(., .id = "source")
rio::export(comp_all_cell_mean_loc,"data/intermediate_output/vegetation/comp_all_cell_mean_loc.csv")



#Group 4 bioregion
comp_all_cell_bioregiongr4 <- comp_all_cell %>% 
  dplyr::filter(bioregion_name %in% c("Atlantic", "Boreal", "Continental", "Mediterranean"))

comp_all_cell_bioregiongr4_median <- comp_all_cell_bioregiongr4 %>% 
  dplyr::group_by(bin_age, source, bioregion_name) %>% 
  dplyr::summarise(tree_cover = median(tree_cover)) %>% 
  dplyr::ungroup()

comp_all_cell_bioregiongr4_mean <- comp_all_cell_bioregiongr4 %>% 
  dplyr::group_by(bin_age, source, bioregion_name) %>% 
  dplyr::summarise(tree_cover = mean(tree_cover)) %>% 
  dplyr::ungroup()

comp_all_cell_bioregiongr4_n <- comp_all_cell_bioregiongr4 %>% 
  dplyr::group_by(bin_age, source, bioregion_name) %>% 
  dplyr::summarise(number_records = n()) %>% 
  dplyr::ungroup()

comp_all_cell_bioregiongr4_median_ls <- comp_all_cell_bioregiongr4_median %>% 
  dplyr::select(bin_age, bioregion_name, tree_cover, source, bioregion_name) %>% 
  dplyr::mutate(grouping = paste0(source, "_", bioregion_name)) %>% 
  tidytable::group_split(grouping, .keep = TRUE, .named = TRUE)

comp_all_cell_bioregiongr4_mean_ls <- comp_all_cell_bioregiongr4_mean %>% 
  dplyr::select(bin_age, bioregion_name,tree_cover, source, bioregion_name) %>% 
  dplyr::mutate(grouping = paste0(source, "_", bioregion_name)) %>% 
  tidytable::group_split(grouping, .keep = TRUE, .named = TRUE)

comp_all_cell_bioregiongr4_median_loc <- lapply(comp_all_cell_bioregiongr4_median_ls, f_locfit_tree_average_bioregion) %>% 
  dplyr::bind_rows(., .id = "grouping") %>% 
  dplyr::mutate(bioregion_name = stringr::str_extract(grouping, "(?<=_).*"), source = stringr::str_extract(grouping, ".*(?=_)"))
comp_all_cell_bioregiongr4_mean_ls <- lapply(comp_all_cell_bioregiongr4_mean_ls, f_locfit_tree_average_bioregion) %>% 
  dplyr::bind_rows(., .id = "grouping") %>% 
  dplyr::mutate(bioregion_name = stringr::str_extract(grouping, "(?<=_).*"), source = stringr::str_extract(grouping, ".*(?=_)"))



##Modern model comparison
pollen_records2 <- rio::import("data/intermediate_output/vegetation/pollen_records2.csv") %>% 
  dplyr::arrange(ID_ENTITY)

pollen_records_xy3 <- pollen_records2 %>% #This is for the Copernicus tree cover CRS
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>%
  sf::st_transform(crs = 3035)

pollen_records_xy_buffer4 <- pollen_records_xy3 %>% 
  sf::st_buffer(., pollen_records_xy3$i0.75) %>% #based on basin_size and a median FSP from Githumbi
  dplyr::arrange(ID_ENTITY) %>% 
  dplyr::select(-c(i0.75)) %>%
  sf::st_as_sf(crs = 3035) #Ensure in sf


#Extract values for first bin Zanon
first_bin_zan <- terra::rast("./data/input/ZANON/forest_cover/forest_cover_0.grd")
terra::crs(first_bin_zan) <- "EPSG:4258"

first_bin_zan_sf <- terra::as.data.frame(first_bin_zan, xy = TRUE) %>% 
  sf::st_as_sf(coords = c("x", "y"), agr = "constant",crs = 4258) %>%
  sf::st_transform(crs = 3035) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1], y = sf::st_coordinates(.)[,2]) %>%  
  sf::st_drop_geometry() %>% 
  dplyr::select(x,y,z) %>% 
  dplyr::arrange(y) %>% 
  dplyr::mutate(ylag = y-(lag(y)))

first_bin_zan_3035_7km <- raster::rasterize(dplyr::select(first_bin_zan_sf,x,y), euro_map_3035_rast_7, field = first_bin_zan_sf$z, fun = mean) #approx resolution at 44º

first_bin_zan_extract <- exactextractr::exact_extract(first_bin_zan_3035_7km, pollen_records_xy_buffer4,'median', force_df = TRUE) %>% 
  dplyr::rename(zan_tree = 1) %>% 
  dplyr::mutate(ID_ENTITY = pollen_records_xy3$ID_ENTITY, .before = 1)
first_bin_zan_extract_ID <- terra::cellFromXY(first_bin_zan_3035_7km, sf::st_coordinates(pollen_records_xy3))  %>% 
  dplyr::as_tibble() %>% 
  dplyr::rename(cell_ID_zan = 1)   %>% 
  dplyr::mutate(ID_ENTITY = pollen_records_xy3$ID_ENTITY, .before = 1)


#Extract values for first bin Serge
first_bin_serge <- serge_tree_rast[[1]]
terra::crs(first_bin_serge) <- "EPSG:4258"

first_bin_serge_sf <- terra::as.data.frame(first_bin_serge, xy = TRUE) %>% 
  sf::st_as_sf(coords = c("x", "y"), agr = "constant",crs = 4258) %>%
  sf::st_transform(crs = 3035) %>% 
  tidyr::drop_na()  

first_bin_serge_3035_80km <- terra::rasterize(first_bin_serge_sf, euro_map_3035_rast_80, field = first_bin_serge_sf$serge_tree, fun = mean)

first_bin_serge_extract <- exactextractr::exact_extract(first_bin_serge_3035_80km, pollen_records_xy_buffer4,'median', force_df = TRUE) %>% 
  dplyr::rename(serge_tree = 1) %>% 
  dplyr::mutate(ID_ENTITY = pollen_records_xy3$ID_ENTITY, .before = 1)
first_bin_serge_extract_ID <- terra::cellFromXY(first_bin_serge_3035_80km, sf::st_coordinates(pollen_records_xy3)) %>% 
  dplyr::as_tibble() %>% 
  dplyr::rename(cell_ID_serge = 1) %>% 
  dplyr::mutate(ID_ENTITY = pollen_records_xy3$ID_ENTITY, .before = 1)

#Modern comparison
modern_comp <- ap_tree_cover %>% 
  dplyr::filter(tree>0) %>% 
  dplyr::left_join(first_bin_zan_extract, by = "ID_ENTITY") %>% 
  dplyr::left_join(first_bin_zan_extract_ID, by = "ID_ENTITY") %>% 
  dplyr::left_join(first_bin_serge_extract, by = "ID_ENTITY") %>% 
  dplyr::left_join(first_bin_serge_extract_ID, by = "ID_ENTITY") %>%
  dplyr::left_join(dplyr::select(fitted_tree_plot_qmap, ID_ENTITY, refitted), by = "ID_ENTITY")

modern_comp_zan <- modern_comp %>% 
  dplyr::select(zan_tree, tree, refitted) %>% 
  dplyr::mutate(tree = 100*tree) %>% 
  dplyr::mutate(refitted = 100*refitted) %>%
  tidyr::drop_na()
rio::export(modern_comp_zan,"data/intermediate_output/vegetation/modern_comp_zan.csv")

modern_comp_serge <- modern_comp %>% 
  dplyr::select(serge_tree, tree, refitted) %>% 
  dplyr::mutate(tree = 100*tree) %>% 
  dplyr::mutate(refitted = 100*refitted) %>%
  tidyr::drop_na()
rio::export(modern_comp_serge,"data/intermediate_output/vegetation/modern_comp_serge.csv")

cor(modern_comp_zan$zan_tree,modern_comp_zan$tree)
cor(modern_comp_serge$serge_tree,modern_comp_serge$tree)

cor(modern_comp_zan$zan_tree,modern_comp_zan$refitted)
cor(modern_comp_serge$serge_tree,modern_comp_serge$refitted)


  
#Modern comparison, single cell
modern_comp_zanID <- modern_comp %>% 
  dplyr::group_by(cell_ID_zan) %>% 
  dplyr::summarise(zan_tree = mean(zan_tree), tree = 100*mean(tree), refitted = 100*mean(refitted)) %>% 
  dplyr::ungroup() %>% 
  tidyr::drop_na()
rio::export(modern_comp_zanID, "data/intermediate_output/vegetation/modern_comp_zanID.csv")

modern_comp_sergeID <- modern_comp %>% 
  dplyr::group_by(cell_ID_serge) %>% 
  dplyr::summarise(serge_tree = mean(serge_tree), tree = 100*mean(tree), refitted = 100*mean(refitted)) %>% 
  dplyr::ungroup()  %>% 
  dplyr::filter(!is.na(serge_tree))
rio::export(modern_comp_sergeID, "data/intermediate_output/vegetation/modern_comp_sergeID.csv")

cor(modern_comp_zanID$zan_tree,modern_comp_zanID$tree)
cor(modern_comp_sergeID$serge_tree,modern_comp_sergeID$tree)

cor(modern_comp_zanID$zan_tree,modern_comp_zanID$refitted)
cor(modern_comp_sergeID$serge_tree,modern_comp_sergeID$refitted)



# ---------------------------------------------------------


# 5. (as 2,3,and partially 4 above, but where no shannon index included)
# ---------------------------------------------------------
# ---------------------------------------------------------

# 5.22. Load in data and quantile mapping adjustment model

#Input files
ap_tree_cover <- rio::import("data/intermediate_output/vegetation/ap_tree_cover.csv")
tree_model_ns <- readRDS("data/intermediate_output/vegetation/tree_beta_no_tree_shannon.rda") #modern tree cover model
boot_tree_model_ns <- readRDS("data/intermediate_output/vegetation/boot_tree_beta_no_tree_shannon.rda")
taxa_cat <- rio::import("data/input/taxa_cat.csv") %>% #import information regarding the taxa
  dplyr::arrange(taxon_name) 
taxa_cat_single <- taxa_cat %>% #for single clean_taxon_name
  dplyr::select(-taxon_name) %>% 
  dplyr::distinct()
loocv_predict_ns <- rio::import( "data/intermediate_output/vegetation/loocv_no_tree_shannon.csv")
boot_predict_ns <- readRDS( "data/intermediate_output/vegetation/boot_dataset_pred_df_no_tree_shannon.rda")

##Set map extents (as 1.Map_setup)
euro_extent_4258_ymin <- 34 
euro_extent_4258_ymax <- 73
euro_extent_4258_xmin <- -12
euro_extent_4258_xmax <- 45

#Load in necessary maps and rasters
euro_map_3035 <- sf::st_read("data/intermediate_output/vegetation/euro_map/3035/euro_map_3035.shp") #load map (1.Map_setup)
euro_map_3035_rast_50 <- raster::raster("data/intermediate_output/vegetation/euro_map/euro_map_3035_rast_50.tif") #(1.Map_setup)
euro_map_3035_rast_7 <- raster::raster("data/intermediate_output/vegetation/euro_map/euro_map_3035_rast_7.tif") #(1.Map_setup) - for use with Zanon data
euro_map_3035_rast_80 <- raster::raster("data/intermediate_output/vegetation/euro_map/euro_map_3035_rast_80.tif") #(1.Map_setup) - for use with Serge data

#Analyse relationships observed tree cover against LOOCV predictions 
fitted_tree_plot_ns <- loocv_predict_ns %>% 
  dplyr::select(ID_ENTITY, tree, prediction, E, latitude, longitude) %>% 
  dplyr::rename(observed = tree, predicted = prediction, difference = E) %>% #rename variables
  dplyr::mutate(difference = - difference) %>% #observation minus predictions
  dplyr::mutate(observed100 = observed*100, predicted100 = predicted*100) %>% #convert to precentage
  dplyr::mutate(observed_group = ggplot2::cut_width(observed100, width = 10, center = 5)) %>% 
  dplyr::arrange(observed_group)

rio::export(fitted_tree_plot_ns, "data/intermediate_output/vegetation/fitted_tree_plot_ns.csv")

fitted_tree_plot_nonaobs_ns <- dplyr::filter(fitted_tree_plot_ns, !is.na(observed)) #exclude records with  no observed cover
round(cor(fitted_tree_plot_nonaobs_ns$observed, fitted_tree_plot_nonaobs_ns$predicted),3) #correlation between obersevations and predictions
lm(fitted_tree_plot_nonaobs_ns$predicted~fitted_tree_plot_nonaobs_ns$observed)
round(max(fitted_tree_plot_ns$predicted, na.rm = TRUE),2) #investigate maximum predicted values
round(max(fitted_tree_plot_ns$observed, na.rm = TRUE),2) #investigate maximum observed values



#Empirical adjustment quantile mapping
tree_qmap_model_ssplin_ns <- qmap::fitQmapSSPLIN(fitted_tree_plot_nonaobs_ns$observed, fitted_tree_plot_nonaobs_ns$predicted, wet.day = TRUE, qstep = 0.001) #quantile fit using smoothing spline
fitted_tree_plot_qmap_ns <- fitted_tree_plot_nonaobs_ns %>% 
  dplyr::mutate(refitted = qmap::doQmapSSPLIN(fitted_tree_plot_nonaobs_ns$predicted, tree_qmap_model_ssplin_ns)) %>% #re-fit values
  dplyr::mutate(difference2 = refitted-observed) %>%   #refitted minus observed
  dplyr::mutate(refitted100 = refitted*100) %>% 
  dplyr::arrange(observed_group)

rio::export(fitted_tree_plot_qmap_ns, "data/intermediate_output/vegetation/fitted_tree_plot_qmap_ns.csv")

refitted_tree_plot_nonaobs_ns <- dplyr::filter(fitted_tree_plot_qmap_ns, !is.na(observed))
round(cor(refitted_tree_plot_nonaobs_ns$observed, refitted_tree_plot_nonaobs_ns$refitted),3) #correlation between observations and adjusted predictions

refitted_model_3035_ns <- fitted_tree_plot_qmap_ns %>%
  dplyr::select(longitude, latitude, difference2) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>%
  sf::st_transform(crs = 3035)

refitted_model_rast_ns <- terra::rasterize(refitted_model_3035_ns, euro_map_3035_rast_50, fun = mean)
refitted_model_df_ns <- terra::as.data.frame(refitted_model_rast_ns, xy = TRUE, na.rm = TRUE) %>% 
  dplyr::mutate(difference3 = ggplot2::cut_width(difference2, width = 0.4)) %>% 
  dplyr::mutate(Percentage_difference = dplyr::if_else(difference3 == "(-0.2,0.2]", "-20 to 20", difference3)) %>% 
  dplyr::mutate(Percentage_difference = dplyr::if_else(difference3 == "(-0.6,-0.2]", "-60 to -20", Percentage_difference)) %>%
  dplyr::mutate(Percentage_difference = dplyr::if_else(difference3 == "[-1,-0.6]", "-100 to -60", Percentage_difference)) %>%
  dplyr::mutate(Percentage_difference = dplyr::if_else(difference3 == "(0.2,0.6]", "20 to 60", Percentage_difference)) %>%
  dplyr::mutate(Percentage_difference = dplyr::if_else(difference3 == "(0.6,1]", "60 to 100", Percentage_difference)) 

refitted_model_df_ns$Percentage_difference <- factor(refitted_model_df_ns$Percentage_difference, levels = c("-100 to -60",
                                                                                                            "-60 to -20",
                                                                                                            "-20 to 20",
                                                                                                            "20 to 60",
                                                                                                            "60 to 100"))
rio::export(refitted_model_df_ns, "data/intermediate_output/vegetation/refitted_model_df_ns.csv")



#Analyse relationships observed tree cover against AP values
ap_tree_plot_ns <- loocv_predict_ns %>% 
  dplyr::select(ID_ENTITY, tree, ap_cover, latitude, longitude) %>% 
  dplyr::rename(observed = tree) %>% #rename variables
  dplyr::mutate(difference = observed - ap_cover) %>% 
  dplyr::mutate(difference = - difference) %>% # predictions minus observations
  dplyr::mutate(observed100 = observed*100, ap100 = ap_cover*100) %>% #convert to precentage
  dplyr::mutate(observed_group = ggplot2::cut_width(observed100, width = 10, center = 5)) %>% 
  dplyr::arrange(observed_group)

rio::export(ap_tree_plot_ns, "data/intermediate_output/vegetation/ap_tree_plot_ns.csv")

ap_tree_plot_nonaobs_ns <- dplyr::filter(ap_tree_plot_ns, !is.na(observed)) #exclude records with  no observed cover
round(cor(ap_tree_plot_nonaobs_ns$observed, ap_tree_plot_nonaobs_ns$ap_cover),3) #correlation between obersevations and predictions
round(max(ap_tree_plot_ns$ap_cover, na.rm = TRUE),2) #investigate maximum ap values
round(max(ap_tree_plot_ns$observed, na.rm = TRUE),2) #investigate maximum observed values


#Bootstraps
#Analyse relationships observed tree cover against LOOCV predictions
boot_fitted_tree_plot_ns <- boot_predict_ns %>%
  dplyr::mutate(E = tree - prediction) %>%
  dplyr::select(ID_ENTITY, tree, prediction, E, latitude, longitude, boot) %>%
  dplyr::rename(observed = tree, predicted = prediction, difference = E) %>% #rename variables
  dplyr::mutate(difference = - difference) %>% #observation minus predictions
  dplyr::mutate(observed100 = observed*100, predicted100 = predicted*100) %>% #convert to precentage
  dplyr::mutate(observed_group = ggplot2::cut_width(observed100, width = 10, center = 5)) %>%
  dplyr::arrange(boot,observed_group)

rio::export(boot_fitted_tree_plot_ns, "data/intermediate_output/vegetation/boot_fitted_tree_plot_ns.csv")
boot_fitted_tree_plot_nonaobs_ns <- dplyr::filter(boot_fitted_tree_plot_ns, !is.na(observed)) #exclude records with  no observed cover
boot_fitted_tree_plot_nonaobs_ls_ns <- boot_fitted_tree_plot_nonaobs_ns %>%
  dplyr::group_by(boot) %>%
  dplyr::group_split()

#Empirical adjustment quantile mapping
boot_tree_qmap_model_ssplin_ns <- lapply(boot_fitted_tree_plot_nonaobs_ls_ns, function(x){
  qmap::fitQmapSSPLIN(x$observed, x$predicted, wet.day = TRUE, qstep = 0.001) #quantile fit using smoothing spline
})



# ---------------------------------------------------------




# 5.3. Import fossil pollen data and reshape
# ---------------------------------------------------------
# ---------------------------------------------------------

#Load special EPD data 
age_model <- rio::import("data/input/SPECIAL-EPD/age_model.csv") %>% 
  dplyr::filter(median != "unknown") %>% 
  dplyr::mutate(median = as.numeric(median))
dates <- rio::import("data/input/SPECIAL-EPD/dates.csv") %>% 
  dplyr::rename(depth = `depth (cm)`)
entity <- rio::import("data/input/SPECIAL-EPD/metadata.csv")
pollen_count <- rio::import("data/input/SPECIAL-EPD/pollen_counts_amalgamated.csv")
sample <- rio::import("data/input/SPECIAL-EPD/samples.csv") %>% 
  dplyr::rename(depth = `depth (cm)`)

#Rename and reshape data
multiple_entity_site_check <- sample %>% #limit entities to one per site
  dplyr::left_join(dplyr::select(age_model, ID_SAMPLE, median), by = "ID_SAMPLE") %>%
  dplyr::filter(median <= 12000) %>%
  dplyr::group_by(ID_ENTITY) %>%
  dplyr::summarise(number_per_entity = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(entity, ID_ENTITY, ID_SITE), by = "ID_ENTITY")  %>%
  dplyr::group_by(ID_SITE)  %>%
  dplyr::slice_max(number_per_entity, with_ties = FALSE) %>% #if same number then just one allowed
  dplyr::ungroup()

lq_multiple_entity_site_check <- sample %>% #limit entities to one per site
  dplyr::left_join(dplyr::select(age_model, ID_SAMPLE, median, UNCERT_25), by = "ID_SAMPLE")  %>%
  dplyr::filter(UNCERT_25 != "unknown") %>% #for quartile analysis of age model
  dplyr::mutate(lowerq = median + as.numeric(UNCERT_25)) %>%
  dplyr::filter(lowerq >= -70) %>%  #To prevent dates into the future
  dplyr::filter(lowerq <= 12000) %>%
  dplyr::group_by(ID_ENTITY) %>%
  dplyr::summarise(number_per_entity = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(entity, ID_ENTITY, ID_SITE), by = "ID_ENTITY")  %>%
  dplyr::group_by(ID_SITE) %>%
  dplyr::slice_max(number_per_entity, with_ties = FALSE) %>% #if same number then just one allowed
  dplyr::ungroup()

uq_multiple_entity_site_check <- sample %>% #limit entities to one per site
  dplyr::left_join(dplyr::select(age_model, ID_SAMPLE, median, UNCERT_75), by = "ID_SAMPLE")  %>%
  dplyr::filter(UNCERT_75 != "unknown") %>% #for quartile analysis of age model
  dplyr::mutate(upperq = median - as.numeric(UNCERT_75)) %>%
  dplyr::filter(upperq >= -70) %>%  #To prevent dates into the future
  dplyr::filter(upperq <= 12000) %>%
  dplyr::group_by(ID_ENTITY) %>%
  dplyr::summarise(number_per_entity = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(entity, ID_ENTITY, ID_SITE), by = "ID_ENTITY")  %>%
  dplyr::group_by(ID_SITE) %>%
  dplyr::slice_max(number_per_entity, with_ties = FALSE) %>% #if same number then just one allowed
  dplyr::ungroup()

entity_epd <- entity %>% 
  dplyr::mutate(handle = entity_name) %>% 
  dplyr::filter(ID_ENTITY %in% c(multiple_entity_site_check$ID_ENTITY)) 

lq_entity_epd <- entity %>% 
  dplyr::mutate(handle = entity_name) %>% 
  dplyr::filter(ID_ENTITY %in% c(lq_multiple_entity_site_check$ID_ENTITY)) 

uq_entity_epd <- entity %>% 
  dplyr::mutate(handle = entity_name) %>% 
  dplyr::filter(ID_ENTITY %in% c(uq_multiple_entity_site_check$ID_ENTITY)) 

sample_epd <- sample
age_model_epd <- age_model
pollen_count_epd <- pollen_count

pollen_sample_ages <- sample_epd %>% #Depth and date info regarding samples
  dplyr::left_join(age_model_epd, by = "ID_SAMPLE") %>% 
  dplyr::left_join(entity_epd, by = "ID_ENTITY") %>% 
  dplyr::filter(!site_type %in% c("archaeological site", 
                                  "cave",
                                  "coastal, estuarine",
                                  "coastal, lagoon",
                                  "glacial",
                                  "fluvial",
                                  "marine"))  

lq_pollen_sample_ages <- sample_epd %>% #Depth and date info regarding samples
  dplyr::left_join(age_model_epd, by = "ID_SAMPLE")  %>% 
  dplyr::filter(UNCERT_25 != "unknown") %>% 
  dplyr::mutate(lowerq = median + as.numeric(UNCERT_25)) %>% 
  dplyr::filter(lowerq >= -70) %>%  #To prevent dates into the future
  dplyr::left_join(lq_entity_epd, by = "ID_ENTITY") %>% 
  dplyr::filter(!site_type %in% c("archaeological site", 
                                  "cave",
                                  "coastal, estuarine",
                                  "coastal, lagoon",
                                  "glacial",
                                  "fluvial",
                                  "marine"))   

uq_pollen_sample_ages <- sample_epd %>% #Depth and date info regarding samples
  dplyr::left_join(age_model_epd, by = "ID_SAMPLE")  %>% 
  dplyr::filter(UNCERT_75 != "unknown") %>% 
  dplyr::mutate(upperq = median - as.numeric(UNCERT_75)) %>% 
  dplyr::filter(upperq >= -70) %>%  #To prevent dates into the future
  dplyr::left_join(uq_entity_epd, by = "ID_ENTITY") %>% 
  dplyr::filter(!site_type %in% c("archaeological site", 
                                  "cave",
                                  "coastal, estuarine",
                                  "coastal, lagoon",
                                  "glacial",
                                  "fluvial",
                                  "marine"))  

pollen_sample_age_nona <- pollen_sample_ages %>% #filter to records with a median date
  dplyr::filter(!is.na(median))
total_epd_records <- length(unique(pollen_sample_age_nona$handle))

lq_pollen_sample_age_nona <- lq_pollen_sample_ages %>% 
  dplyr::filter(!is.na(lowerq))
lq_total_epd_records <- length(unique(lq_pollen_sample_age_nona$handle))

uq_pollen_sample_age_nona <- uq_pollen_sample_ages %>% 
  dplyr::filter(!is.na(upperq))
uq_total_epd_records <- length(unique(uq_pollen_sample_age_nona$handle))


#Select data
epd_counts_taxon_am2 <- pollen_count_epd %>% 
  tidyr::pivot_longer(!ID_SAMPLE, names_to = "taxon_name", values_to = "count") %>% 
  tidyr::drop_na() %>% 
  dplyr::left_join(dplyr::select(pollen_sample_ages, ID_SAMPLE, depth, median, handle), by = "ID_SAMPLE") %>% #Ensure required details added
  dplyr::rename(quantity = count) %>% #Rename for later
  dplyr::mutate(Notes = 0) %>% #Add for later combining
  dplyr::select(ID_SAMPLE, handle, depth, median, taxon_name, quantity, Notes) %>% 
  dplyr::filter(!is.na(median)) #To ensure all count info has a date

epd_pollen_counts <- epd_counts_taxon_am2 %>% #All of the pollen data, with dates
  dplyr::filter(median <= 14000) %>% #drop older data
  dplyr::filter(!is.na(quantity)) %>%  #drop na quantities
  dplyr::filter(quantity >= 0) %>% #get rid of negative quantities
  dplyr::arrange(handle, depth, taxon_name) %>% 
  dplyr::left_join(dplyr::select(entity_epd, handle, latitude, longitude, elevation), by = "handle") %>%  #add in coordinate and elevation info
  dplyr::filter(dplyr::between(latitude, euro_extent_4258_ymin, euro_extent_4258_ymax)) %>% #Limit to European range
  dplyr::filter(dplyr::between(longitude, euro_extent_4258_xmin, euro_extent_4258_xmax))

lq_epd_counts_taxon_am2 <- pollen_count_epd %>% 
  tidyr::pivot_longer(!ID_SAMPLE, names_to = "taxon_name", values_to = "count") %>%  #CSVs
  tidyr::drop_na() %>% #CSVs
  dplyr::left_join(dplyr::select(lq_pollen_sample_ages, ID_SAMPLE, depth, lowerq, handle), by = "ID_SAMPLE") %>% #Ensure required details added
  dplyr::rename(quantity = count) %>% #Rename for later
  dplyr::mutate(Notes = 0) %>% #Add for later combining
  dplyr::select(ID_SAMPLE, handle, depth, lowerq, taxon_name, quantity, Notes) %>% 
  dplyr::filter(!is.na(lowerq)) #To ensure all count info has a date

lq_epd_pollen_counts <- lq_epd_counts_taxon_am2 %>% #All of the pollen data, with dates
  dplyr::filter(lowerq <= 14000) %>% #drop older data
  dplyr::filter(!is.na(quantity)) %>%  #drop na quantities
  dplyr::filter(quantity >= 0) %>% #get rid of negative quantities
  dplyr::arrange(handle, depth, taxon_name) %>% 
  dplyr::left_join(dplyr::select(lq_entity_epd, handle, latitude, longitude, elevation), by = "handle") %>%  #add in coordinate and elevation info
  dplyr::filter(dplyr::between(latitude, euro_extent_4258_ymin, euro_extent_4258_ymax)) %>% #Limit to European range
  dplyr::filter(dplyr::between(longitude, euro_extent_4258_xmin, euro_extent_4258_xmax))

uq_epd_counts_taxon_am2 <- pollen_count_epd %>% 
  tidyr::pivot_longer(!ID_SAMPLE, names_to = "taxon_name", values_to = "count") %>%  #CSVs
  tidyr::drop_na() %>% #CSVs
  dplyr::left_join(dplyr::select(uq_pollen_sample_ages, ID_SAMPLE, depth, upperq, handle), by = "ID_SAMPLE") %>% #Ensure required details added
  dplyr::rename(quantity = count) %>% #Rename for later
  dplyr::mutate(Notes = 0) %>% #Add for later combining
  dplyr::select(ID_SAMPLE, handle, depth, upperq, taxon_name, quantity, Notes) %>% 
  dplyr::filter(!is.na(upperq))   #To ensure all count info has a date

uq_epd_pollen_counts <- uq_epd_counts_taxon_am2 %>% #All of the pollen data, with dates
  dplyr::filter(upperq <= 14000) %>% #drop older data
  dplyr::filter(!is.na(quantity)) %>%  #drop na quantities
  dplyr::filter(quantity >= 0) %>% #get rid of negative quantities
  dplyr::arrange(handle, depth, taxon_name) %>% 
  dplyr::left_join(dplyr::select(uq_entity_epd, handle, latitude, longitude, elevation), by = "handle") %>%  #add in coordinate and elevation info
  dplyr::filter(dplyr::between(latitude, euro_extent_4258_ymin, euro_extent_4258_ymax)) %>% #Limit to European range
  dplyr::filter(dplyr::between(longitude, euro_extent_4258_xmin, euro_extent_4258_xmax))

#Add taxon cleaner data and filter
epd_pollen_counts_cleaner <- epd_pollen_counts %>% 
  dplyr::left_join(taxa_cat, by = c("taxon_name")) %>% 
  dplyr::rename(date_value = median)
lq_epd_pollen_counts_cleaner <- lq_epd_pollen_counts %>% 
  dplyr::left_join(taxa_cat, by = c("taxon_name")) %>%
  dplyr::rename(date_value = lowerq)
uq_epd_pollen_counts_cleaner <- uq_epd_pollen_counts %>% 
  dplyr::left_join(taxa_cat, by = c("taxon_name")) %>%
  dplyr::rename(date_value = upperq)

epd_pollen_counts_cleaner_ls <- list(epd_pollen_counts_cleaner,lq_epd_pollen_counts_cleaner,uq_epd_pollen_counts_cleaner)
names(epd_pollen_counts_cleaner_ls) <- c("median", "lq", "uq")

epd_pollen_counts_tps_ls <- parallel::mclapply(epd_pollen_counts_cleaner_ls, function(x){
  x %>% 
    dplyr::filter(terrestrial_pollen_sum == "yes") %>% #just those within TPS
    dplyr::select(-taxon_name) %>% 
    dplyr::distinct() %>% #remove potential duplicate rows caused by clean_taxa_name being the same for some taxa_name
    dplyr::mutate(europe = dplyr::if_else(is.na(europe), "none", europe)) %>% 
    dplyr::filter(europe != "not native to Europe") %>%  #remove non native species
    dplyr::select(-europe) %>% 
    dplyr::select(ID_SAMPLE, handle, depth, date_value, latitude, longitude, elevation, clean_taxon_name, quantity) %>%
    dplyr::mutate(clean_taxon_name = stringr::str_replace_all(clean_taxon_name, " type", "")) %>% #to stop type being a seperate species
    dplyr::group_by(ID_SAMPLE, handle, depth, date_value, latitude, longitude, elevation, clean_taxon_name) %>% #to sum quantities over same species
    dplyr::summarise(quantity = sum(quantity)) %>% 
    dplyr::ungroup() 
}, mc.cores = 6)


##Calculate indices
epd_pollen_wider_ls <- parallel::mclapply(epd_pollen_counts_tps_ls, function(x){
  x %>% 
    dplyr::arrange(ID_SAMPLE, clean_taxon_name) %>% 
    dplyr::select(ID_SAMPLE, clean_taxon_name, quantity) %>% 
    tidyr::pivot_wider(names_from = clean_taxon_name, values_from = quantity, values_fill = 0) 
},mc.cores = 6)

#Pollen percentages
epd_pollen_wider_per_ls <- parallel::mclapply(epd_pollen_wider_ls, function(x){
  x %>% 
    smpds::normalise_taxa(cols = 1)  #calculate percentages, ignoring first column
},mc.cores = 6)


#Shannon index
# epd_shannon_ls <- parallel::mclapply(epd_pollen_wider_per_ls, function(x){
#   x %>% 
#     dplyr::select(-ID_SAMPLE) %>% 
#     as.matrix() %>% 
#     vegan::diversity(., index = "shannon") %>% 
#     dplyr::as_tibble() %>%   
#     dplyr::mutate(ID_SAMPLE = x$ID_SAMPLE, .before = 1) %>% 
#     dplyr::rename(shannon = value)
# }, mc.cores = 6)
# epd_tree_shannon_ls <- parallel::mclapply(epd_pollen_wider_per_ls, function(x){
#   x %>% 
#     tidyr::pivot_longer(!c(ID_SAMPLE), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% #re-shape 
#       dplyr::left_join(dplyr::select(taxa_cat_single, clean_taxon_name, ap_sp_hp), by = "clean_taxon_name") %>%
#       dplyr::filter(ap_sp_hp == "AP") %>%
#       dplyr::select(ID_SAMPLE, clean_taxon_name, percentage_cover) %>%
#       tidyr::pivot_wider(names_from = clean_taxon_name, values_from = percentage_cover) %>%
#       dplyr::select(-ID_SAMPLE) %>%
#       as.matrix() %>%
#       vegan::diversity(., index = "shannon") %>%
#       dplyr::as_tibble() %>%
#       dplyr::mutate(ID_SAMPLE = x$ID_SAMPLE, .before = 1) %>%
#       dplyr::rename(tree_shannon = value)
# }, mc.cores = 6)


#Hills N2 filter
epd_hillsN2_ls <- parallel::mclapply(epd_pollen_wider_ls, function(x){
  x %>% 
    dplyr::select(-ID_SAMPLE) %>% 
    analogue::n2(., "sites") %>%
    dplyr::as_tibble() %>%
    dplyr::rename(HillsN2 = value) %>% 
    dplyr::mutate(ID_SAMPLE = x$ID_SAMPLE, .before = 1)
},mc.cores = 6)


#AP percentages
epd_pollen_counts_APinfo_ls <- parallel::mclapply(epd_pollen_wider_per_ls, function(x){
  x %>% 
    tidyr::pivot_longer(!c(ID_SAMPLE), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% 
    dplyr::left_join(taxa_cat_single, by = "clean_taxon_name") 
},mc.cores = 6)

epd_pollen_ap_per_ls <- parallel::mclapply(epd_pollen_counts_APinfo_ls, function(x){
  x %>% 
    dplyr::filter(ap_sp_hp == "AP") %>% #Just AP
    dplyr::group_by(ID_SAMPLE) %>% 
    dplyr::summarise(ap_cover = sum(percentage_cover)) %>% 
    dplyr::ungroup()
},mc.cores = 6)


#Needleshare 
epd_pollen_ap_needleshare_per_ls1 <- lapply(epd_pollen_counts_APinfo_ls,function(x){
  x %>%
    dplyr::filter(ap_sp_hp == "AP") %>%
    dplyr::mutate(needle = dplyr::if_else(is.na(ap_needle_broad), "unknown", ap_needle_broad)) %>%   #need to remove NA
    dplyr::filter(needle == "needle") %>%
    dplyr::group_by(ID_SAMPLE) %>%
    dplyr::summarise(ap_cover_needle = sum(percentage_cover)) %>%
    dplyr::ungroup()
})
epd_pollen_ap_needleshare_per_ls2 <- purrr::map2(epd_pollen_ap_needleshare_per_ls1, epd_pollen_ap_per_ls, dplyr::left_join)
epd_pollen_ap_needleshare_per_ls3 <- lapply(epd_pollen_ap_needleshare_per_ls2, function(x){
  x %>% 
    dplyr::mutate(needle_share = ap_cover_needle/ap_cover)  %>%
    dplyr::mutate(needle_share = tidyr::replace_na(needle_share, 0))
})


#SP cover
epd_pollen_counts_SPinfo_ls <- parallel::mclapply(epd_pollen_wider_per_ls, function(x){
  x %>% 
    tidyr::pivot_longer(!c(ID_SAMPLE), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% 
    dplyr::left_join(taxa_cat_single, by = "clean_taxon_name") 
}, mc.cores = 6)

epd_pollen_sp_per_ls <- parallel::mclapply(epd_pollen_counts_SPinfo_ls, function(x){
  x %>% 
    dplyr::filter(ap_sp_hp == "SP") %>% #Just SP
    dplyr::group_by(ID_SAMPLE) %>% 
    dplyr::summarise(sp_cover = sum(percentage_cover)) %>% 
    dplyr::ungroup()
}, mc.cores = 6)


##Input for downcore
epd_pollen_counts_tps <- epd_pollen_counts_tps_ls$median
# epd_shannon <- epd_shannon_ls$median
# epd_tree_shannon <- epd_tree_shannon_ls$median
epd_hillsN2 <- epd_hillsN2_ls$median
epd_pollen_ap_per <- epd_pollen_ap_per_ls$median
epd_pollen_sp_per <- epd_pollen_sp_per_ls$median
epd_pollen_ap_needleshare_per <- epd_pollen_ap_needleshare_per_ls3$median

lq_epd_pollen_counts_tps <- epd_pollen_counts_tps_ls$lq
# lq_epd_shannon <- epd_shannon_ls$lq
# lq_epd_tree_shannon <- epd_tree_shannon_ls$lq
lq_epd_hillsN2 <- epd_hillsN2_ls$lq
lq_epd_pollen_ap_per <- epd_pollen_ap_per_ls$lq
lq_epd_pollen_sp_per <- epd_pollen_sp_per_ls$lq
lq_epd_pollen_ap_needleshare_per <- epd_pollen_ap_needleshare_per_ls3$lq

uq_epd_pollen_counts_tps <- epd_pollen_counts_tps_ls$uq
# uq_epd_shannon <- epd_shannon_ls$uq
# uq_epd_tree_shannon <- epd_tree_shannon_ls$uq
uq_epd_hillsN2 <- epd_hillsN2_ls$uq
uq_epd_pollen_ap_per <- epd_pollen_ap_per_ls$uq
uq_epd_pollen_sp_per <- epd_pollen_sp_per_ls$uq
uq_epd_pollen_ap_needleshare_per <- epd_pollen_ap_needleshare_per_ls3$uq


epd_downcore_input_ns <- epd_pollen_counts_tps %>% 
  dplyr::select(ID_SAMPLE, elevation) %>%
  dplyr::distinct() %>% 
  # dplyr::left_join(epd_shannon, by = "ID_SAMPLE") %>%
  # dplyr::left_join(epd_tree_shannon, by = "ID_SAMPLE") %>%
  dplyr::left_join(epd_hillsN2, by = "ID_SAMPLE") %>% 
  dplyr::left_join(epd_pollen_ap_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(ap_cover = ap_cover / 100) %>% 
  dplyr::left_join(epd_pollen_sp_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(sp_cover = sp_cover / 100) %>% 
  dplyr::left_join(dplyr::select(epd_pollen_ap_needleshare_per, ID_SAMPLE, needle_share), by = "ID_SAMPLE") %>% 
  dplyr::filter(HillsN2 >= 2) %>% #Hills filter
  dplyr::select(-HillsN2) %>% 
  dplyr::left_join(dplyr::select(pollen_sample_ages, ID_SAMPLE, site_type), by = "ID_SAMPLE") %>% 
  dplyr::mutate(site_model = dplyr::if_else(site_type == "lake", 1, 0)) %>% 
  dplyr::select(-site_type)

rio::export(epd_downcore_input_ns, "data/intermediate_output/vegetation/epd_downcore_input_ns.csv")

lq_epd_downcore_input_ns <- lq_epd_pollen_counts_tps %>% 
  dplyr::select(ID_SAMPLE, elevation) %>%
  dplyr::distinct() %>% 
  # dplyr::left_join(lq_epd_shannon, by = "ID_SAMPLE") %>%
  # dplyr::left_join(lq_epd_tree_shannon, by = "ID_SAMPLE") %>%
  dplyr::left_join(lq_epd_hillsN2, by = "ID_SAMPLE") %>% 
  dplyr::left_join(lq_epd_pollen_ap_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(ap_cover = ap_cover / 100) %>% 
  dplyr::left_join(lq_epd_pollen_sp_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(sp_cover = sp_cover / 100) %>% 
  dplyr::left_join(dplyr::select(lq_epd_pollen_ap_needleshare_per, ID_SAMPLE, needle_share), by = "ID_SAMPLE") %>% 
  dplyr::filter(HillsN2 >= 2) %>% #Hills filter
  dplyr::select(-HillsN2) %>% 
  dplyr::left_join(dplyr::select(lq_pollen_sample_ages, ID_SAMPLE, site_type), by = "ID_SAMPLE") %>% 
  dplyr::mutate(site_model = dplyr::if_else(site_type == "lake", 1, 0)) %>% 
  dplyr::select(-site_type)

rio::export(lq_epd_downcore_input_ns, "data/intermediate_output/vegetation/lq_epd_downcore_input_ns.csv")

uq_epd_downcore_input_ns <- uq_epd_pollen_counts_tps %>% 
  dplyr::select(ID_SAMPLE, elevation) %>%
  dplyr::distinct() %>% 
  # dplyr::left_join(uq_epd_shannon, by = "ID_SAMPLE") %>%
  # dplyr::left_join(uq_epd_tree_shannon, by = "ID_SAMPLE") %>%
  dplyr::left_join(uq_epd_hillsN2, by = "ID_SAMPLE") %>% 
  dplyr::left_join(uq_epd_pollen_ap_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(ap_cover = ap_cover / 100) %>% 
  dplyr::left_join(uq_epd_pollen_sp_per, by = "ID_SAMPLE") %>% 
  dplyr::mutate(sp_cover = sp_cover / 100) %>% 
  dplyr::left_join(dplyr::select(uq_epd_pollen_ap_needleshare_per, ID_SAMPLE, needle_share), by = "ID_SAMPLE") %>% 
  dplyr::filter(HillsN2 >= 2) %>% #Hills filter
  dplyr::select(-HillsN2) %>% 
  dplyr::left_join(dplyr::select(uq_pollen_sample_ages, ID_SAMPLE, site_type), by = "ID_SAMPLE") %>% 
  dplyr::mutate(site_model = dplyr::if_else(site_type == "lake", 1, 0)) %>% 
  dplyr::select(-site_type)

rio::export(uq_epd_downcore_input_ns, "data/intermediate_output/vegetation/uq_epd_downcore_input_ns.csv")

epd_downcore_input_records_ns <- epd_downcore_input_ns %>% 
  dplyr::left_join(dplyr::select(pollen_sample_ages, ID_SAMPLE, ID_ENTITY, latitude, longitude, site_type), by = "ID_SAMPLE") %>% 
  dplyr::filter(site_type %in% c("lake", "terrestrial bog/mire/fen", "terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, bog/mire/fen", "terrestrial, marsh")) %>% 
  dplyr::filter(elevation < 1000) %>% #Limit to records below 1000m
  dplyr::select(ID_ENTITY, latitude, longitude) %>% 
  dplyr::distinct()

pollen_records2 <- rio::import("data/intermediate_output/vegetation/pollen_records2.csv")
pollen_records2_filt <- pollen_records2 %>% 
  dplyr::filter(ID_ENTITY %in% c(ap_tree_cover$ID_ENTITY))

epd_downcore_input_records_vect_ns <- epd_downcore_input_records_ns %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>% 
  sf::st_transform(crs = 3035) %>% 
  terra::vect()
epd_downcore_input_records_buf5k_ns <- epd_downcore_input_records_ns %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>% 
  sf::st_transform(crs = 3035) %>% 
  sf::st_buffer(5000)
epd_downcore_input_records_buf_sourcemedian_ns <- epd_downcore_input_records_ns %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>% 
  sf::st_transform(crs = 3035) %>% 
  sf::st_buffer(median(pollen_records2_filt$i0.75))

cop_tree_cover_masked_3035 <-  terra::rast("data/intermediate_output/vegetation/copernicus/cop_masked_tree_cover_3035_100m.tif")
epd_downcore_input_observedtree_ns <- terra::extract(cop_tree_cover_masked_3035, epd_downcore_input_records_vect_ns)
epd_downcore_input_observedtree_5k_ns <- exactextractr::exact_extract(cop_tree_cover_masked_3035, epd_downcore_input_records_buf5k_ns, 'mean', force_df = TRUE)
epd_downcore_input_observedtree_sourcemedian_ns <- exactextractr::exact_extract(cop_tree_cover_masked_3035, epd_downcore_input_records_buf_sourcemedian_ns, 'mean', force_df = TRUE)
median(epd_downcore_input_observedtree_ns$mean, na.rm = TRUE)
median(epd_downcore_input_observedtree_5k_ns$mean, na.rm = TRUE)
median(epd_downcore_input_observedtree_sourcemedian_ns$mean, na.rm = TRUE)

# ---------------------------------------------------------




# 5.4. Reconstructions
# ---------------------------------------------------------
# ---------------------------------------------------------
## Modelled fit
epd_downcore_input_ls_ns <- list(rio::import("data/intermediate_output/vegetation/epd_downcore_input_ns.csv"),
                                 rio::import("data/intermediate_output/vegetation/lq_epd_downcore_input_ns.csv"),
                                 rio::import("data/intermediate_output/vegetation/uq_epd_downcore_input_ns.csv"))
names(epd_downcore_input_ls_ns) <- c("median", "lq", "uq")
pred_EPD_tree_ls_ns <- lapply(epd_downcore_input_ls_ns, function(x){
  betareg::predict(tree_model_ns, newdata = dplyr::select(x, ap_cover, elevation, needle_share, site_model, sp_cover),type="response") #modelled
})

pred_EPD_tree_refit_ls_ns <- lapply(pred_EPD_tree_ls_ns, function(x){
  qmap::doQmapSSPLIN(as.numeric(x), tree_qmap_model_ssplin_ns)
})

recon_tree_ns <- cbind(epd_downcore_input_ls_ns$median, pred_EPD_tree_refit_ls_ns$median)  %>%
  dplyr::rename(tree_cover = dplyr::last_col()) %>%
  dplyr::mutate(tree_cover = dplyr::if_else(tree_cover>1, 1, tree_cover)) %>% #Because otherwise we may have re-fitted values greater than 100%
  dplyr::mutate(tree_cover = 100*tree_cover) %>% 
  dplyr::left_join(dplyr::select(pollen_sample_ages, ID_SAMPLE, entity_name, latitude, longitude, median, site_type), by = "ID_SAMPLE")  %>% 
  dplyr::filter(site_type %in% c("lake", "terrestrial bog/mire/fen", "terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, bog/mire/fen", "terrestrial, marsh")) %>% 
  dplyr::filter(elevation < 1000) #Limit to records below 1000m
rio::export(recon_tree_ns, "data/intermediate_output/vegetation/recon_tree_ns.csv")

lq_recon_tree_ns <- cbind(epd_downcore_input_ls_ns$lq, pred_EPD_tree_refit_ls_ns$lq)  %>%
  dplyr::rename(tree_cover = dplyr::last_col()) %>%
  dplyr::mutate(tree_cover = dplyr::if_else(tree_cover>1, 1, tree_cover)) %>% #Because otherwise we may have re-fitted values greater than 100%
  dplyr::mutate(tree_cover = 100*tree_cover) %>% 
  dplyr::left_join(dplyr::select(lq_pollen_sample_ages, ID_SAMPLE, entity_name, latitude, longitude, lowerq, site_type), by = "ID_SAMPLE")  %>% 
  dplyr::filter(site_type %in% c("lake", "terrestrial bog/mire/fen", "terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, bog/mire/fen", "terrestrial, marsh")) %>% 
  dplyr::filter(elevation < 1000) #Limit to records below 1000m
rio::export(lq_recon_tree_ns, "data/intermediate_output/vegetation/lq_recon_tree_ns.csv")

uq_recon_tree_ns <- cbind(epd_downcore_input_ls_ns$uq, pred_EPD_tree_refit_ls_ns$uq)  %>%
  dplyr::rename(tree_cover = dplyr::last_col()) %>%
  dplyr::mutate(tree_cover = dplyr::if_else(tree_cover>1, 1, tree_cover)) %>% #Because otherwise we may have re-fitted values greater than 100%
  dplyr::mutate(tree_cover = 100*tree_cover) %>% 
  dplyr::left_join(dplyr::select(uq_pollen_sample_ages, ID_SAMPLE, entity_name, latitude, longitude, upperq, site_type), by = "ID_SAMPLE")  %>% 
  dplyr::filter(site_type %in% c("lake", "terrestrial bog/mire/fen", "terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, bog/mire/fen", "terrestrial, marsh")) %>% 
  dplyr::filter(elevation < 1000) #Limit to records below 1000m
rio::export(uq_recon_tree_ns, "data/intermediate_output/vegetation/uq_recon_tree_ns.csv")

ap_cover_tree_ns <- recon_tree_ns %>% 
  dplyr::mutate(tree_cover = ap_cover*100) #convert to ap cover subsequent analysis
rio::export(ap_cover_tree_ns, "data/intermediate_output/vegetation/ap_cover_tree_ns.csv")

#Bin data into 200yrs bins
recon_tree_ls_ns <- list(recon_tree_ns, lq_recon_tree_ns, uq_recon_tree_ns, ap_cover_tree_ns)
names(recon_tree_ls_ns) <- c("median", "lq", "uq", "ap")

recon_tree_binned_200_ls_ns <- lapply(recon_tree_ls_ns, function(x){
  x %>% 
    dplyr::rename(date_value = 11) %>% 
    dplyr::mutate(bin = ggplot2::cut_width(date_value, width = 200, center = 0, labels = F)) %>%    #place in 200 year bins
    dplyr::mutate(bin_age = (bin * 200) - 200) %>%
    dplyr::filter(bin_age <=14000)  %>% 
    dplyr::group_by(entity_name, bin_age) %>% 
    dplyr::summarise(tree_cover = mean(tree_cover), number = n()) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(dplyr::select(pollen_sample_ages, entity_name, latitude, longitude), by = "entity_name") %>% 
    dplyr::distinct()
})

bin_200_mean_recon_tree_ls_ns <- lapply(recon_tree_binned_200_ls_ns, function(x){ #Calculate mean through time
  x %>% 
    dplyr::group_by(bin_age) %>% 
    dplyr::summarise(tree_cover_mean = mean(tree_cover), tree_cover_sd = sd(tree_cover)) %>% 
    dplyr::ungroup()
}) 

bin_200_median_recon_tree_ls_ns <- lapply(recon_tree_binned_200_ls_ns, function(x){ #Calculate median through time
  x %>% 
    dplyr::group_by(bin_age) %>% 
    dplyr::summarise(tree_cover_median = median(tree_cover), tree_cover_lower = quantile(tree_cover,probs = c(0.25)),tree_cover_higher = quantile(tree_cover,probs = c(0.75))) %>% 
    dplyr::ungroup()
}) 

bin_200_max_recon_tree_ls_ns <- lapply(recon_tree_binned_200_ls_ns, function(x){ #Calculate max through time
  x %>% 
    dplyr::group_by(bin_age) %>% 
    dplyr::summarise(tree_cover_max = max(tree_cover)) %>% 
    dplyr::ungroup()
}) 

bin_200_n_recon_tree_ls_ns <- lapply(recon_tree_binned_200_ls_ns, function(x){ #Calculate number through time
  x %>% 
    dplyr::select(-tree_cover, -number) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(bin_age) %>% 
    dplyr::summarise(number_records = n()) %>% 
    dplyr::ungroup()
}) 

rio::export(recon_tree_binned_200_ls_ns$median, "data/intermediate_output/vegetation/recon_tree_binned_200_ns.csv")
rio::export(bin_200_mean_recon_tree_ls_ns$median, "data/intermediate_output/vegetation/bin_200_mean_recon_tree_ns.csv")
rio::export(bin_200_median_recon_tree_ls_ns$median, "data/intermediate_output/vegetation/bin_200_median_recon_tree_ns.csv")
rio::export(bin_200_max_recon_tree_ls_ns$median, "data/intermediate_output/vegetation/bin_200_max_recon_tree_ns.csv")
rio::export(bin_200_n_recon_tree_ls_ns$median, "data/intermediate_output/vegetation/bin_200_n_recon_tree_ns.csv")

rio::export(recon_tree_binned_200_ls_ns$lq, "data/intermediate_output/vegetation/lq_recon_tree_binned_200_ns.csv")
rio::export(bin_200_mean_recon_tree_ls_ns$lq, "data/intermediate_output/vegetation/lq_bin_200_mean_recon_tree_ns.csv")
rio::export(bin_200_median_recon_tree_ls_ns$lq, "data/intermediate_output/vegetation/lq_bin_200_median_recon_tree_ns.csv")
rio::export(bin_200_max_recon_tree_ls_ns$lq, "data/intermediate_output/vegetation/lq_bin_200_max_recon_tree_ns.csv")
rio::export(bin_200_n_recon_tree_ls_ns$lq, "data/intermediate_output/vegetation/lq_bin_200_n_recon_tree_ns.csv")

rio::export(recon_tree_binned_200_ls_ns$uq, "data/intermediate_output/vegetation/uq_recon_tree_binned_200_ns.csv")
rio::export(bin_200_mean_recon_tree_ls_ns$uq, "data/intermediate_output/vegetation/uq_bin_200_mean_recon_tree_ns.csv")
rio::export(bin_200_median_recon_tree_ls_ns$uq, "data/intermediate_output/vegetation/uq_bin_200_median_recon_tree_ns.csv")
rio::export(bin_200_max_recon_tree_ls_ns$uq, "data/intermediate_output/vegetation/uq_bin_200_max_recon_tree_ns.csv")
rio::export(bin_200_n_recon_tree_ls_ns$uq, "data/intermediate_output/vegetation/uq_bin_200_n_recon_tree_ns.csv")

rio::export(recon_tree_binned_200_ls_ns$ap, "data/intermediate_output/vegetation/ap_tree_binned_200_ns.csv")
rio::export(bin_200_mean_recon_tree_ls_ns$ap, "data/intermediate_output/vegetation/bin_200_mean_ap_tree_ns.csv")
rio::export(bin_200_median_recon_tree_ls_ns$ap, "data/intermediate_output/vegetation/bin_200_median_ap_tree_ns.csv")
rio::export(bin_200_max_recon_tree_ls_ns$ap, "data/intermediate_output/vegetation/bin_200_max_ap_tree_ns.csv")
rio::export(bin_200_n_recon_tree_ls_ns$ap, "data/intermediate_output/vegetation/bin_200_n_ap_tree_ns.csv")


##Reconstruction analysis
recon_tree_ns <- rio::import("data/intermediate_output/vegetation/recon_tree_ns.csv")
recon_tree_binned_200_ns <- recon_tree_binned_200_ls_ns$median
bin_200_mean_recon_tree_ns <- bin_200_mean_recon_tree_ls_ns$median
bin_200_median_recon_tree_ns <- bin_200_median_recon_tree_ls_ns$median
bin_200_max_recon_tree_ns <- bin_200_max_recon_tree_ls_ns$median
bin_200_n_recon_tree_ns <- bin_200_n_recon_tree_ls_ns$median

binned_200yr_reconstructions_ns <- recon_tree_binned_200_ns %>% 
  dplyr::rename(bin_centre = bin_age, number_samples = number) %>% 
  dplyr::select(entity_name, longitude, latitude, bin_centre, number_samples, tree_cover)
rio::export(binned_200yr_reconstructions_ns, "figs/binned_200yr_reconstructions_ns.csv")

#Bootstraps medians
nrecords_ns <- length(unique(recon_tree_binned_200_ns$entity_name)) #Number of records
recon_tree_binned_200_records_ns <- recon_tree_binned_200_ns %>% #Add row number of each entity
  dplyr::select(entity_name) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(entity_row = dplyr::row_number())
recon_tree_binned_200_num_ns <- recon_tree_binned_200_ns %>% 
  dplyr::left_join(recon_tree_binned_200_records_ns, by = "entity_name")

boot_ls_ns <- list()
set.seed(42)
nreps <- 1000 #number of reps
for (i in 1:nreps){ #generate list of boostrapped resampled dfs
  print(i)
  record_boot_ns <- sample(seq(1:nrecords_ns), nrecords_ns, replace = TRUE) #select records, with replacement
  bin_200_recon_tree_boot_ns <- dplyr::tibble(entity_row = record_boot_ns) %>% 
    dplyr::left_join(recon_tree_binned_200_num_ns, by = "entity_row", relationship = "many-to-many") %>% #add info by entity
    dplyr::group_by(bin_age) %>% 
    dplyr::summarise(tree_cover_mean = mean(tree_cover), tree_cover_median = median(tree_cover),tree_cover_max = max(tree_cover) ) %>% #calculate bin mean, median and max
    dplyr::ungroup() 
  boot_ls_ns[[i]] <- bin_200_recon_tree_boot_ns #add to list
}
bin_200_recon_tree_boots_ns <- dplyr::bind_rows(boot_ls_ns, .id = "column_label") #reduce to single df
rio::export(bin_200_recon_tree_boots_ns, "data/intermediate_output/vegetation/bin_200_recon_tree_boots_ns.csv")

bin_200_recon_tree_boots_5_95_ns <- bin_200_recon_tree_boots_ns %>% #Calculate 95% CIs
  dplyr::group_by(bin_age) %>% 
  dplyr::summarise(tree_cover_lower_mean = quantile(tree_cover_mean,probs = c(0.025)),
                   tree_cover_higher_mean = quantile(tree_cover_mean,probs = c(0.975)),
                   tree_cover_lower_median = quantile(tree_cover_median,probs = c(0.025)),
                   tree_cover_higher_median = quantile(tree_cover_median,probs = c(0.975)),
                   tree_cover_lower_max = quantile(tree_cover_max,probs = c(0.025)),
                   tree_cover_higher_max = quantile(tree_cover_max,probs = c(0.975)),
  ) %>% 
  dplyr::ungroup()  %>% 
  tidyr::pivot_longer(!bin_age, names_to = "lower_upper", values_to = "tree_cover_5_95") %>%  
  dplyr::mutate(column_label = dplyr::if_else(stringr::str_detect(lower_upper, "tree_cover_lower"), nreps+2, nreps+3)) #need to add for plotting
rio::export(bin_200_recon_tree_boots_5_95_ns, "data/intermediate_output/vegetation/bin_200_recon_tree_boots_5_95_ns.csv")

bin_200_mean_recon_tree_plot_ns <- bin_200_mean_recon_tree_ns %>% 
  dplyr::mutate(column_label = nreps+1)
rio::export(bin_200_mean_recon_tree_plot_ns, "data/intermediate_output/vegetation/bin_200_mean_recon_tree_plot_ns.csv")

bin_200_median_recon_tree_plot_ns <- bin_200_median_recon_tree_ns %>% 
  dplyr::mutate(column_label = nreps+1) #need to include within plot
rio::export(bin_200_median_recon_tree_plot_ns, "data/intermediate_output/vegetation/bin_200_median_recon_tree_plot_ns.csv")

bin_200_max_recon_tree_plot_ns <- bin_200_max_recon_tree_ns %>% 
  dplyr::mutate(column_label = nreps+1)
rio::export(bin_200_max_recon_tree_plot_ns, "data/intermediate_output/vegetation/bin_200_max_recon_tree_plot_ns.csv")

bin_200_recon_tree_boots_5_95_mean_ns <- bin_200_recon_tree_boots_5_95_ns %>% 
  dplyr::filter(stringr::str_detect(lower_upper, "mean"))
rio::export(bin_200_recon_tree_boots_5_95_mean_ns, "data/intermediate_output/vegetation/bin_200_recon_tree_boots_5_95_mean_ns.csv")

bin_200_recon_tree_boots_5_95_median_ns <- bin_200_recon_tree_boots_5_95_ns %>% 
  dplyr::filter(stringr::str_detect(lower_upper, "median"))
rio::export(bin_200_recon_tree_boots_5_95_median_ns, "data/intermediate_output/vegetation/bin_200_recon_tree_boots_5_95_median_ns.csv")

bin_200_recon_tree_boots_5_95_max_ns <- bin_200_recon_tree_boots_5_95_ns %>% 
  dplyr::filter(stringr::str_detect(lower_upper, "max"))
rio::export(bin_200_recon_tree_boots_5_95_max_ns, "data/intermediate_output/vegetation/bin_200_recon_tree_boots_5_95_max_ns.csv")


