#Sweeney, Harrison and Vander Linden 2024. 
#European forest cover during the Holocene reconstructed from pollen records
# ---------------------------------------------------------

#There are six seperate scripts associated with this research:
# 1. Map_setup.R
# 2. Modern_pollen.R 
# 3. Model_import_data.R (this script)
# 4. Tree_model.R
# 5. Tree_recon.R
# 6. Figures.R

# ---------------------------------------------------------

############ Model import data ############ 

# This script generates the data required to build the model of modern
# tree cover


# ---------------------------------------------------------
# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Load in modern pollen data and calculate source area
# 3. Shannon index calculation
# 4. Extract tree cover values based on pollen site locations

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
# dir.create("data/intermediate_output/vegetation/SMPDS")
# dir.create("data/intermediate_output/vegetation/raster_tree")
# dir.create("figs")
# dir.create("figs/supplement")
# dir.create("figs/tree_cover_rasters")


# ---------------------------------------------------------




# 2. Load in modern pollen data and calculate source area
# ---------------------------------------------------------
# ---------------------------------------------------------

#Pollen data
pollen <- rio::import("data/intermediate_output/vegetation/SMPDS/SMPDS_am2_per.csv") #all of the pollen data by record, just TPS and native

taxa_names <- pollen %>% #taxa names for pollen records selected
  dplyr::select(Abies: ncol(.)) %>% 
  colnames(.)

#Load in taxa information
taxa_cat <- rio::import("data/input/taxa_cat.csv") %>% #import information regarding the taxa
  dplyr::arrange(taxon_name) 
taxa_cat_single <- taxa_cat %>% #for single clean_taxon_name
  dplyr::select(-taxon_name) %>% 
  dplyr::distinct()

#Recalculate values / percentages with one entity per site/per location
pollen_basin <- pollen %>% #records with basin information 
  dplyr::select(ID_ENTITY, basin_size) %>%
  dplyr::distinct() %>% #single row per entity
  dplyr::filter(!(basin_size %in% c("not applicable", "not known", "not recorded"))) %>%  #leave out records without basin info
  dplyr::mutate(basin_size = dplyr::if_else(basin_size == "large (50.1-500 km2)", "275", basin_size)) %>% #standard values for categories
  dplyr::mutate(basin_size = dplyr::if_else(basin_size == "medium (1.1-50 km2)", "25.5", basin_size)) %>% #standard values for categories
  dplyr::mutate(basin_size = dplyr::if_else(basin_size == "small (0.01-1 km2)", "0.5", basin_size)) %>% #standard values for categories
  dplyr::mutate(basin_size = as.numeric(basin_size)) %>% 
  dplyr::filter(basin_size != 39000) %>%  #this is the Sea of Azov
  dplyr::filter(basin_size != 18000) #area refers to entire hydrological basin in SW Turkey

pollen_single <- pollen %>% #single record per site meta info
  dplyr::select(ID_ENTITY, site_name, elevation, latitude, longitude) %>% 
  dplyr::left_join(pollen_basin, by = "ID_ENTITY") %>% #add basin info
  dplyr::filter(!is.na(basin_size)) %>% #just records with basin info
  dplyr::group_by(site_name) %>% 
  dplyr::summarise(ID_ENTITY = min(ID_ENTITY), elevation = mean(elevation), basin_size = mean(basin_size), latitude = mean(latitude), longitude = mean(longitude)) %>% #same for each site
  dplyr::ungroup() 

rio::export(pollen_single, "data/intermediate_output/vegetation/pollen_single.csv")
  
pollen_single_per <- pollen %>% #single record per site and then by coords, percentages recalculation
  dplyr::select(site_name, ID_ENTITY, taxa_names[1]:last(taxa_names)) %>%
  dplyr::left_join(pollen_basin, by = "ID_ENTITY") %>% #add basin info
  dplyr::filter(!is.na(basin_size)) %>% #just records with basin info
  dplyr::select(-basin_size) %>%  
  dplyr::group_by(site_name) %>% 
  dplyr::summarise(ID_ENTITY = min(ID_ENTITY), dplyr::across(taxa_names[1]:last(taxa_names), sum)) %>% #average across records at same site
  dplyr::ungroup() %>% 
  dplyr::left_join(dplyr::select(pollen, ID_ENTITY, latitude, longitude), by = "ID_ENTITY") %>% 
  dplyr::group_by(latitude, longitude) %>% 
  dplyr::summarise(ID_ENTITY = min(ID_ENTITY), dplyr::across(taxa_names[1]:last(taxa_names), sum)) %>% #average across all locations with same location info
  dplyr::ungroup() %>% 
  dplyr::select(ID_ENTITY, taxa_names[1]:last(taxa_names)) %>%     
  smpds::normalise_taxa(cols = 1)

rio::export(pollen_single_per, "data/intermediate_output/vegetation/pollen_single_per.csv")


#Meta info about the pollen records
pollen_records <- pollen_single %>% 
  dplyr::left_join(dplyr::select(pollen, ID_ENTITY, ID_SAMPLE, entity_name, site_type, entity_type, age_BP, publication, doi), by = "ID_ENTITY") %>% #add extra meta info
  dplyr::mutate(age_BP = dplyr::if_else(age_BP == -200, "modern", as.character(age_BP))) %>% #Change to numeric
  dplyr::select(ID_SAMPLE, site_name, ID_ENTITY, entity_name, latitude, longitude, elevation, basin_size, site_type, entity_type, basin_size, age_BP, publication, doi)

rio::export(pollen_records, "data/intermediate_output/vegetation/pollen_records_meta.xlsx")



# For each basin size, calculate the correct buffer

# This is based on Prentice's formula 1988, and uses average tree pollen fall speeds published in Githumbi. The range for these fall speeds are 
# a)0.01 for Castanea sativa to b)0.12 for Abies alba, with a mean of c)0.034 and median of d)0.030. 
# Based on 75*FSP/3, with 3 representing the wind speed, these fall speeds generate bi values of a)0.25; b)3; c)0.850; and d)0.75
# Prentice's forumula calculated the 70% characteristic source as: (radius of basin^(1/8) - (ln 0.3/bi))^(1/1/8).
# Radius is calculated as (basin area / pi)^(1/2)

#Function to calculate appropriate source area
f_buffer <- function(data,b,area){
  data_out <- data %>% 
  dplyr::mutate("i{{b}}" := (((1000000*area/pi)^(1/2))^(0.125) - (log(0.3)/b))^(1/0.125)) #need to convert area to m; epsilon 0.3 (i.e 70% radius)
  assign(paste0("pollen_buffer_", b), data_out, envir = parent.frame()) #store output
}

pollen_buffer_calc <- pollen_single  %>% #select pollen sites with some info on basin size. All pollen data. 
  dplyr::select(ID_ENTITY, latitude, longitude, basin_size)  

f_buffer(pollen_buffer_calc,0.75,pollen_buffer_calc$basin_size) #this is the median based on Githumbi (and RPP values 0.0305)

pollen_buffer_variable <- pollen_buffer_0.75 %>% #Merge information for each FSP
  dplyr::mutate(basin_radius = (basin_size/pi)^(1/2), .before = i0.75) #km, add basin radius



# ---------------------------------------------------------


# 3. Shannon index
# ---------------------------------------------------------
# ---------------------------------------------------------

#Shannon index calc (only those in TPS)
veg_shannon <- pollen_single_per %>% 
  tidyr::pivot_longer(!c(ID_ENTITY), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% #re-shape 
  dplyr::left_join(taxa_cat_single, by = "clean_taxon_name") %>%     #all species tps needed
  dplyr::filter(terrestrial_pollen_sum == "yes")  %>%   #select only those records within TPS
  dplyr::mutate(europe = dplyr::if_else(is.na(europe), "none", europe)) %>% 
  dplyr::filter(europe != "not native to Europe") %>%    #take out species that are not native
  dplyr::select(ID_ENTITY, clean_taxon_name, percentage_cover) %>%   
  tidyr::pivot_wider(names_from = clean_taxon_name, values_from = percentage_cover) 

veg_shannon1 <- pollen_single_per %>%
  dplyr::select(ID_ENTITY, taxa_names[1]:last(taxa_names)) %>% 
  tidyr::pivot_longer(!c(ID_ENTITY), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% #re-shape 
  dplyr::left_join(taxa_cat_single, by = "clean_taxon_name") %>%     #all species tps needed
  dplyr::filter(terrestrial_pollen_sum == "yes")  %>%   #select only those records within TPS
  dplyr::mutate(europe = dplyr::if_else(is.na(europe), "none", europe)) %>% 
  dplyr::filter(europe != "not native to Europe") %>%    #take out species that are not native
  dplyr::select(ID_ENTITY, clean_taxon_name, percentage_cover) %>%   
  tidyr::pivot_wider(names_from = clean_taxon_name, values_from = percentage_cover) %>% 
  dplyr::select(-ID_ENTITY) %>% 
  as.matrix()

veg_shannon2 <- vegan::diversity(veg_shannon1, index = "shannon") %>% 
  dplyr::as_tibble()

veg_shannon3 <- veg_shannon %>% 
  dplyr::select(ID_ENTITY) %>% 
  dplyr::mutate(shannon = veg_shannon2$value)

rio::export(veg_shannon3, "data/intermediate_output/vegetation/veg_shannon3.csv")
rm(veg_shannon3,veg_shannon2,veg_shannon1,veg_shannon)

tree_shannon <- pollen_single_per %>% 
  tidyr::pivot_longer(!c(ID_ENTITY), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% #re-shape 
  dplyr::left_join(taxa_cat_single, by = "clean_taxon_name") %>%     #all species tps needed
  dplyr::filter(terrestrial_pollen_sum == "yes")  %>%   #select only those records within TPS
  dplyr::mutate(europe = dplyr::if_else(is.na(europe), "none", europe)) %>% 
  dplyr::filter(europe != "not native to Europe") %>%     #take out species that are not native
  dplyr::filter(ap_sp_hp == "AP") %>% #limit to trees  
  dplyr::select(ID_ENTITY, clean_taxon_name, percentage_cover) %>%   
  tidyr::pivot_wider(names_from = clean_taxon_name, values_from = percentage_cover) 

tree_shannon1 <- pollen_single_per %>%
  dplyr::select(ID_ENTITY, taxa_names[1]:last(taxa_names)) %>% 
  tidyr::pivot_longer(!c(ID_ENTITY), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% #re-shape 
  dplyr::left_join(taxa_cat_single, by = "clean_taxon_name") %>%     #all species tps needed
  dplyr::filter(terrestrial_pollen_sum == "yes")  %>%   #select only those records within TPS
  dplyr::mutate(europe = dplyr::if_else(is.na(europe), "none", europe)) %>% 
  dplyr::filter(europe != "not native to Europe") %>%    #take out species that are not native
  dplyr::filter(ap_sp_hp == "AP") %>% #limit to trees 
  dplyr::select(ID_ENTITY, clean_taxon_name, percentage_cover) %>%   
  tidyr::pivot_wider(names_from = clean_taxon_name, values_from = percentage_cover) %>% 
  dplyr::select(-ID_ENTITY) %>% 
  as.matrix()

tree_shannon2 <- vegan::diversity(tree_shannon1, index = "shannon") %>% 
  dplyr::as_tibble()

tree_shannon3 <- tree_shannon %>% 
  dplyr::select(ID_ENTITY) %>% 
  dplyr::mutate(tree_shannon = tree_shannon2$value)

rio::export(tree_shannon3, "data/intermediate_output/vegetation/tree_shannon3.csv")
rm(tree_shannon3,tree_shannon2,tree_shannon1,tree_shannon)



# ---------------------------------------------------------
# ---------------------------------------------------------





# 4. Extract tree cover values based on pollen site locations
# ---------------------------------------------------------

#Pollen from SMPDS analysis
pollen2 <- pollen_single %>% 
  dplyr::left_join(dplyr::select(pollen_buffer_variable, ID_ENTITY, i0.75), by = c("ID_ENTITY")) %>% #Incorporate basin size info
  dplyr::distinct() #To ensure no clear duplicates
  
pollen_records2 <- pollen2 %>% #Get site based info about long and lat. Needed because of multi-record sites
  dplyr::select(ID_ENTITY, latitude, longitude, i0.75) %>%
  dplyr::arrange(ID_ENTITY)  #Standardise ordering

rio::export(pollen_records2, "data/intermediate_output/vegetation/pollen_records2.csv")

pollen_records_xy <- pollen_records2 %>% #To use in extraction
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258)

pollen_records_xy2 <- pollen_records2 %>% #To use in extraction
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>%
  sf::st_transform(crs = 25830)

pollen_records_xy3 <- pollen_records2 %>% #This is for the Copernicus tree cover CRS
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>%
  sf::st_transform(crs = 3035)

pollen_records_xy_df2 <- pollen_records2 %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>%
  sf::st_transform(crs = 25830) %>%
  sf::st_coordinates() %>%
  dplyr::as_tibble() %>%  
  dplyr::mutate(ID_ENTITY = pollen_records2$ID_ENTITY)

pollen_records_xy_df3 <- pollen_records2 %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4258) %>%
  sf::st_transform(crs = 3035) %>%
  sf::st_coordinates() %>%
  dplyr::as_tibble() %>% 
  dplyr::rename(X3035 = X, Y3035 = Y)  %>%  
  dplyr::mutate(ID_ENTITY = pollen_records2$ID_ENTITY)

pollen_records_xyboth <- pollen_records_xy2 %>% #All crs of pollen data
  sf::st_coordinates() %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(ID_ENTITY = pollen_records2$ID_ENTITY) %>%
  dplyr::rename(X25830 = X) %>% 
  dplyr::rename(Y25830 = Y) %>% 
  dplyr::mutate(longitude = pollen_records2$longitude) %>%
  dplyr::mutate(latitude = pollen_records2$latitude) %>%
  dplyr::mutate(X3035 = pollen_records_xy_df3$X3035) %>% 
  dplyr::mutate(Y3035 = pollen_records_xy_df3$Y3035) 

euro_map_3035 <- sf::st_read("data/intermediate_output/vegetation/euro_map/3035/euro_map_3035.shp") #load map (1.Map_setup)
ggplot2::ggplot(data = euro_map_3035)+
  geom_sf()+
  geom_sf(fill = "seashell") +
  theme(panel.background = element_rect(fill = "aliceblue"))+
  geom_point(pollen_records_xy_df3, mapping = aes(x = X3035, y = Y3035), color = "red")

pollen_records_xy_buffer4 <- pollen_records_xy3 %>% #variable buffer
  sf::st_buffer(., pollen_records_xy3$i0.75) %>% #based on basin_size and a median FSP from Githumbi
  dplyr::arrange(ID_ENTITY) %>% 
  dplyr::select(-c(i0.75)) %>%
  sf::st_as_sf(crs = 3035) #Ensure in sf

saveRDS(pollen_records_xy_buffer4, "data/intermediate_output/vegetation/pollen_records_xy_buffer4.rds")
pollen_records_xy_buffer4 <- readRDS("data/intermediate_output/vegetation/pollen_records_xy_buffer4.rds")

pollen_records_xy_buffer4_area <- sf::st_area(pollen_records_xy_buffer4) %>% 
  dplyr::as_tibble()


#Extract data from land cover data
# Define function: Value - mean
f_extract_veg_mean <- function(h){
  extracted <- exactextractr::exact_extract(h, pollen_records_xy_buffer4, "mean", force_df = TRUE) %>%  #with buffer
    dplyr::mutate(ID_ENTITY = pollen_records_xyboth$ID_ENTITY) %>% 
    dplyr::mutate(X25830 = pollen_records_xyboth$X25830) %>% 
    dplyr::mutate(Y25830 = pollen_records_xyboth$Y25830)  %>%
    dplyr::mutate(X3035 = pollen_records_xyboth$X3035) %>% 
    dplyr::mutate(Y3035 = pollen_records_xyboth$Y3035)  %>%
    dplyr::mutate(latitude = pollen_records2$latitude) %>% 
    dplyr::mutate(longitude = pollen_records2$longitude) 
  extracted_all <- extracted %>% 
    dplyr::filter(!is.na(mean))
}

# Define function: Cells and na counts
f_extract_na_prop <- function(h){
  extracted_tree_cells_num <- exactextractr::exact_extract(h, pollen_records_xy_buffer4, function(values, coverage_fraction) length(values))
  extracted_tree_cells_numna <- exactextractr::exact_extract(h, pollen_records_xy_buffer4, function(values, coverage_fraction) sum(is.na(values)))
  extracted_tree_cells_propna <- extracted_tree_cells_numna/extracted_tree_cells_num %>% 
    dplyr::as_tibble() %>% 
    dplyr::rename(prop_na = value) 
  extracted_cells_df_na <- extracted_tree_cells_propna %>% 
    dplyr::mutate(ID_ENTITY = pollen_records_xyboth$ID_ENTITY)
}

#load in rasters for extraction
cop_tree_cover_masked_3035 <-  terra::rast("data/intermediate_output/vegetation/copernicus/cop_masked_tree_cover_3035_100m.tif") #masked tree cover

na_prop_extracted <- f_extract_na_prop(cop_tree_cover_masked_3035)
rio::export(na_prop_extracted, "data/intermediate_output/vegetation/na_prop_extracted.csv")
rm(na_prop_extracted)

mean_extracted <- f_extract_veg_mean(cop_tree_cover_masked_3035)
rio::export(mean_extracted, "data/intermediate_output/vegetation/mean_extracted.csv")
rm(mean_extracted)

#All pollen
mean_extracted <- rio::import("data/intermediate_output/vegetation/mean_extracted.csv")
na_prop_extracted <- rio::import("data/intermediate_output/vegetation/na_prop_extracted.csv")

pollen_extracted_info <- pollen_single_per %>% #start with pollen info with basin values  
  dplyr::left_join(pollen_single, by = "ID_ENTITY") %>% 
  dplyr::left_join(dplyr::select(mean_extracted, -c(latitude, longitude)), by = "ID_ENTITY") %>%
  dplyr::left_join(na_prop_extracted, by = "ID_ENTITY") 

rio::export(pollen_extracted_info, "data/intermediate_output/vegetation/pollen_extracted_info.csv")

