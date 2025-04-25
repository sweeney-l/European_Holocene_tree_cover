#Sweeney, Harrison and Vander Linden 2024. 
#European forest cover during the Holocene reconstructed from pollen records
# ---------------------------------------------------------

#There are six seperate scripts associated with this research:
# 1. Map_setup.R
# 2. Modern_pollen.R (this script)
# 3. Model_import_data.R
# 4. Tree_model.R
# 5. Tree_recon.R
# 6. Figures.R

# ---------------------------------------------------------

############ Modern pollen ############ 

# This script imports the modern pollen data and reshapes the data for later
# use


# ---------------------------------------------------------
# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Set map extents
# 3. Import  and treat SMPDS data

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

# ---------------------------------------------------------




# 3. Import  and treat SMPDS data
# ---------------------------------------------------------
# ---------------------------------------------------------
data("entity", package = "smpds")
data("pollen_count", package = "smpds")
data("taxon_name", package = "smpds")

#Rename data
SMPDS_entity <- entity
SMPDS_pollen_count <- pollen_count %>% 
  dplyr::filter(amalgamation_level == 2) #limit to higher level species grouping
SMPDS_taxon_name <- taxon_name

taxa_cat <- rio::import("data/input/taxa_cat.csv") %>% #import information regarding the taxa
  dplyr::arrange(taxon_name) 
taxa_cat_single <- taxa_cat %>% #for single clean_taxon_name
  dplyr::select(-taxon_name) %>% 
  dplyr::distinct()

SMPDS_combined <- SMPDS_pollen_count %>% #merge datatables
  dplyr::left_join(SMPDS_taxon_name, by = "ID_TAXON") %>%
  dplyr::left_join(SMPDS_entity, by = "ID_SAMPLE")
 
SMPDS_meta <- SMPDS_combined %>% #record metadata
  dplyr::select(ID_SAMPLE, site_name, ID_ENTITY, entity_name, latitude, longitude, elevation, entity_type, site_type, age_BP, basin_size, publication, doi) %>%
  dplyr::distinct()

#Identify modern records
SMPDS_modern <- SMPDS_combined %>% #select those listed as modern and assign youngest date
  dplyr::filter(age_BP %in%  c("modern", "assumed modern", "Modern")) %>% 
  dplyr::mutate(age_BP = -200) #set this as the youngest to use

SMPDS_50 <- SMPDS_combined %>% #Select those with an age since 1950
  dplyr::filter(!age_BP %in% c("modern","assumed modern", "Modern")) %>% #exclude SMPDS_modern
  dplyr::mutate(age_BP = as.numeric(age_BP)) %>% 
  dplyr::filter(age_BP <= 0) #1950

SMPDS_modern50 <- rbind(SMPDS_modern, SMPDS_50) %>%  #join those with dates in last 50 years and those called "modern"
  dplyr::arrange(ID_SAMPLE)

#Add information about the taxa
SMPDS_cleaner <- SMPDS_modern50 %>% #Add info about the taxa
  dplyr::left_join(taxa_cat, by = "taxon_name")

SMPDS_mod_euro <- SMPDS_cleaner %>% #limit to area of interest
  dplyr::filter(dplyr::between(longitude,euro_extent_4258_xmin,euro_extent_4258_xmax) & #set for the long lat
                  dplyr::between(latitude,euro_extent_4258_ymin, euro_extent_4258_ymax)) %>% 
  dplyr::arrange(clean_taxon_name)

missing_taxon_cleaner <- SMPDS_mod_euro %>% #check whether there were any missing taxa from catagorisation
  dplyr::filter(is.na(clean_taxon_name)) %>% 
  dplyr::select(taxon_name) %>% 
  dplyr::distinct()

#Filter taxa to be highest amalgamation level and part of TPS
SMPDS_mod_euro_meta <- SMPDS_mod_euro %>% 
  dplyr::select(ID_SAMPLE, site_name, ID_ENTITY, entity_name, latitude, longitude, elevation, entity_type, site_type, age_BP, basin_size, publication, doi) %>% 
  dplyr::distinct()

SMPDS_am2 <- SMPDS_mod_euro %>% 
  dplyr::filter(terrestrial_pollen_sum == "yes") %>% #only select those considered as part of the terrestrial pollen sum
  dplyr::group_by(ID_ENTITY, clean_taxon_name) %>% 
  dplyr::summarise(count = mean(count), ID_SAMPLE = min(ID_SAMPLE)) %>% #To ensure duplicate record counts are averaged
  dplyr::ungroup() %>% 
  dplyr::group_by(ID_SAMPLE, clean_taxon_name) %>%    
  dplyr::summarise(count = sum(count)) %>%   #To make sure that any duplicate higher level taxa are summed here
  dplyr::ungroup() %>% 
  dplyr::select(ID_SAMPLE, count, clean_taxon_name) %>%    
  dplyr::left_join(dplyr::select(taxa_cat_single, clean_taxon_name, europe), by = "clean_taxon_name") %>%   
  dplyr::filter(europe != "not native to Europe") %>%  #remove non native species (but include ferns etc.)
  dplyr::select(-europe)  %>% 
  tidyr::pivot_wider(names_from = clean_taxon_name, values_from = count, values_fill = 0) %>% #set to zero any NA values
  dplyr::select(ID_SAMPLE, sort(colnames(.))) #make sure in order

#Hills check
Hills_am2 <- SMPDS_am2 %>% 
  dplyr::select(Abies:dplyr::last_col()) %>% 
  analogue::n2(., "sites")


#Finalise data
updated_basin_info <- rio::import("data/input/SMPDS/SMPDSv2_updated_meta_info.csv") %>% #Load in updates to SMPDS
  dplyr::rename(basin_size1 = basin_size, elevation1 = elevation, site_name1 = site_name, site_type1 = site_type, entity_type1 = entity_type)
updated_basin_info[updated_basin_info == ""] <- NA #convert all empty cells to NA
SMPDS_am2_per <- SMPDS_am2 %>% #Calculate percentage covers
  dplyr::mutate(hills_N2 = Hills_am2) %>%   #add Hills information
  dplyr::filter(hills_N2 >= 2) %>% #Filter based on Hills as per Wei 2021
  dplyr::select(!hills_N2) %>% #remove Hills info
  smpds::normalise_taxa(cols = 1) %>% #calculate percentages, ignoring first column 
  dplyr::left_join(SMPDS_mod_euro_meta, by = "ID_SAMPLE") %>% #add meta data
  dplyr::select(ID_SAMPLE,
                site_name,
                ID_ENTITY, 
                entity_name, 
                latitude, 
                longitude, 
                elevation, 
                site_type,
                entity_type,
                basin_size,
                age_BP,
                publication,
                doi,
                dplyr::everything()) %>% 
  dplyr::mutate(site_type = dplyr::if_else(is.na(site_type), "not known", site_type)) %>% 
  dplyr::filter(!site_type %in% c("archaeological site", 
                                  "cave",
                                  "coastal, estuarine",
                                  "coastal, lagoon",
                                  "glacial",
                                  "fluvial",
                                  "marine")) %>%     #take out definitively marine and other sites
  dplyr::left_join(dplyr::select(updated_basin_info, ID_SAMPLE, basin_size1, site_type1, elevation1, entity_type1, site_name1, New_lat, New_lon), by = "ID_SAMPLE") %>% #Add basin info
  dplyr::mutate(basin_size = dplyr::if_else(!is.na(basin_size1), basin_size1, basin_size)) %>% 
  dplyr::mutate(site_type = dplyr::if_else(!is.na(site_type1), site_type1, site_type)) %>%
  dplyr::mutate(elevation = dplyr::if_else(!is.na(elevation1), elevation1, elevation)) %>%
  dplyr::mutate(entity_type = dplyr::if_else(!is.na(entity_type1), entity_type1, entity_type)) %>%
  dplyr::mutate(site_name = dplyr::if_else(!is.na(site_name1), site_name1, site_name)) %>%
  dplyr::mutate(latitude = dplyr::if_else(!is.na(New_lat), New_lat, latitude)) %>% #Update lats and lons
  dplyr::mutate(longitude = dplyr::if_else(!is.na(New_lon), New_lon, longitude)) %>%
  dplyr::select(-c(basin_size1, site_type1, elevation1, entity_type1, site_name1, New_lat, New_lon)) %>%
  dplyr::mutate(basin_size = dplyr::if_else(entity_type %in% c("moss polster or moss", "pollen trap"), "0.00005", basin_size)) #add 50m2 value for all moss polster and pollen trap

rio::export(SMPDS_am2_per, "data/intermediate_output/vegetation/SMPDS/SMPDS_am2_per.csv") #Final version for use
