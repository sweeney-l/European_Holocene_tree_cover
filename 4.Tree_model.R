#Sweeney, Harrison and Vander Linden 2024. 
#European forest cover during the Holocene reconstructed from pollen records
# ---------------------------------------------------------

#There are six seperate scripts associated with this research:
# 1. Map_setup.R 
# 2. Modern_pollen.R 
# 3. Model_import_data.R
# 4. Tree_model.R (this script)
# 5. Tree_recon.R
# 6. Figures.R

# ---------------------------------------------------------

############ Tree model ############ 

# This script sets develops the modern tree model


# ---------------------------------------------------------
# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Import and set up data
# 3. Run model and analysis 
# 4. LOOCV

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




# 2. Import and set up data
# ---------------------------------------------------------
# ---------------------------------------------------------
pollen <- rio::import("data/intermediate_output/vegetation/SMPDS/SMPDS_am2_per.csv") #all of the pollen data
pollen_records <- rio::import("data/intermediate_output/vegetation/pollen_records_meta.xlsx") #meta info about records
pollen_single_per <- rio::import("data/intermediate_output/vegetation/pollen_single_per.csv")
pollen_buffer_variable <- rio::import("data/intermediate_output/vegetation/pollen_buffer_variable.xlsx")
veg_shannon3 <- rio::import("data/intermediate_output/vegetation/veg_shannon3.csv")
tree_shannon3 <- rio::import("data/intermediate_output/vegetation/tree_shannon3.csv")
pollen_extracted_info <- rio::import("data/intermediate_output/vegetation/pollen_extracted_info.csv")

taxa_names <- pollen %>% #taxa names for pollen records selected
  dplyr::select(Abies: ncol(.)) %>% 
  colnames(.)

taxa_cat <- rio::import("data/input/taxa_cat.csv") %>% #import information regarding the taxa
  dplyr::arrange(taxon_name) 
taxa_cat_single <- taxa_cat %>% #for single clean_taxon_name
  dplyr::select(-taxon_name) %>% 
  dplyr::distinct()

pollen_3035 <- pollen_extracted_info %>% #location of all data, by type
  dplyr::select("ID_ENTITY", "X3035", "Y3035") %>% 
  dplyr::left_join(dplyr::select(pollen_records, ID_ENTITY, site_type, entity_type), by = "ID_ENTITY")

euro_map_3035 <- sf::st_read("data/intermediate_output/vegetation/euro_map/3035/euro_map_3035.shp") #load map (1.Map_setup)
ggplot2::ggplot(data = euro_map_3035)+ #map of all pollen info, by type of site
  geom_sf()+
  geom_sf(fill = "seashell") +
  theme(panel.background = element_rect(fill = "aliceblue"))+
  geom_point(pollen_3035, mapping = aes(x = X3035, y = Y3035, color = site_type))

pollen_records_xy_buffer4 <- readRDS("data/intermediate_output/vegetation/pollen_records_xy_buffer4.rds")
ggplot2::ggplot(data = euro_map_3035)+ #map of all pollen info, by type of site
  geom_sf()+
  geom_sf(fill = "seashell") +
  theme(panel.background = element_rect(fill = "aliceblue"))+
  geom_sf(pollen_records_xy_buffer4, mapping = aes(color = "red"))

#Build tree model
sp_cover <- pollen_extracted_info %>% #SP cover
  dplyr::select(ID_ENTITY, mean, prop_na, taxa_names[1]:last(taxa_names)) %>%
  dplyr::rename(tree = 2, prop_na = 3) %>%
  tidyr::drop_na() %>% 
  tidyr::pivot_longer(!c(ID_ENTITY, tree, prop_na), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% 
  dplyr::left_join(taxa_cat_single, by = "clean_taxon_name") %>% 
  dplyr::filter(terrestrial_pollen_sum == "yes") %>% #select only those taxa within TPS
  dplyr::filter(ap_sp_hp == "SP") %>%      #just SP for % calculation
  dplyr::group_by(ID_ENTITY) %>% 
  dplyr::summarise(sp_cover = sum(percentage_cover)) %>% #calculate % cover
  dplyr::ungroup() %>% 
  dplyr::mutate(sp_cover = sp_cover/100) 

rio::export(sp_cover, "data/intermediate_output/vegetation/sp_cover.csv")

ap_tree_cover <- pollen_extracted_info %>% 
  dplyr::select(ID_ENTITY, mean, prop_na, taxa_names[1]:last(taxa_names)) %>%
  dplyr::rename(tree = 2, prop_na = 3) %>% 
  tidyr::drop_na() %>%     
  tidyr::pivot_longer(!c(ID_ENTITY, tree, prop_na), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% 
  dplyr::left_join(taxa_cat_single, by = "clean_taxon_name") %>%
  dplyr::filter(terrestrial_pollen_sum == "yes") %>% #select only those taxa within TPS
  dplyr::mutate(europe = dplyr::if_else(is.na(europe), "none", europe)) %>%
  dplyr::filter(europe != "not native to Europe") %>% #take out non-natives
  dplyr::filter(ap_sp_hp == "AP") %>%    #just AP for total AP% calculation
  dplyr::mutate(needle = dplyr::if_else(is.na(ap_needle_broad), "unknown", ap_needle_broad)) %>%   #need to remove NA
  dplyr::left_join(dplyr::select(pollen_records, ID_ENTITY, site_type, entity_type), by = "ID_ENTITY") %>% #add site meta info
  dplyr::filter(site_type %in% c("lake", "terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, marsh")) %>% #Limit analysis to Lakes and bog(ish)
  dplyr::filter(!entity_type %in% c("litter", "soil sample","pollen trap", "moss polster or moss")) %>%
  dplyr::group_by(ID_ENTITY) %>% 
  dplyr::summarise(tree = max(tree), prop_na = max(prop_na), ap_cover = sum(percentage_cover), needle = sum(percentage_cover[needle == "needle"])) %>% #calculate % cover
  dplyr::ungroup() %>% 
  dplyr::mutate(difference = ap_cover - tree) %>%
  dplyr::mutate(tree = tree/100, ap_cover = ap_cover/100, needle = needle/100) %>% 
  dplyr::mutate(needle_share = needle/ap_cover) %>% #share of total AP that is needleleaf
  dplyr::left_join(dplyr::select(pollen_extracted_info, ID_ENTITY, X3035, Y3035, elevation, latitude,longitude, basin_size), by = "ID_ENTITY") %>% #add coordinate info
  dplyr::left_join(dplyr::select(pollen_records, ID_ENTITY, site_type, entity_type), by = "ID_ENTITY") %>% #add site meta info
  dplyr::filter(!(basin_size > 0.5 & site_type %in% c("terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, marsh"))) %>% #filter bogs larger than 0.5 (400m radius as per Githumbi)
  dplyr::filter(elevation <1000) %>% #add elevation filter
  dplyr::filter(elevation <500) %>% #add elevation filter
  dplyr::mutate(site_model = dplyr::if_else(site_type == "lake", 1, 0)) %>% 
  dplyr::left_join(tree_shannon3, by = "ID_ENTITY") %>%  #add info about species diversity in pollen record (TPS only)
  dplyr::left_join(sp_cover, by = "ID_ENTITY") %>% 
  dplyr::filter(ap_cover > 0) %>%            #because needleleaf share will be NaN for single record
  dplyr::filter(prop_na < 0.5)      #limit associated with no of NA in extracted tree cover

rio::export(ap_tree_cover, "data/intermediate_output/vegetation/ap_tree_cover.csv")

variable_buff_utilsed <- ap_tree_cover %>% 
  dplyr::left_join(dplyr::select(pollen_buffer_variable, ID_ENTITY, i0.75), by = "ID_ENTITY")
max(variable_buff_utilsed$i0.75)
min(variable_buff_utilsed$i0.75)
median(variable_buff_utilsed$i0.75)

#Lake only
ap_tree_cover_lake <- pollen_extracted_info %>% 
  dplyr::select(ID_ENTITY, mean, prop_na, taxa_names[1]:last(taxa_names)) %>%
  dplyr::rename(tree = 2, prop_na = 3) %>% 
  tidyr::drop_na() %>%     
  tidyr::pivot_longer(!c(ID_ENTITY, tree, prop_na), names_to = "clean_taxon_name", values_to = "percentage_cover") %>% 
  dplyr::left_join(taxa_cat_single, by = "clean_taxon_name") %>%
  dplyr::filter(terrestrial_pollen_sum == "yes") %>% #select only those taxa within TPS
  dplyr::mutate(europe = dplyr::if_else(is.na(europe), "none", europe)) %>%
  dplyr::filter(europe != "not native to Europe") %>% #take out non-natives
  dplyr::filter(ap_sp_hp == "AP") %>%    #just AP for total AP% calculation
  dplyr::mutate(needle = dplyr::if_else(is.na(ap_needle_broad), "unknown", ap_needle_broad)) %>%   #need to remove NA
  dplyr::left_join(dplyr::select(pollen_records, ID_ENTITY, site_type, entity_type), by = "ID_ENTITY") %>% #add site meta info
  dplyr::filter(site_type %in% c("lake")) %>% #Limit analysis to Lakes
  dplyr::filter(!entity_type %in% c("litter", "soil sample","pollen trap", "moss polster or moss")) %>%
  dplyr::group_by(ID_ENTITY) %>% 
  dplyr::summarise(tree = max(tree), prop_na = max(prop_na), ap_cover = sum(percentage_cover), needle = sum(percentage_cover[needle == "needle"])) %>% #calculate % cover
  dplyr::ungroup() %>% 
  dplyr::mutate(difference = ap_cover - tree) %>%
  dplyr::mutate(tree = tree/100, ap_cover = ap_cover/100, needle = needle/100) %>% 
  dplyr::mutate(needle_share = needle/ap_cover) %>% #share of total AP that is needleleaf
  dplyr::left_join(dplyr::select(pollen_extracted_info, ID_ENTITY, X3035, Y3035, elevation, latitude,longitude, basin_size), by = "ID_ENTITY") %>% #add coordinate info
  dplyr::left_join(dplyr::select(pollen_records, ID_ENTITY, site_type, entity_type), by = "ID_ENTITY") %>% #add site meta info
  dplyr::filter(!(basin_size > 0.5 & site_type %in% c("terrestrial, blanket bog", "terrestrial, bog/fen/swamp", "terrestrial, bog/lake", "terrestrial, marsh"))) %>% #filter bogs larger than 0.5 (400m radius as per Githumbi)
  dplyr::filter(elevation <1000) %>% #add elevation filter
  dplyr::mutate(site_model = dplyr::if_else(site_type == "lake", 1, 0)) %>% 
  dplyr::left_join(tree_shannon3, by = "ID_ENTITY") %>%  #add info about species diversity in pollen record (TPS only)
  dplyr::left_join(sp_cover, by = "ID_ENTITY") %>% 
  dplyr::filter(ap_cover > 0) %>%         #because needleleaf share will be NaN for single record
  dplyr::filter(prop_na < 0.5)      #limit associated with no of NA in extracted tree cover

# ---------------------------------------------------------


# 3. Run model and analysis
# ---------------------------------------------------------
# ---------------------------------------------------------
summary(ap_tree_cover$tree)
hist(ap_tree_cover$tree)
tree_ap_correlation <- cor(ap_tree_cover$tree, ap_tree_cover$ap_cover)

tree_beta <- betareg::betareg(tree~ap_cover*elevation+needle_share+I(needle_share^2)+tree_shannon*elevation+I(tree_shannon^2)*elevation+site_model*elevation+sp_cover*elevation|needle_share+tree_shannon+site_model, data = dplyr::filter(ap_tree_cover, tree > 0)) #final normal
tree_beta_stat <- betareg::betareg(tree~ap_cover*elevation+needle_share+I(needle_share^2)+tree_shannon*elevation+I(tree_shannon^2)*elevation+site_model*elevation+sp_cover*elevation, data = dplyr::filter(ap_tree_cover, tree > 0)) #without variable precision
tree_beta_poly <- betareg::betareg(tree~ap_cover*elevation+poly(needle_share,2)+poly(tree_shannon,2)*elevation+site_model*elevation+sp_cover*elevation|needle_share+tree_shannon+site_model, data = dplyr::filter(ap_tree_cover, tree > 0)) #final normal poly
tree_beta_no_interaction <- betareg::betareg(tree~ap_cover+poly(needle_share,2)+poly(tree_shannon,2)+site_model+sp_cover+elevation|needle_share+tree_shannon+site_model, data = dplyr::filter(ap_tree_cover, tree > 0)) #model with no interaction effects
tree_beta_lake <- betareg::betareg(tree~ap_cover*elevation+needle_share+I(needle_share^2)+tree_shannon*elevation+I(tree_shannon^2)*elevation+sp_cover*elevation|needle_share+tree_shannon,data = dplyr::filter(ap_tree_cover_lake, tree > 0)) #final normal with just lake sites

tree_beta_no_ap_cover <- betareg::betareg(tree~needle_share+I(needle_share^2)+tree_shannon*elevation+I(tree_shannon^2)*elevation+site_model*elevation+sp_cover*elevation|needle_share+tree_shannon+site_model, data = dplyr::filter(ap_tree_cover, tree > 0)) #final normal no ap
tree_beta_no_sp_cover <- betareg::betareg(tree~ap_cover*elevation+needle_share+I(needle_share^2)+tree_shannon*elevation+I(tree_shannon^2)*elevation+site_model*elevation|needle_share+tree_shannon+site_model, data = dplyr::filter(ap_tree_cover, tree > 0)) #final normal no sp
tree_beta_no_elevation <- betareg::betareg(tree~ap_cover+needle_share+I(needle_share^2)+tree_shannon+I(tree_shannon^2)+site_model+sp_cover|needle_share+tree_shannon+site_model, data = dplyr::filter(ap_tree_cover, tree > 0)) #final normal no elevation
tree_beta_no_needle_share <- betareg::betareg(tree~ap_cover*elevation+tree_shannon*elevation+I(tree_shannon^2)*elevation+site_model*elevation+sp_cover*elevation|tree_shannon+site_model, data = dplyr::filter(ap_tree_cover, tree > 0)) #final normal no needle share
tree_beta_no_tree_shannon <- betareg::betareg(tree~ap_cover*elevation+needle_share+I(needle_share^2)+site_model*elevation+sp_cover*elevation|needle_share+site_model, data = dplyr::filter(ap_tree_cover, tree > 0)) #final normal no tree shannon
tree_beta_no_site_model <- betareg::betareg(tree~ap_cover*elevation+needle_share+I(needle_share^2)+tree_shannon*elevation+I(tree_shannon^2)*elevation+sp_cover*elevation|needle_share+tree_shannon, data = dplyr::filter(ap_tree_cover, tree > 0)) #final normal no site model
tree_beta_only_ap_sp <- betareg::betareg(tree~ap_cover+sp_cover, data = dplyr::filter(ap_tree_cover, tree > 0)) #final normal


saveRDS(tree_beta, file = "data/intermediate_output/vegetation/tree_beta.rda") #store model
saveRDS(tree_beta_stat, file = "data/intermediate_output/vegetation/tree_beta_stat.rda") #store model without precision variables
saveRDS(tree_beta_poly, file = "data/intermediate_output/vegetation/tree_beta_poly.rda") #store model with polnomial fit
saveRDS(tree_beta_no_interaction, file = "data/intermediate_output/vegetation/tree_beta_no_interaction.rda") #store model with no interaction efefcts
saveRDS(tree_beta_lake, file = "data/intermediate_output/vegetation/tree_beta_lake.rda") #store lake model
saveRDS(tree_beta_no_ap_cover, file = "data/intermediate_output/vegetation/tree_beta_no_ap_cover.rda") #store model without ap
saveRDS(tree_beta_no_sp_cover, file = "data/intermediate_output/vegetation/tree_beta_no_sp_cover.rda") #store model without sp
saveRDS(tree_beta_no_elevation, file = "data/intermediate_output/vegetation/tree_beta_no_elevation.rda") #store model without elevation
saveRDS(tree_beta_no_needle_share, file = "data/intermediate_output/vegetation/tree_beta_no_needle_share.rda") #store model without needle share
saveRDS(tree_beta_no_tree_shannon, file = "data/intermediate_output/vegetation/tree_beta_no_tree_shannon.rda") #store model without tree shannon
saveRDS(tree_beta_no_site_model, file = "data/intermediate_output/vegetation/tree_beta_no_site_model.rda") #store model without site type
saveRDS(tree_beta_only_ap_sp, file = "data/intermediate_output/vegetation/tree_beta_only_ap_sp.rda") #store model with only AP and SP


#Model diagnostics
summary(tree_beta)
summary(tree_beta_stat)
summary(tree_beta_poly)
summary(tree_beta_no_interaction)
summary(tree_beta_lake)
summary(tree_beta_no_ap_cover)
summary(tree_beta_no_sp_cover)
summary(tree_beta_no_elevation)
summary(tree_beta_no_needle_share)
summary(tree_beta_no_tree_shannon)
summary(tree_beta_no_site_model)
summary(tree_beta_only_ap_sp)

AIC(tree_beta)
AIC(tree_beta_stat)
AIC(tree_beta_poly)
AIC(tree_beta_no_interaction)
AIC(tree_beta_lake)
AIC(tree_beta_no_ap_cover)
AIC(tree_beta_no_sp_cover)
AIC(tree_beta_no_elevation)
AIC(tree_beta_no_needle_share)
AIC(tree_beta_no_tree_shannon)
AIC(tree_beta_no_site_model)
AIC(tree_beta_only_ap_sp)

1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta, ~1))[1] - logLik(tree_beta)[1])) #Cox_snell R2
1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_stat, ~1))[1] - logLik(tree_beta_stat)[1])) #Cox_snell R2
1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_poly, ~1))[1] - logLik(tree_beta_poly)[1])) #Cox_snell R2
1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_interaction, ~1))[1] - logLik(tree_beta_no_interaction)[1])) #Cox_snell R2
1 - exp((2/nrow(dplyr::filter(ap_tree_cover_lake, tree >= 0))) * (logLik(update(tree_beta_lake, ~1))[1] - logLik(tree_beta_lake)[1])) #Cox_snell R2
1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_ap_cover, ~1))[1] - logLik(tree_beta_no_ap_cover)[1])) #Cox_snell R2
1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_sp_cover, ~1))[1] - logLik(tree_beta_no_sp_cover)[1])) #Cox_snell R2
1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_elevation, ~1))[1] - logLik(tree_beta_no_elevation)[1])) #Cox_snell R2
1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_needle_share, ~1))[1] - logLik(tree_beta_no_needle_share)[1])) #Cox_snell R2
1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_tree_shannon, ~1))[1] - logLik(tree_beta_no_tree_shannon)[1])) #Cox_snell R2
1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_site_model, ~1))[1] - logLik(tree_beta_no_site_model)[1])) #Cox_snell R2
1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_only_ap_sp, ~1))[1] - logLik(tree_beta_only_ap_sp)[1])) #Cox_snell R2


lmtest::lrtest(tree_beta)
lmtest::lrtest(tree_beta_stat,tree_beta) #LR test for variable precision

car::vif(tree_beta)
car::vif(tree_beta_no_interaction) #VIF values for no interaction

plot(resid(tree_beta) ~ fitted(tree_beta)) #investigate residuals
lines(lowess(resid(tree_beta) ~ fitted(tree_beta)))
dat.resid <- sum(resid(tree_beta, type = "pearson")^2)
1 - pchisq(dat.resid, tree_beta$df.resid) #check for fit
dat.resid/tree_beta$df.resid
link <-  tree_beta$link$mean

broom::glance(tree_beta) #model summary info
tree_beta$link$mean$linkinv(coef(tree_beta))
model_resid <- dplyr::as_tibble(residuals(tree_beta)) 

broom::augment(tree_beta)
broom::tidy(tree_beta)
broom::tidy(tree_beta, conf.int = TRUE)

plot(tree_beta, which = 1:4, type = "sweighted2")
plot(tree_beta, which = 5, type = "deviance", sub.caption = "")
plot(tree_beta, which = 1, type = "deviance", sub.caption = "")

#Model bootstraps and predictions
n_boots <- 1000 #Number of resamples
boot_ap_tree_cover <- list()
for (i in seq_len(n_boots)){ #resample mode training data
  boot_ap_tree_cover[[i]]<-ap_tree_cover[sample(seq_len(nrow(ap_tree_cover)), nrow(ap_tree_cover), replace=TRUE), ]
}

boot_tree_beta <- parallel::mclapply(boot_ap_tree_cover, function(x){ #reun betareg regressions on bootstraps
  tree_beta <- betareg::betareg(tree~ap_cover*elevation+needle_share+I(needle_share^2)+tree_shannon*elevation+I(tree_shannon^2)*elevation+site_model*elevation+sp_cover*elevation|needle_share+tree_shannon+site_model, data = dplyr::filter(x, tree > 0)) #final normal
}, mc.cores = 6)
saveRDS(boot_tree_beta, file = "data/intermediate_output/vegetation/boot_tree_beta.rda") #store bootstrapped models for downcore predictions

boot_dataset_pred <- lapply(boot_tree_beta, function(x){ #generate predictions based on model bootstraps and orginal data
  betareg::predict(x, ap_tree_cover_nona, type = "response") 
})

boot_dataset_pred_df <- dplyr::bind_rows(boot_dataset_pred, .id = "boot") %>%   #reshape data
  tidyr::pivot_longer(-boot, names_to = "ID", values_to = "prediction") %>% 
  dplyr::mutate(ID = as.numeric(ID)) %>% 
  dplyr::left_join(ap_tree_cover_nona, by = "ID")
saveRDS(boot_dataset_pred_df, file = "data/intermediate_output/vegetation/boot_dataset_pred_df.rda") #store bootstrapped predictions


# ---------------------------------------------------------


# 4. Leave one out cross validation
# ---------------------------------------------------------
# ---------------------------------------------------------

ap_tree_cover_nona <- ap_tree_cover %>% 
  dplyr::filter(tree > 0) %>% 
  dplyr::mutate(ID = dplyr::row_number()) 
dataset_test_ls <- split(ap_tree_cover_nona,1:nrow(ap_tree_cover_nona))

dataset_train_ls <- lapply(unique(ap_tree_cover_nona$ID), function(x){
  dplyr::filter(ap_tree_cover_nona, ID != x)
})
dataset_train_mod_ls <- parallel::mclapply(dataset_train_ls, function(x){
  betareg::betareg(tree~ap_cover*elevation+needle_share+I(needle_share^2)+tree_shannon*elevation+I(tree_shannon^2)*elevation+site_model*elevation+sp_cover*elevation|tree_shannon+site_model+needle_share, data = x,link = "logit") #Model choice
}, mc.cores = 8)
dataset_pred_ls <- purrr::map2(dataset_train_mod_ls,dataset_test_ls, betareg::predict, type = "response")
dataset_pred_df <- dplyr::bind_rows(dataset_pred_ls) %>% 
  dplyr::mutate(ID = dplyr::row_number())

dataset_pred_test <- ap_tree_cover_nona %>% #Calculate error
  dplyr::mutate(prediction = dataset_pred_df$`1`) %>% 
  dplyr::mutate(E = tree - prediction) %>% #observation minus prediction
  dplyr::mutate(SE = E^2) 

dataset_MAE <- mean(abs(dataset_pred_test$E))
dataset_MSE <- mean(dataset_pred_test$SE)
dataset_RMSE <- sqrt(dataset_MSE)
dataset_r2 <- cor(dataset_pred_test$tree, dataset_pred_test$prediction)^2

rio::export(dataset_pred_test, "data/intermediate_output/vegetation/loocv.csv")


ap_tree_cover_lake_nona <- ap_tree_cover_lake %>% 
  dplyr::filter(tree > 0) %>% 
  dplyr::mutate(ID = dplyr::row_number()) 
dataset_test_lake_ls <- split(ap_tree_cover_lake_nona,1:nrow(ap_tree_cover_lake_nona))

dataset_train_lake_ls <- lapply(unique(ap_tree_cover_lake_nona$ID), function(x){
  dplyr::filter(ap_tree_cover_lake_nona, ID != x)
})
dataset_train_mod_lake_ls <- parallel::mclapply(dataset_train_lake_ls, function(x){
  betareg::betareg(tree~ap_cover*elevation+needle_share+I(needle_share^2)+tree_shannon*elevation+I(tree_shannon^2)*elevation+sp_cover*elevation|tree_shannon+needle_share, data = x,link = "logit") #Model choice
}, mc.cores = 8)
dataset_pred_lake_ls <- purrr::map2(dataset_train_mod_lake_ls,dataset_test_lake_ls, betareg::predict, type = "response")
dataset_pred_lake_df <- dplyr::bind_rows(dataset_pred_lake_ls) %>% 
  dplyr::mutate(ID = dplyr::row_number())

dataset_pred_lake_test <- ap_tree_cover_lake_nona %>% #Calculate error
  dplyr::mutate(prediction = dataset_pred_lake_df$`1`) %>% 
  dplyr::mutate(E = tree - prediction) %>% #observation minus prediction
  dplyr::mutate(SE = E^2) 

dataset_lake_MAE <- mean(abs(dataset_pred_lake_test$E))
dataset_lake_MSE <- mean(dataset_pred_lake_test$SE)
dataset_lake_RMSE <- sqrt(dataset_lake_MSE)
dataset_lake_r2 <- cor(dataset_pred_lake_test$tree, dataset_pred_lake_test$prediction)^2


