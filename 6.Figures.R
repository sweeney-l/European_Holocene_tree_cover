#Sweeney, Harrison and Vander Linden 2024. 
#European forest cover during the Holocene reconstructed from pollen records
# ---------------------------------------------------------

#There are six seperate scripts associated with this research:
# 1. Map_setup.R 
# 2. Modern_pollen.R 
# 3. Model_import_data.R
# 4. Tree_model.R
# 5. Tree_recon.R 
# 6. Figures.R (this script)

# ---------------------------------------------------------

############ Figures ############ 

# Figures included within the paper and supplementary material


# ---------------------------------------------------------
# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Load in necessary data
# 3. Paper figures 
# 4. Supplementary figures

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
# dir.create("data/input/SMPDS")
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



# 2. Load in necessary data
# ---------------------------------------------------------
# ---------------------------------------------------------
ap_tree_cover <- rio::import("data/intermediate_output/vegetation/ap_tree_cover.csv")
ap_tree_plot <- rio::import("data/intermediate_output/vegetation/ap_tree_plot.csv")
bin_200_median_ap_tree_plot <- rio::import("data/intermediate_output/vegetation/bin_200_median_ap_tree_plot.csv")
bin_200_ap_tree_boots <- rio::import("data/intermediate_output/vegetation/bin_200_ap_tree_boots.csv")
bin_200_ap_tree_boots_5_95_median <- rio::import("data/intermediate_output/vegetation/bin_200_ap_tree_boots_5_95_median.csv")
bin_200_mean_recon_tree_plot <- rio::import("data/intermediate_output/vegetation/bin_200_mean_recon_tree_plot.csv")
bin_200_median_recon_tree_loc <- rio::import("data/intermediate_output/vegetation/bin_200_median_recon_tree_loc.csv")
bin_200_median_recon_tree_plot <- rio::import("data/intermediate_output/vegetation/bin_200_median_recon_tree_plot.csv")
bin_200_recon_tree_boots <- rio::import("data/intermediate_output/vegetation/bin_200_recon_tree_boots.csv")
bin_200_recon_tree_boots_5_95_mean <- rio::import("data/intermediate_output/vegetation/bin_200_recon_tree_boots_5_95_mean.csv")
bin_200_recon_tree_boots_5_95_median <- rio::import("data/intermediate_output/vegetation/bin_200_recon_tree_boots_5_95_median.csv")
bioregion_bin_200_median_recon_tree_loc_gr4 <- rio::import("data/intermediate_output/vegetation/bioregion_bin_200_median_recon_tree_loc_gr4.csv")
bioregion_recon_tree_binned_200 <- rio::import("data/intermediate_output/vegetation/bioregion_recon_tree_binned_200.csv")
comp_all_bioregiongr4_eqw_median_EU_loc <- rio::import("data/intermediate_output/vegetation/comp_all_bioregiongr4_eqw_median_EU_loc.csv")
comp_all_bioregiongr4_median_EU_loc <- rio::import("data/intermediate_output/vegetation/comp_all_bioregiongr4_median_EU_loc.csv")
comp_all_bioregiongr4_median_loc <- rio::import("data/intermediate_output/vegetation/comp_all_bioregiongr4_median_loc.csv")
comp_all_cell_median_loc <- rio::import("data/intermediate_output/vegetation/comp_all_cell_median_loc.csv")
comp_all_median_loc <- rio::import("data/intermediate_output/vegetation/comp_all_median_loc.csv")
comp_serge_median_all_locations_plot <- rio::import("data/intermediate_output/vegetation/comp_serge_median_all_locations.csv")
comp_zanon_median_all_locations_plot <- rio::import("data/intermediate_output/vegetation/comp_zanon_median_all_locations.csv")
cop_tree_cover_masked <- terra::rast("data/intermediate_output/vegetation/copernicus/cop_masked_tree_cover_100m.tif")
entity_bioregions <- rio::import("data/intermediate_output/vegetation/entity_bioregions.csv")
epd_bioregion_map_gr4 <- sf::read_sf("data/intermediate_output/vegetation/epd_bioregion_map_gr4.csv")
euro_map_3035 <- sf::st_read("data/intermediate_output/vegetation/euro_map/3035/euro_map_3035.shp")
euro_map_4258 <- sf::st_read("data/intermediate_output/vegetation/euro_map/4258/euro_map_4258.shp") %>% sf::st_set_crs(., 4258)
euro_map_4326 <- sf::st_read("data/intermediate_output/vegetation/euro_map/4326/euro_map_4326.shp")  %>% sf::st_set_crs(., 4326)
fitted_tree_plot <- rio::import("data/intermediate_output/vegetation/fitted_tree_plot.csv")
fitted_tree_plot_qmap <- rio::import("data/intermediate_output/vegetation/fitted_tree_plot_qmap.csv")
lq_bin_200_median_recon_tree <- rio::import("data/intermediate_output/vegetation/lq_bin_200_median_recon_tree.csv")
modern_comp_serge <- rio::import("data/intermediate_output/vegetation/modern_comp_serge.csv")
modern_comp_sergeID <- rio::import("data/intermediate_output/vegetation/modern_comp_sergeID.csv")
modern_comp_zan <- rio::import("data/intermediate_output/vegetation/modern_comp_zan.csv")
modern_comp_zanID <- rio::import("data/intermediate_output/vegetation/modern_comp_zanID.csv")
peak_range_comp_all_bioregiongr4_median_loc <- rio::import("data/intermediate_output/vegetation/peak_range_comp_all_bioregiongr4_median_loc.csv")
recon_tree_locations <- rio::import("data/intermediate_output/vegetation/recon_tree.csv")
recon_tree_rast_200_class_plot <-  readRDS("data/intermediate_output/vegetation/raster_tree/recon_tree_rast_200_class_plot.rda")
refitted_model_df <- rio::import("data/intermediate_output/vegetation/refitted_model_df.csv")
taxa_cat <- rio::import("data/input/taxa_cat.csv")
tree_beta <- readRDS("data/intermediate_output/vegetation/tree_beta.rda") #modern tree cover model
tree_beta_lake <- readRDS("data/intermediate_output/vegetation/tree_beta_lake.rda") #modern tree cover model
tree_beta_no_ap_cover <-  readRDS("data/intermediate_output/vegetation/tree_beta_no_ap_cover.rda") #modern tree cover model without ap
tree_beta_no_sp_cover <- readRDS("data/intermediate_output/vegetation/tree_beta_no_sp_cover.rda") #modern tree cover model without sp
tree_beta_no_elevation <- readRDS("data/intermediate_output/vegetation/tree_beta_no_elevation.rda") #modern tree cover model without elevation
tree_beta_no_needle_share <- readRDS("data/intermediate_output/vegetation/tree_beta_no_needle_share.rda") #modern tree cover model without needle share
tree_beta_no_tree_shannon <- readRDS("data/intermediate_output/vegetation/tree_beta_no_tree_shannon.rda") #modern tree cover model without tree shannon
tree_beta_no_site_model <- readRDS("data/intermediate_output/vegetation/tree_beta_no_site_model.rda") #modern tree cover model without site type
tree_beta_only_ap_sp <- readRDS("data/intermediate_output/vegetation/tree_beta_only_ap_sp.rda") #modern tree cover model without only AP and SP
uq_bin_200_median_recon_tree <- rio::import("data/intermediate_output/vegetation/uq_bin_200_median_recon_tree.csv")


# ---------------------------------------------------------




# 3. Paper figures
# ---------------------------------------------------------
# ---------------------------------------------------------

#Figure 1: Methodology to reconstruct European tree cover during the Holocene
#Manually constructed


# Figure 2: A - Observed tree cover based on compositing annual tree cover maps from the Copernicus land cover data sets (2015-2019) and screening out cells where the dominant land cover was not natural; 
#           B - Modern pollen records used for model fitting; 
#           C - Fossil pollen sites used for tree cover reconstructions; 
#           D -  Classification of the fossil pollen sites into climatic sub-regions

cop_masked_tree_cover_plot <- ggplotify::as.ggplot(~terra::plot(cop_tree_cover_masked, 
            col=grDevices::hcl.colors(50, palette = "viridis", rev = TRUE),
            axes = FALSE,
            mar = c(1.5,4,1.5,6),
            plg=list(title = "Tree cover\n%", title.cex = 0.55, cex = 0.55, size = c(0.75,1)))) # Legend text size

model_locations_plot <- ggplot2::ggplot(data = euro_map_4326)+
  geom_sf(fill = "seashell") +
  theme(panel.background = element_rect(fill = "aliceblue"))+
  geom_point(data = ap_tree_cover, mapping = aes(x = longitude, y = latitude), size = .5)+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())+
  xlab(NULL)+
  ylab(NULL)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  scale_color_discrete(name = "Site type")

recon_tree_locations2 <- recon_tree_locations %>% 
  dplyr::select(entity_name, latitude, longitude) %>% 
  dplyr::distinct() %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), agr = "constant",crs = 4326) %>% 
  sf::st_coordinates()
  
downcore_locations_plot<- ggplot2::ggplot(data = euro_map_4326)+
  geom_sf(fill = "seashell") +
  theme(panel.background = element_rect(fill = "aliceblue"))+
  geom_point(data = recon_tree_locations2, mapping = aes(x = X, y = Y), size = .5, color = "red")+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  xlab(NULL)+
  ylab(NULL)

group.colors <- c(Alpine = "#F8766D", 
                  Anatolian = "#D89000",
                  Arctic = "#FF62BC",
                  Atlantic = "#39B600",
                  BlackSea = "#00BF7D",
                  Boreal = "#00BFC4", 
                  Continental ="#A3A500",
                  Mediterranean = "#9590FF",
                  Pannonian = "#E76BF3",
                  Steppic = "#00B0F6")

epd_bioregion_map_gr4_1 <- epd_bioregion_map_gr4 %>% 
  dplyr::mutate(longitude = as.numeric(longitude), latitude = as.numeric(latitude))

gr4_bioregion_downcore_plot <- ggplot2::ggplot(data = euro_map_4326)+
  geom_sf(fill = "seashell") +
  theme(panel.background = element_rect(fill = "aliceblue"))+
  geom_point(epd_bioregion_map_gr4_1, mapping = aes(x = longitude, y = latitude, color = bioregion_name), size = 0.5)+
  labs(color = "Bioregion name")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  scale_color_manual(values=group.colors)+
  theme(legend.text=element_text(size=7))+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())+
  xlab(NULL)+
  ylab(NULL)+
  theme(legend.position = c(0.80, 0.90), legend.title = element_blank(), legend.key.size = unit(0.2, "cm"), legend.key = element_rect(color = NA, fill = NA))

ggpubr::ggarrange(cop_masked_tree_cover_plot,model_locations_plot,
                  downcore_locations_plot,gr4_bioregion_downcore_plot, 
                  nrow = 2, ncol = 2, 
                  labels = c("A", "B", "C", "D"), hjust = -1.5, vjust = 3, font.label = list(family = "Times New Roman"))
ggsave("figs/2.maps_tree_modern_downcore_bioregion.png",height = 12, width = 16, unit = "cm")



# Table 1: Modern tree cover model coefficients
results_model <- summary(tree_beta)
results_table_mean <- dplyr::as_tibble(results_model[["coefficients"]][["mean"]], rownames = "coefficients")
results_table_precision <- dplyr::as_tibble(results_model[["coefficients"]][["precision"]], rownames = "coefficients")
rio::export(list(results_table_mean,results_table_precision), "figs/1.model_results_table.xlsx")

# Table 2: Change in modern model AIC values when excluding specific variables (exclusion includes interactions, polynomials and precision variables)
modern_AIC_values <- dplyr::tibble(Model = c("Final model", 
                                             "       excluding AP",
                                             "       excluding SP",
                                             "       excluding %needleshare",
                                             "       excluding AP Shannon index",
                                             "       excluding lake or bog site",
                                             "       excluding elevation",
                                             "AP and SP model"),
                                   AIC_value = c(AIC(tree_beta),
                                                     AIC(tree_beta_no_ap_cover),
                                                     AIC(tree_beta_no_sp_cover),
                                                     AIC(tree_beta_no_needle_share),
                                                     AIC(tree_beta_no_tree_shannon),
                                                     AIC(tree_beta_no_site_model),
                                                     AIC(tree_beta_no_elevation),
                                                 AIC(tree_beta_only_ap_sp)),
                                   Cox_Snell = c((1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta, ~1))[1] - logLik(tree_beta)[1]))), 
                                                 (1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_ap_cover, ~1))[1] - logLik(tree_beta_no_ap_cover)[1]))), 
                                                 (1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_sp_cover, ~1))[1] - logLik(tree_beta_no_sp_cover)[1]))), 
                                                 (1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_needle_share, ~1))[1] - logLik(tree_beta_no_needle_share)[1]))), 
                                                 (1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_tree_shannon, ~1))[1] - logLik(tree_beta_no_tree_shannon)[1]))), 
                                                 (1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_site_model, ~1))[1] - logLik(tree_beta_no_site_model)[1]))), 
                                                 (1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_no_elevation, ~1))[1] - logLik(tree_beta_no_elevation)[1]))), 
                                                 (1 - exp((2/nrow(dplyr::filter(ap_tree_cover, tree >= 0))) * (logLik(update(tree_beta_only_ap_sp, ~1))[1] - logLik(tree_beta_only_ap_sp)[1]))))) 
modern_AIC_values_change <- modern_AIC_values %>% 
  dplyr::mutate(Change_in_AIC = AIC_value - as.numeric(modern_AIC_values[1,2]))
rio::export(modern_AIC_values_change, "figs/2.model_AIC_table.xlsx")



# Figure 3: Evaluation of model performance. 
#           A – Differences between predictions and observations (residual), in bins of observed tree cover percentage; 
#           B – Predictions of tree cover compared to observed tree cover
fitted_tree_plot_qmap_order <- fitted_tree_plot_qmap %>%
  dplyr::arrange(observed_group) %>% 
  dplyr::mutate(observed_group = factor(observed_group, level = c("[0,10]","(10,20]","(20,30]","(30,40]","(40,50]","(50,60]", "(60,70]","(70,80]","(80,90]","(90,100]")))

refit_obs_scatter_plot <- ggplot2::ggplot(data = tidyr::drop_na(fitted_tree_plot_qmap), mapping = aes(x = refitted100, y = observed100))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Predicted tree cover (%)")+
  ylab("Observed tree cover (%)")+
  scale_y_continuous(limit = c(-0,100))+
  scale_x_continuous(limit = c(-0,100))+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))

refit_obs_plot <- ggplot2::ggplot(data = tidyr::drop_na(fitted_tree_plot_qmap_order), mapping = aes(y = difference2*100, x = observed_group))+
  geom_boxplot()+
  theme_bw()+
  scale_y_continuous(limit = c(-75,75))+
  scale_x_discrete(labels = c("[0,10]"="<10","(10,20]"="10-19","(20,30]"="20-29","(30,40]"="30-39","(40,50]"="40-49",
                              "(50,60]"="50-59", "(60,70]"="60-69","(70,80]"="70-79","(80,90]"="80-89","(90,100]"="90-100"),
                   guide = guide_axis(n.dodge = 2))+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Observed tree cover (%)")+
  ylab("Residual (%)")+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))


ggpubr::ggarrange(refit_obs_plot, refit_obs_scatter_plot, nrow = 1,
                  labels = c("A", "B"), hjust = -4.5, vjust = 3, font.label = list(family = "Times New Roman"))
ggsave("figs/3.predicted_refitted_obsevations_scatter_qmap.png", height = 6, width = 16, unit = "cm")





# Figure 4: Modern tree cover from Zanon et al (2018) compared to 
#                                                                 (A) observed tree cover and 
#                                                                 (B) our predicted tree cover. 
#           Modern tree cover from Serge et al. (2023) compared to 
#                                                                 (C) observed tree cover values and 
#                                                                 (D) our predicted tree cover


modern_obs_zanon_plot <- ggplot2::ggplot(data = modern_comp_zan, mapping = aes(x = zan_tree, y = tree))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Zanon tree cover (%)")+
  ylab("Observed tree cover (%)")+
  theme(axis.text=element_text(size=9),axis.title=element_text(size=9))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))

modern_obs_serge_plot <- ggplot2::ggplot(data = modern_comp_serge, mapping = aes(x = serge_tree, y = tree))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Serge tree cover (%)")+
  ylab("Observed tree cover (%)")+
  theme(axis.text=element_text(size=9),axis.title=element_text(size=9))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))

modern_pred_zanon_plot <- ggplot2::ggplot(data = modern_comp_zan, mapping = aes(x = zan_tree, y = refitted))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Zanon tree cover (%)")+
  ylab("Predicted tree cover (%)")+
  theme(axis.text=element_text(size=9),axis.title=element_text(size=9))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))

modern_pred_serge_plot <- ggplot2::ggplot(data = modern_comp_serge, mapping = aes(x = serge_tree, y = refitted))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Serge tree cover (%)")+
  ylab("Predicted tree cover (%)")+
  theme(axis.text=element_text(size=9),axis.title=element_text(size=9))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))

ggpubr::ggarrange(modern_obs_zanon_plot,modern_pred_zanon_plot,
                  modern_obs_serge_plot,modern_pred_serge_plot, 
                  nrow = 2, ncol = 2,
                  labels = c("A", "B", "C", "D"), hjust = -4.5, vjust = 3, font.label = list(family = "Times New Roman"))
ggsave("figs/4.modern_comparisons.png",height = 12, width = 16, unit = "cm")




# Figure 5: A - Median reconstructed tree cover for Europe from 12,000 to 0 cal. BP, with 95% confidence intervals for 1000 bootstrap resampling of records; 
#           B - Median reconstructed tree cover for Europe from 12,000 to 0 cal. BP, with differing LOESS regression smoothing half-widths

plot_binned_median_recon_boots <- ggplot2::ggplot(data = dplyr::filter(bin_200_recon_tree_boots, bin_age <= 12000), 
                mapping = aes(y = tree_cover_median, x = bin_age, group = column_label))+
  geom_line(color = "light grey", size = 0.3)+
  geom_line(data = dplyr::filter(bin_200_median_recon_tree_plot, bin_age <= 12000),
            mapping = aes(y = tree_cover_median, x = bin_age), color = "yellow", size = 1)+
  geom_line(data = dplyr::filter(bin_200_recon_tree_boots_5_95_median, bin_age <= 12000),
            mapping = aes(y = tree_cover_5_95, x = bin_age, group = lower_upper), color = "red",size = 0.2)+
  ylim(25,75)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  scale_x_reverse(breaks = seq(12000, 0, -2000))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))

plot_bin_200_median_recon_tree_loc <- ggplot2::ggplot(data = dplyr::filter(bin_200_median_recon_tree_loc, bin_age <= 12000 & half_width %in% c("None", "500-year")), 
                                                      mapping = aes(x = bin_age, y = median, color = half_width))+
  geom_line()+
  ylim(25,75)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  scale_color_discrete(name = "Smoothing half-width", 
                       breaks = c("None", "500-year"))+
  scale_x_reverse(breaks = seq(12000, 0, -2000))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"), 
        legend.title=element_text(size=9), 
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.2, "cm"),
        legend.position = c(0.7, 0.15),
        legend.direction = "vertical")+
  guides(color = guide_legend(nrow = 2))

ggpubr::ggarrange(plot_binned_median_recon_boots,plot_bin_200_median_recon_tree_loc,
                  nrow = 1, ncol = 2,
                  labels = c("A", "B", "C"), hjust = -4.5, vjust = 3, font.label = list(family = "Times New Roman"))
ggsave("figs/5.median_tree_cover_plots.png",height = 6, width = 16, unit = "cm")


# Figure 6: Gridded maps of average reconstructed tree cover for selected periods, for 50km2 grid cells. Bin ages are 200-years in width, with ages referring to mid-point of each bin

recon_tree_rast_200_overview <- gridExtra::grid.arrange(recon_tree_rast_200_class_plot$`11000`, 
                                                        recon_tree_rast_200_class_plot$`9000`,
                                                        recon_tree_rast_200_class_plot$`6000`,
                                                        recon_tree_rast_200_class_plot$`2000`, ncol = 2)
ggsave("figs/6.recon_tree_rast_200_overview.png",recon_tree_rast_200_overview, height = 12, width = 16, unit = "cm")
  

# Figure 7: Median tree cover values, for selected biogeographical regions. Smoothed lines reflect LOESS fitted regression with 1000-year halfwidth

ggplot2::ggplot(data = dplyr::filter(bioregion_bin_200_median_recon_tree_loc_gr4, half_width == "None"), mapping = aes(y = average, x = bin_age, color = bioregion_name))+
  geom_line()+
  geom_line(data = dplyr::filter(bioregion_bin_200_median_recon_tree_loc_gr4, half_width == "1000-year"), mapping = aes(y = average, x = bin_age, color = bioregion_name), size = 1)+
  ylim(10,80)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  labs(color = "Bioregion")+
  scale_color_manual(values=group.colors)+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10), legend.title=element_text(size=9), legend.text=element_text(size=8), legend.box.spacing = unit(0, "pt"))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_x_reverse(limits = c(12000, 0), breaks = seq(12000, 0, -2000))
ggsave("figs/7.binned_median_recon_bioregion_gr4.png", height = 6, width = 10, unit = "cm")



# Figure 8: Reconstructed median tree cover compared to equivalent extracted tree cover medians for Serge et al. (2023) and Zanon et al (2018). Smoothed lines reflect LOESS fitted regression with 1000-year halfwidth 

ggplot2::ggplot(data = dplyr::filter(comp_all_median_loc, half_width == "None" & source != "githumbi"), 
                mapping = aes(y = average, x = bin_age, color = source))+
  geom_line()+
  geom_line(data = dplyr::filter(comp_all_median_loc, half_width == "1000-year" & source != "githumbi"), 
            mapping = aes(y = average, x = bin_age, color = source), size =1)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  ggthemes::scale_color_colorblind(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  # scale_color_discrete(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10), legend.title=element_text(size=9), legend.text=element_text(size=8))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_x_reverse(breaks = seq(12000, 0, -2000))
ggsave("figs/8.comp_all_median_loc.png", height = 6, width = 10, unit = "cm")



# Figure 9: Reconstructed median tree cover compared to equivalent extracted tree cover medians for Serge et al. (2023) and Zanon et al (2018), for selected modern biogeographical regions. Smoothed lines reflect LOESS fitted regression with 1000-year halfwidth. 
#           A – Atlantic; 
#           B- Boreal; 
#           C- Continental; 
#           D – Mediterranean

plot_comp_all_median_loc_atlantic <- ggplot2::ggplot(data = dplyr::filter(comp_all_bioregiongr4_median_loc, 
                                     bioregion_name == "Atlantic" & 
                                       half_width == "None" &
                                       source != "githumbi"), 
                mapping = aes(y = average, x = bin_age, color = source))+
  geom_line()+
  geom_line(data = dplyr::filter(comp_all_bioregiongr4_median_loc, 
                                 bioregion_name == "Atlantic" & 
                                   half_width == "1000-year" &
                                   source != "githumbi"), 
            mapping = aes(y = average, x = bin_age, color = source), size =1)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  ggthemes::scale_color_colorblind(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  # scale_color_discrete(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10), legend.title=element_text(size=9), legend.text=element_text(size=8))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_x_reverse(breaks = seq(12000, 0, -2000))

plot_comp_all_median_loc_boreal <- ggplot2::ggplot(data = dplyr::filter(comp_all_bioregiongr4_median_loc, 
                                     bioregion_name == "Boreal" & 
                                       half_width == "None" &
                                       source != "githumbi"), 
                mapping = aes(y = average, x = bin_age, color = source))+
  geom_line()+
  geom_line(data = dplyr::filter(comp_all_bioregiongr4_median_loc, 
                                 bioregion_name == "Boreal" & 
                                   half_width == "1000-year" &
                                   source != "githumbi"), 
            mapping = aes(y = average, x = bin_age, color = source), size =1)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  ggthemes::scale_color_colorblind(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  # scale_color_discrete(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10), legend.title=element_text(size=9), legend.text=element_text(size=8))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_x_reverse(breaks = seq(12000, 0, -2000))

plot_comp_all_median_loc_continental <- ggplot2::ggplot(data = dplyr::filter(comp_all_bioregiongr4_median_loc, 
                                     bioregion_name == "Continental" & 
                                       half_width == "None" &
                                       source != "githumbi"), 
                mapping = aes(y = average, x = bin_age, color = source))+
  geom_line()+
  geom_line(data = dplyr::filter(comp_all_bioregiongr4_median_loc, 
                                 bioregion_name == "Continental" & 
                                   half_width == "1000-year" &
                                   source != "githumbi"), 
            mapping = aes(y = average, x = bin_age, color = source), size =1)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  ggthemes::scale_color_colorblind(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  # scale_color_discrete(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10), legend.title=element_text(size=9), legend.text=element_text(size=8))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_x_reverse(breaks = seq(12000, 0, -2000))

plot_comp_all_median_loc_mediterranean <- ggplot2::ggplot(data = dplyr::filter(comp_all_bioregiongr4_median_loc, 
                                     bioregion_name == "Mediterranean" & 
                                       half_width == "None" &
                                       source != "githumbi"), 
                mapping = aes(y = average, x = bin_age, color = source))+
  geom_line()+
  geom_line(data = dplyr::filter(comp_all_bioregiongr4_median_loc, 
                                 bioregion_name == "Mediterranean" & 
                                   half_width == "1000-year" &
                                   source != "githumbi"), 
            mapping = aes(y = average, x = bin_age, color = source), size =1)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  ggthemes::scale_color_colorblind(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  # scale_color_discrete(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10), legend.title=element_text(size=9), legend.text=element_text(size=8))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_x_reverse(breaks = seq(12000, 0, -2000))

ggpubr::ggarrange(plot_comp_all_median_loc_atlantic,plot_comp_all_median_loc_boreal,
                  plot_comp_all_median_loc_continental,plot_comp_all_median_loc_mediterranean, 
                  nrow = 2, ncol = 2,
                  common.legend = TRUE, legend="bottom",
                  labels = c("A", "B", "C", "D"), hjust = -4.5, vjust = 3, font.label = list(family = "Times New Roman"))
ggsave("figs/9.bioregion_comparisons_median.png",height = 12, width = 16, unit = "cm")






# ---------------------------------------------------------



# 4. Supplementary figures
# ---------------------------------------------------------
# ---------------------------------------------------------

# Supplementary Table 1: Sources for SMPDSv2 data
#Manually created

# Supplementary Table 2: Higher level grouping of taxa considered native to Europe and included within the Total Terrestrial Pollen Sum
taxa_cat_single <- taxa_cat %>% #for single clean_taxon_name
  dplyr::select(-taxon_name) %>% 
  dplyr::distinct()
taxa_cat_tps_native <- taxa_cat_single %>% 
  dplyr::filter(terrestrial_pollen_sum == "yes") %>% 
  dplyr::filter(europe == "native") %>% 
  dplyr::select(-c(terrestrial_pollen_sum,europe)) 
taxa_cat_broad <- taxa_cat_tps_native %>% 
  dplyr::filter(ap_needle_broad == "broad") %>% 
  dplyr::arrange(clean_taxon_name)
taxa_cat_needle <- taxa_cat_tps_native %>% 
  dplyr::filter(ap_needle_broad == "needle") %>% 
  dplyr::arrange(clean_taxon_name)
taxa_cat_shrub <- taxa_cat_tps_native %>% 
  dplyr::filter(ap_sp_hp == "SP") %>% 
  dplyr::arrange(clean_taxon_name)
taxa_cat_hp <- taxa_cat_tps_native %>% 
  dplyr::filter(ap_sp_hp == "HP") %>% 
  dplyr::arrange(clean_taxon_name)
taxa_cat_table <- tibble("Tree/Shrub/Non-arboreal" = c("Tree","Tree",
                                                       "Shrub","Shrub",
                                                       "Non-arboreal","Non-arboreal"),
                         "Tree sub-division" = c("Broadleaf", "Needleleaf",
                                                 NA,NA,
                                                 NA,NA),
                         "Taxa grouping" = c(toString(c(taxa_cat_broad$clean_taxon_name)),toString(c(taxa_cat_needle$clean_taxon_name)),
                                             toString(c(taxa_cat_shrub$clean_taxon_name)),NA,
                                             toString(c(taxa_cat_hp$clean_taxon_name)),NA))
rio::export(taxa_cat_table,"figs/supplement/ST2.taxa_cat_table.xlsx")


# Supplementary Table 3: Modern tree cover model coefficients: Lake sites
results_model_lake <- summary(tree_beta_lake)
results_table_lake_mean <- dplyr::as_tibble(results_model_lake[["coefficients"]][["mean"]], rownames = "coefficients")
results_table_lake_precision <- dplyr::as_tibble(results_model_lake[["coefficients"]][["precision"]], rownames = "coefficients")
rio::export(list(results_table_lake_mean,results_table_lake_precision), "figs/supplement/ST3.model_results_lake_table.xlsx")


# Supplementary Table 4: Modern tree cover model coefficients: Including higher elevation sites* (n.b. run tree beta without elevation filter)
results_model <- summary(tree_beta)
results_table_mean <- dplyr::as_tibble(results_model[["coefficients"]][["mean"]], rownames = "coefficients")
results_table_precision <- dplyr::as_tibble(results_model[["coefficients"]][["precision"]], rownames = "coefficients")
rio::export(list(results_table_mean,results_table_precision), "figs/supplement/ST4.model_results_table.xlsx")

# Supplementary Table 4: Modern tree cover model coefficients: Stricter 500m limit* (n.b. run tree beta with 500m elevation filter)
results_model <- summary(tree_beta)
results_table_mean <- dplyr::as_tibble(results_model[["coefficients"]][["mean"]], rownames = "coefficients")
results_table_precision <- dplyr::as_tibble(results_model[["coefficients"]][["precision"]], rownames = "coefficients")
rio::export(list(results_table_mean,results_table_precision), "figs/supplement/ST5.model_results_table.xlsx")


# Supplementary Figure 1: Model performance: 
#                                           A – Non-adjusted predictions of tree cover compared to observed tree cover; 
#                                           B - Predictions of tree cover compared with observed tree cover; 
#                                           C - Differences between non-adjusted predictions and observations (residual), in bins of observed tree cover percentage; 
#                                           D - Differences between predictions and observations (residual), in bins of observed tree cover percentage

fitted_tree_plot_qmap_order <- fitted_tree_plot_qmap %>%
  dplyr::arrange(observed_group) %>% 
  dplyr::mutate(observed_group = factor(observed_group, level = c("[0,10]","(10,20]","(20,30]","(30,40]","(40,50]","(50,60]", "(60,70]","(70,80]","(80,90]","(90,100]")))

fitted_tree_plot_order <- fitted_tree_plot %>%
  dplyr::arrange(observed_group) %>% 
  dplyr::mutate(observed_group = factor(observed_group, level = c("[0,10]","(10,20]","(20,30]","(30,40]","(40,50]","(50,60]", "(60,70]","(70,80]","(80,90]","(90,100]")))

pred_obs_scatter_plot <- ggplot2::ggplot(data = tidyr::drop_na(fitted_tree_plot), mapping = aes(x = predicted100, y = observed100))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Non-adjusted predicted tree cover (%)")+
  ylab("Observed tree cover (%)")+
  scale_y_continuous(limit = c(-0,100))+
  scale_x_continuous(limit = c(-0,100))+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))

refit_obs_scatter_plot <- ggplot2::ggplot(data = tidyr::drop_na(fitted_tree_plot_qmap), mapping = aes(x = refitted100, y = observed100))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Predicted tree cover (%)")+
  ylab("Observed tree cover (%)")+
  scale_y_continuous(limit = c(-0,100))+
  scale_x_continuous(limit = c(-0,100))+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))

pred_obs_plot <- ggplot2::ggplot(data = tidyr::drop_na(fitted_tree_plot_order), mapping = aes(y = difference*100, x = observed_group))+
  geom_boxplot()+
  theme_bw()+
  scale_y_continuous(limit = c(-75,75))+
  scale_x_discrete(labels = c("[0,10]"="<10","(10,20]"="10-19","(20,30]"="20-29","(30,40]"="30-39","(40,50]"="40-49",
                              "(50,60]"="50-59", "(60,70]"="60-69","(70,80]"="70-79","(80,90]"="80-89","(90,100]"="90-100"),
                   guide = guide_axis(n.dodge = 2))+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Observed tree cover (%)")+
  ylab("Residual (%)")+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))

refit_obs_plot <- ggplot2::ggplot(data = tidyr::drop_na(fitted_tree_plot_qmap_order), mapping = aes(y = difference2*100, x = observed_group))+
  geom_boxplot()+
  theme_bw()+
  scale_y_continuous(limit = c(-75,75))+
  scale_x_discrete(labels = c("[0,10]"="<10","(10,20]"="10-19","(20,30]"="20-29","(30,40]"="30-39","(40,50]"="40-49",
                              "(50,60]"="50-59", "(60,70]"="60-69","(70,80]"="70-79","(80,90]"="80-89","(90,100]"="90-100"),
                   guide = guide_axis(n.dodge = 2))+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Observed tree cover (%)")+
  ylab("Residual (%)")+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))

ggpubr::ggarrange(pred_obs_scatter_plot, refit_obs_scatter_plot, pred_obs_plot, refit_obs_plot, nrow = 2, ncol = 2,
                  labels = c("A", "B", "C", "D"), hjust = -4.5, vjust = 3, font.label = list(family = "Times New Roman"))
ggsave("figs/supplement/SF1.predicted_refitted_obsevations_scatter_qmap.png", height = 12, width = 16, unit = "cm")



# Supplementary Figure 2: Spatial structure of differences between quantile mapped adjusted predictions and observations
ggplot(data = euro_map_3035)+
  geom_sf()+
  geom_sf(fill = "seashell") +
  theme(panel.background = element_rect(fill = "aliceblue"))+
  geom_tile(refitted_model_df, mapping = aes(x = x, y = y, fill = Percentage_difference))+
  labs(fill="Percentage difference")+
  ggthemes::scale_fill_colorblind()+
  theme(axis.title.y=element_blank(),  
        axis.title.x=element_blank())+
  # theme(legend.position = c(.15,.85))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))+
  theme(legend.title = element_text(size=10))
ggsave("figs/supplement/SF2.map_differences_model.png",height = 10, width = 16, unit = "cm")


# Supplementary Figure 3: Arboreal pollen percentage compared to observed tree cover: 
#                                                                                   A - AP% compared to observed tree cover for each record;
#                                                                                   B - Differences between AP% and observations (residual), in bins of observed tree cover percentage
ap_tree_plot_order <- ap_tree_plot %>%
  dplyr::arrange(observed_group) %>% 
  dplyr::mutate(observed_group = factor(observed_group, level = c("[0,10]","(10,20]","(20,30]","(30,40]","(40,50]","(50,60]", "(60,70]","(70,80]","(80,90]","(90,100]")))

ap_obs_scatter_plot <- ggplot2::ggplot(data = tidyr::drop_na(ap_tree_plot), mapping = aes(x = ap100, y = observed100))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Arboreal pollen (%)")+
  ylab("Observed tree cover (%)")+
  scale_y_continuous(limit = c(-0,100))+
  scale_x_continuous(limit = c(-0,100))+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))

ap_obs_plot <- ggplot2::ggplot(data = tidyr::drop_na(ap_tree_plot_order), mapping = aes(y = difference*100, x = observed_group))+
  geom_boxplot()+
  theme_bw()+
  scale_y_continuous(limit = c(-75,75))+
  scale_x_discrete(labels = c("[0,10]"="<10","(10,20]"="10-19","(20,30]"="20-29","(30,40]"="30-39","(40,50]"="40-49",
                              "(50,60]"="50-59", "(60,70]"="60-69","(70,80]"="70-79","(80,90]"="80-89","(90,100]"="90-100"),
                   guide = guide_axis(n.dodge = 2))+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Observed tree cover (%)")+
  ylab("Residual (%)")+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))

ggpubr::ggarrange(ap_obs_scatter_plot, ap_obs_plot, nrow = 1, ncol = 2,
                  labels = c("A", "B", "C", "D"), hjust = -4.5, vjust = 3, font.label = list(family = "Times New Roman"))
ggsave("figs/supplement/SF3.ap_obsevations_scatter.png", height = 6, width = 16, unit = "cm")


# Supplementary Figure 4: Extracted tree cover values averaged per gridcell for the first bin of Zanon et al (2018) compared to 
#                             (A) observed tree cover values and 
#                             (B) predicted model tree cover values. 
#                         Extracted tree cover values averaged per gridcell for the first bin of Serge et al. (2023) compared to 
#                             (C) observed tree cover values and 
#                             (D) predicted model tree cover values

modern_obs_zanon_plot_ID <- ggplot2::ggplot(data = modern_comp_zanID, mapping = aes(x = zan_tree, y = tree))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Zanon tree cover (%)")+
  ylab("Observed tree cover (%)")+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))

modern_obs_serge_plot_ID <- ggplot2::ggplot(data = modern_comp_sergeID, mapping = aes(x = serge_tree, y = tree))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Serge tree cover (%)")+
  ylab("Observed tree cover (%)")+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))

modern_pred_zanon_plot_ID <- ggplot2::ggplot(data = modern_comp_zanID, mapping = aes(x = zan_tree, y = refitted))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Zanon tree cover (%)")+
  ylab("Predicted tree cover (%)")+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))

modern_pred_serge_plot_ID <- ggplot2::ggplot(data = modern_comp_sergeID, mapping = aes(x = serge_tree, y = refitted))+
  geom_point(size = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  xlab("Serge tree cover (%)")+
  ylab("Predicted tree cover (%)")+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))

ggpubr::ggarrange(modern_obs_zanon_plot_ID,modern_pred_zanon_plot_ID,
                  modern_obs_serge_plot_ID,modern_pred_serge_plot_ID, 
                  nrow = 2, ncol = 2,
                  labels = c("A", "B", "C", "D"), hjust = -4.5, vjust = 3, font.label = list(family = "Times New Roman"))
ggsave("figs/supplement/SF4.modern_comparisons_ID.png",height = 12, width = 16, unit = "cm")



# Supplementary Figure 5: A - Mean reconstructed tree cover for Europe from 12,000 to 0 cal. BP, with 95% confidence intervals for 1000 bootstrap resamples of records; 
#                         B - Mean reconstructed tree cover for Europe from 12,000 to 0 cal. BP, with differing LOESS regression smoothing half-widths applied

plot_binned_mean_recon_boots<- ggplot2::ggplot(data = dplyr::filter(bin_200_recon_tree_boots, bin_age <= 12000), 
                                               mapping = aes(y = tree_cover_mean, x = bin_age, group = column_label))+
  geom_line(color = "light grey", size = 0.3)+
  geom_line(data = dplyr::filter(bin_200_mean_recon_tree_plot, bin_age <= 12000),
            mapping = aes(y = tree_cover_mean, x = bin_age), color = "yellow", size = 1)+
  geom_line(data = dplyr::filter(bin_200_recon_tree_boots_5_95_mean, bin_age <= 12000),
            mapping = aes(y = tree_cover_5_95, x = bin_age, group = lower_upper), color = "red",size = 0.2)+
  ylim(25,75)+
  theme_bw()+
  ylab("Mean tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  scale_x_reverse(breaks = seq(12000, 0, -2000))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))

plot_bin_200_mean_recon_tree_loc <- ggplot2::ggplot(data = filter(bin_200_mean_recon_tree_loc, bin_age <= 12000  & half_width %in% c("None", "500-year")),
                                                    mapping = aes(x = bin_age, y = mean, color = half_width))+
  geom_line()+
  ylim(25,75)+
  theme_bw()+
  ylab("Mean tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  scale_color_discrete(name = "Smoothing half-width", 
                       breaks = c("None", "500-year"))+
  scale_x_reverse(breaks = seq(12000, 0, -2000))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"), 
        legend.title=element_text(size=9), 
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.2, "cm"),
        legend.position = c(0.7, 0.15),
        legend.direction = "vertical")+
  guides(color = guide_legend(nrow = 2))

ggpubr::ggarrange(plot_binned_mean_recon_boots,plot_bin_200_mean_recon_tree_loc,
                  nrow = 1, ncol = 2,
                  labels = c("A", "B"), hjust = -4.5, vjust = 3, font.label = list(family = "Times New Roman"))
ggsave("figs/supplement/SF5.mean_tree_cover_plots.png",height = 6, width = 16, unit = "cm")

# Supplementary Figure 6: Median reconstructed tree cover, with model bootstraps
plot_boot_binned_median <- ggplot2::ggplot(data = dplyr::filter(boot_bin_200_median_recon_tree, bin_age <= 12000), 
                                           mapping = aes(y = median_tree_cover, x = bin_age, group = column_label))+
  geom_line(color = "light grey", size = 0.3)+
  geom_line(data = dplyr::filter(bin_200_median_recon_tree_plot, bin_age <= 12000),
            mapping = aes(y = tree_cover_median, x = bin_age), color = "yellow", size = 1)+
  geom_line(data = dplyr::filter(boot_bin_200_median_recon_tree_5_95, bin_age <= 12000),
            mapping = aes(y = tree_cover_5_95, x = bin_age, group = lower_upper), color = "red",size = 0.2)+
  ylim(25,75)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  scale_x_reverse(breaks = seq(12000, 0, -2000))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))
ggsave("figs/supplement/SF6.plot_boot_binned_median.png",height = 6, width = 8, unit = "cm")

# SF7 - Gridded maps
#See "figs/tree_cover_rasters/"



# Supplementary Table 5: Data by biogeographical region
data_binned_bioregion <- bioregion_recon_tree_binned_200 %>% 
  dplyr::filter(bin_age <= 12000) %>% 
  dplyr::group_by(entity_name) %>% 
  dplyr::summarise(number = n()) %>%  
  dplyr::ungroup() %>% 
  dplyr::left_join(entity_bioregions)  %>% 
  dplyr::group_by(bioregion_name)  %>% 
  dplyr::summarise(number_records = length(unique(entity_name)), number_bins = sum(number))
rio::export(data_binned_bioregion,"figs/supplement/ST5.data_binned_bioregion.csv")  


# Supplementary Figure 7: Median tree cover based on different the binwidths used in other vegetation reconstructions
ggplot2::ggplot(data = dplyr::filter(bin_200_median_recon_tree_plot, bin_age <= 12000),
                                                      mapping = aes(y = tree_cover_median, x = bin_age, color = as.factor(column_label)))+
  geom_line()+
  geom_line(data = filter(comp_serge_median_all_locations_plot, source == "recon" & bin_age <= 12000), mapping = aes(y = tree_cover, x = bin_age, color = as.factor(column_label))) +
  geom_line(data = filter(comp_zanon_median_all_locations_plot, source == "recon" & bin_age <= 12000), mapping = aes(y = tree_cover, x = bin_age, color = as.factor(column_label)))  +
  ylim(25,75)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  scale_x_reverse(breaks = seq(12000, 0, -2000))+
  ggthemes::scale_color_colorblind(name = "Bin type", labels = c("200-years", "Serge", "Zanon"))+
  # scale_color_discrete(name = "Bin type", labels = c("200-years", "Serge", "Zanon"))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10),legend.title=element_text(size=9), legend.text=element_text(size=8))+
  theme(text = element_text(family = "Times New Roman"))
ggsave("figs/supplement/SF7.otherreconbin_median.png",height = 6, width = 10, unit = "cm")


# Supplementary Figure 8: Median tree cover calculated using different fossil sample age model estimates

lq_bin_200_median_recon_tree_plot <- lq_bin_200_median_recon_tree%>% 
  dplyr::mutate(column_label = 2001)
uq_bin_200_median_recon_tree_plot <- uq_bin_200_median_recon_tree%>% 
  dplyr::mutate(column_label = 3001)

ggplot2::ggplot(data = dplyr::filter(bin_200_median_recon_tree_plot, bin_age <= 12000),
                                                     mapping = aes(y = tree_cover_median, x = bin_age, color = as.factor(column_label)))+
  geom_line()+
  geom_line(data = filter(lq_bin_200_median_recon_tree_plot, bin_age <= 12000), mapping = aes(y = tree_cover_median, x = bin_age, color = as.factor(column_label))) +
  geom_line(data = filter(uq_bin_200_median_recon_tree_plot, bin_age <= 12000), mapping = aes(y = tree_cover_median, x = bin_age, color = as.factor(column_label)))  +
  ylim(25,75)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  scale_x_reverse(breaks = seq(12000, 0, -2000))+
  ggthemes::scale_color_colorblind(name = "Age model values", labels = c("Median", "Lower quartile", "Upper quartile"))+
  # scale_color_discrete(name = "Age model values", labels = c("Median", "Lower quartile", "Upper quartile"))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10),legend.title=element_text(size=9), legend.text=element_text(size=8))+
  theme(text = element_text(family = "Times New Roman"))
ggsave("figs/supplement/SF8.intqage_median.png",height = 6, width = 10, unit = "cm")



# Supplementary Figure 9: Median tree cover reconstructions, with reconstructions for Serge and Zanon data based ongridcell averages rather than record values
ggplot2::ggplot(data = dplyr::filter(comp_all_cell_median_loc, half_width == "None" & source != "githumbi"), 
                mapping = aes(y = average, x = bin_age, color = source))+
  geom_line()+
  geom_line(data = dplyr::filter(comp_all_cell_median_loc, half_width == "1000-year" & source != "githumbi"), 
            mapping = aes(y = average, x = bin_age, color = source), size =1)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  ggthemes::scale_color_colorblind(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  # scale_color_discrete(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10), legend.title=element_text(size=9), legend.text=element_text(size=8))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_x_reverse(limits = c(12000, 0), breaks = seq(12000, 0, -2000))
ggsave("figs/supplement/SF9.comp_all_cell_median_loc.png", height = 6, width = 10, unit = "cm")



# Supplementary Figure 10: Median tree cover based on records within the Atlantic, Boreal, Continental and Mediterranean biogeogrpahical regions only 
ggplot2::ggplot(data = dplyr::filter(comp_all_bioregiongr4_median_EU_loc, half_width == "None" & source != "githumbi"), 
                                                   mapping = aes(y = average, x = bin_age, color = source))+
  geom_line()+
  geom_line(data = dplyr::filter(comp_all_bioregiongr4_median_EU_loc, half_width == "1000-year" & source != "githumbi"), 
            mapping = aes(y = average, x = bin_age, color = source), size =1)+
  theme_bw()+
  ylab("Median tree cover (%)")+
  xlab("Bin age (cal. BP)")+
  ggthemes::scale_color_colorblind(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  # scale_color_discrete(name = "Source", labels = c("This study", "Serge", "Zanon"))+
  theme(axis.text=element_text(size=10, ),axis.title=element_text(size=10), legend.title=element_text(size=9), legend.text=element_text(size=8))+
  theme(text = element_text(family = "Times New Roman"))+
  scale_x_reverse(breaks = seq(12000, 0, -2000))
ggsave("figs/supplement/SF10.comp_all_bioregiongr4_median_EU_loc.png", height = 6, width = 10, unit = "cm")




