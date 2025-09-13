library(tidyverse)
library(here)
library(terra)
library(sp)
library(raster)
library(mapview)
library(gstat)
library(sf)
library(shapefiles)
library(ggspatial)
library(ggrepel)


# OBSERVED RICHNESS -------------------------------------------------------

# prepare data for naive occupancy
glgs <- read_tsv(here("data", "OTUtable_12S_toSP_GLG23_20240528.txt"))

glgs <- glgs |> 
  mutate(
    riverOrigin = case_when(
      is.na(DistFromPAedge) ~ "outsideGLGS",
      TRUE ~ "insideGLGS"
    )
  ) |> 
  relocate(riverOrigin, .before = DistFromPAedge)

binarise1 <- function(x, na.rm = FALSE) vegan::decostand(x, method = "pa") 


# SPECIES TO REMOVE

# RETAIN ONLY SPECIES THAT WERE DETECTED AT >1 SITE
spptoremove <- glgs %>% 
  as.data.frame() %>%
  dplyr::select(!contains(c("Spikein"))) %>% 
  dplyr::select(!contains(c("Cololabis_saira", "Sardina_pilchardus",
                            "Sardinella_lemuru", "Engraulis_ringens",
                            "Micromesistius_poutassou", "Scophthalmus_maximus",
                            "Scomberomorus_niphonius", "Trichiurus_haumela","Rana_catesbeiana"))) %>%
  #dplyr::select(!contains(c("Bos_gaurus","Bos_grunniens","Bos_taurus","Bubalus_bubalis","Capra_hircus","Ovis_aries","Cervus_nippon","Sus_scrofa","Canis_lupus","Felis_catus","Oryctolagus_cuniculu","Equus_asinus","Equus_caballus","Anas_platyrhynchos","Gallus_gallus","Melopsittacus_undulatus","Pelodiscus_sinensis","Trachemys_scripta"))) %>%   
  dplyr::select(SampleID, PCR, site, contains("OTU")) %>% 
  dplyr::rename_with(.cols = starts_with(c("Amphibia", "Aves", "Mammalia", 
                                           "Reptilia", "Teleostei", 
                                           "Actinopteri")), 
                     function(x){paste0("OTU_", x)}) |>
  group_by(site) %>% 
  summarise(across(contains("OTU"), function(x){
    sum(x > 0)
  })) %>% 
  arrange(site) %>% 
  dplyr::select(contains("OTU")) %>% 
  dplyr::mutate(across(contains("OTU"), binarise1)) |>
  dplyr::filter(!row_number() %in% c(102)) %>%
  dplyr::select(where( ~ is.numeric(.x) && sum(.x) <2)) %>%
  as.matrix %>%
  colnames()


# binarise read number and sum by site and sample replicate. 
# For instance, if a species is detected in all 3 sample replicates and all 6 
# PCRs per sample replicate, the output will be 6,6,6 for that species and site.
glgs_sum1 <- glgs %>% 
  
  # remove SpikeIn, keep "Actinopteri" because it exits in GLG as invasive sp
  dplyr::select(-starts_with("Spikein"))|>
  dplyr::select(!contains(c("Cololabis_saira", "Sardina_pilchardus",
                            "Sardinella_lemuru", "Engraulis_ringens",
                            "Micromesistius_poutassou", "Scophthalmus_maximus",
                            "Scomberomorus_niphonius", "Trichiurus_haumela","Rana_catesbeiana")))|>
  # remove all the domestic species and invasive herps
  dplyr::select(!contains(c("Bos_gaurus","Bos_grunniens","Bos_taurus","Bubalus_bubalis","Capra_hircus","Ovis_aries","Cervus_nippon","Sus_scrofa","Canis_lupus","Felis_catus","Oryctolagus_cuniculu","Equus_asinus","Equus_caballus","Anas_platyrhynchos","Gallus_gallus","Melopsittacus_undulatus","Pelodiscus_sinensis","Trachemys_scripta"))) %>%   
  # add "OTU_" to the beginning of all OTU columns
  rename_with(.cols = starts_with(c("Amphibia", "Aves", "Mammalia", 
                                    "Reptilia", "Teleostei", "Actinopteri")), 
              function(x){paste0("OTU_", x)}) |>  
  #rename_with(.cols = starts_with(c("Teleostei")), 
  #            function(x){paste0("OTU_", x)}) |> 
  
  # change numbers in OTU columns to 0/1
  mutate(across(starts_with("OTU_"), binarise1)) |> 
  dplyr::select(-SampleID, -PCR, -Rep_of_sample) |>
  
  # sum detections by site and Rep_of_sample (stage 1 replicates)
  summarise(
    across(Filter_vol:DistFromPAedge, \(x) first(x)),
    across(starts_with("OTU_"), \(x) sum(x)), 
    .by = c(site)
  ) |>
  filter(status == "sample")

# prepare data for naive occupancy
glgs_sum2 <- glgs_sum1 |>
  dplyr::select(-one_of(spptoremove)) %>%
  as.data.frame()

# prepare data for naive occupancy
glgs_sum2.1 <- glgs_sum2 |>
  mutate(across(starts_with("OTU_"), binarise1)) |>
  as.data.frame() |>
  mutate(obs_rich = rowSums(across(starts_with("OTU_"))))

glgs_sum2.1$obs_rich
str(glgs_sum2.1)

# calculate richness by taxa
# create new variable "phylum" for plotting purposes
#spp_names <- colnames(glgs_sum2[,22:337])
spp_names <- colnames(glgs_sum2[,22:322])
str(spp_names)
spp_list <- read_tsv(here("data","SP_withBODYSIZE.txt"))
spp_list_native <- spp_list %>% 
  filter(class == "Teleostei") %>% filter(status == "wild")%>% 
  mutate(otuID = paste0("OTU_", otuID))
spp_list_invasive <- spp_list %>% 
  filter(class == "Teleostei") %>% filter(status == "invasive") %>% 
  mutate(otuID = paste0("OTU_", otuID))

classesNames <- sapply(spp_names, function(x){
  strsplit(x, split = "_")[[1]][2]
})

amph <- which(classesNames == "Amphibia" | classesNames == "Reptilia")
str(amph)
glgs_sum2.1$amph_rich <- rowSums(glgs_sum2.1[,21+amph]) # at 1/18 cutoff
glgs_sum2.1$amph_rich3 <- rowSums(glgs_sum2[,21+amph] > 2, na.rm = TRUE) # at 3/18 cutoff

mammal <- which(classesNames == "Mammalia")
str(mammal)
glgs_sum2.1$mammal_rich <- rowSums(glgs_sum2.1[,21+mammal]) # at 1/18 cutoff
glgs_sum2.1$mammal_rich3 <- rowSums(glgs_sum2[,21+mammal] > 2, na.rm = TRUE) # at 3/18 cutoff

ave <- which(classesNames == "Aves")
str(mammal)
glgs_sum2.1$ave_rich <- rowSums(glgs_sum2.1[,21+ave]) # at 1/18 cutoff
glgs_sum2.1$ave_rich3 <- rowSums(glgs_sum2[,21+ave] > 2, na.rm = TRUE) # at 3/18 cutoff

terr <- which(classesNames %in% c("Mammalia", "Amphibia", "Reptilia", "Aves"))
glgs_sum2.1$terr_rich <- rowSums(glgs_sum2[,21+terr] > 0, na.rm = TRUE) # at 1/18 cutoff
glgs_sum2.1$terr_rich3 <- rowSums(glgs_sum2[,21+terr] > 2, na.rm = TRUE) # at 3/18 cutoff

fish <- which(classesNames == "Teleostei" | classesNames == "Actinopteri")
str(fish)
glgs_sum2.1$fish_rich <- rowSums(glgs_sum2[,21+fish] > 0, na.rm = TRUE) # at 1/18 cutoff
glgs_sum2.1$fish_rich3 <- rowSums(glgs_sum2[,21+fish] > 2, na.rm = TRUE) # at 3/18 cutoff

fish_native <- which(spp_names %in% spp_list_native$otuID)
str(fish_native)
glgs_sum2.1$fish_native_rich <- rowSums(glgs_sum2.1[,21+fish_native]) # at 1/18 cutoff
glgs_sum2.1$fish_native_rich3 <- rowSums(glgs_sum2[,21+fish_native] > 2, na.rm = TRUE) # at 3/18 cutoff

fish_nonnative <- which(spp_names %in% spp_list_invasive$otuID)
str(fish_nonnative)
glgs_sum2.1$fish_nonnative_rich <- rowSums(glgs_sum2.1[,21+fish_nonnative]) # at 1/18 cutoff
glgs_sum2.1$fish_nonnative_rich3 <- rowSums(glgs_sum2[,21+fish_nonnative] > 2, na.rm = TRUE) # at 3/18 cutoff

# change names for plotting
oldnames <- spp_names

# extract species names from OTUs
aa <- str_split(oldnames, "_")
spp <- rep(NA, length(aa))

for(i in 1:length(aa)) {
  tmp <- aa[[i]][5:6]
  spp[i] <- paste(tmp[1],".",tmp[2], sep = "")
}

spp

glgs_sum2 <-  glgs_sum2 |>
  rename_with(~ spp, all_of(oldnames))
names(glgs_sum2)

glgs_sum2.1 <-  glgs_sum2.1 |>
  rename_with(~ spp, all_of(oldnames))
names(glgs_sum2.1)


#######################################################################################################
# Plotting the maps of Species richness observations and predictions with the same point size ranges

# 1. SETUP: Prepare spatial objects and define constants

# Load shared spatial data
GLGS <- st_read(here("data", "glgs", "GLG_Region_FULL_UTM.shp"))
PA <- st_read(here("data", "glgs", "GLGS_PAs.shp"))
raster_map <- raster(here("data", "glgs", "forest_1km.tif"))
raster_df <- as.data.frame(raster_map, xy = TRUE)
GLGS <- st_transform(GLGS, 4326)
PA <- st_transform(PA, 4326)

# Convert the main observation data to an sf object
glgs_obs_data_sf <- st_as_sf(glgs_sum2.1, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# List of domestic species to remove from model outputs
domestic_spp_model <- c("Bos.gaurus","Bos.grunniens","Bos.taurus","Bubalus.bubalis",
                        "Capra.hircus","Ovis.aries","Cervus.nippon","Sus.scrofa",
                        "Canis.lupus","Felis.catus","Oryctolagus.cuniculu",
                        "Equus.asinus","Equus.caballus","Anas.platyrhynchos",
                        "Gallus.gallus","Melopsittacus.undulatus","Pelodiscus.sinensis",
                        "Trachemys.scripta","Rana.catesbeiana")

# Helper lists for native/invasive fish models
spp_list_native <- read_tsv(here("data","SP_withBODYSIZE.txt")) %>% 
  filter(class == "Teleostei", status == "wild") %>% 
  mutate(spp = gsub(" ", ".", species), spp = paste0("Mean.", spp))

spp_list_invasive <- read_tsv(here("data","SP_withBODYSIZE.txt")) %>% 
  filter(class == "Teleostei", status == "invasive") %>% 
  mutate(spp = gsub(" ", ".", species), spp = paste0("Mean.", spp))

# 2. FUNCTIONS: Create generic functions for data processing and plotting

#' Process Model Output Data (with native/invasive fish logic)
process_model_data <- function(msom_path, occplus_path, join_data, domestic_spp_to_remove = NULL, is_fish_data = FALSE) {
  
  msom_df <- read.csv(here(msom_path))
  occplus_df <- read.csv(here(occplus_path))
  
  # --- Process MSOM data ---
  if (!is.null(domestic_spp_to_remove)) {
    msom_df <- msom_df %>% dplyr::select(!contains(domestic_spp_to_remove))
  }
  
  # Always calculate total summed probabilities
  msom_df <- msom_df %>% 
    mutate(summed_probs = rowSums(across(starts_with("Mean."))))
  
  if (is_fish_data) {
    col_names_msom <- colnames(msom_df)
    fish_native_msom <- which(col_names_msom %in% spp_list_native$spp)
    msom_df$fish_native_rich <- rowSums(msom_df[fish_native_msom])
    
    fish_nonnative_msom <- which(col_names_msom %in% spp_list_invasive$spp)
    msom_df$fish_nonnative_rich <- rowSums(msom_df[fish_nonnative_msom])
  }
  
  msom_sf <- msom_df %>% st_as_sf(coords = 6:7, crs = 4326)
  
  # --- Process OccPlus data ---
  occplus_df <- inner_join(join_data, occplus_df, by = "site")
  if (!is.null(domestic_spp_to_remove)) {
    occplus_df <- occplus_df %>% dplyr::select(!contains(domestic_spp_to_remove))
  }
  
  # Always calculate total summed probabilities
  occplus_df <- occplus_df %>% 
    mutate(summed_probs = rowSums(across(starts_with("Mean."))))
  
  if (is_fish_data) {
    col_names_occ <- colnames(occplus_df)
    fish_native_occ <- which(col_names_occ %in% spp_list_native$spp)
    occplus_df$fish_native_rich <- rowSums(occplus_df[fish_native_occ])
    
    fish_nonnative_occ <- which(col_names_occ %in% spp_list_invasive$spp)
    occplus_df$fish_nonnative_rich <- rowSums(occplus_df[fish_nonnative_occ])
  }
  
  occplus_sf <- occplus_df %>% st_as_sf(coords = 3:4, crs = 4326)
  
  return(list(msom_sf = msom_sf, occplus_sf = occplus_sf))
}


#' Create and Save a Richness Map
create_richness_map <- function(data_sf, size_var, title, point_fill, size_labels, filename) {
  
  breaks <- as.numeric(gsub(".*-([0-9]+)|≥([0-9]+)|([0-9]+)", "\\1\\2\\3", size_labels))
  breaks <- breaks[!is.na(breaks)]
  
  size_mapping_func <- function(value) {
    if (value <= breaks[1]) return(1)
    if (value <= breaks[2]) return(2)
    if (value <= breaks[3]) return(3)
    if (value <= breaks[4]) return(4)
    return(5)
  }
  
  data_sf$size_mapped <- sapply(data_sf[[size_var]], size_mapping_func)
  
  p <- ggplot() +
    geom_raster(data = raster_df, aes(x = x, y = y, fill = forest_1km)) +
    scale_fill_gradient(low = "white", high = "grey70", na.value = NA) +
    guides(fill = "none") +
    geom_sf(data = GLGS, color = "darkblue", fill = NA, size = 1.2) +
    geom_sf(data = PA, color = "green", fill = "green", alpha = 0.8) +
    geom_sf(data = data_sf,
            aes(size = size_mapped, alpha = riverOrigin),
            shape = 21, color = "black", fill = point_fill,
            show.legend = 'point', inherit.aes = FALSE) +
    labs(size = title) + # Use labs() for the title to respect newlines
    theme(
      legend.position = c(0.35, 0.65),
      legend.title = element_text(face = "bold", hjust = 0, size = 12, lineheight = 1.1),
      legend.spacing = unit(0, 'cm'),
      legend.background = element_rect(fill = NA, color = NA),
      axis.title = element_blank()
    ) +
    scale_alpha_manual(name = NULL, values = c("insideGLGS" = 1, "outsideGLGS" = 0.5), labels = c("Inside PAs", "Outside PAs")) +
    guides(alpha = guide_legend(override.aes = list(size = 4), keywidth = 1, keyheight = 1),
           size = guide_legend(order = 1)) +
    coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
    scale_size_continuous( # The name argument here is overridden by labs(), but we still need the scale for labels
      name = title, 
      limits = c(1, 5), 
      breaks = 1:5, 
      labels = size_labels
    ) +
    annotation_scale(location = "bl", width_hint = 0.2, style = "ticks") +
    annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
    geom_tile(data = data.frame(x = 97.1, y = 25.7), aes(x, y), fill = "green", width = 0.22, height = 0.14, alpha = 0.8) +
    annotate("text", x = 97.25, y = 25.7, label = "Protected Areas", hjust = 0, size = 4)
  
  ggsave(filename = here("outputs", "maps1018", filename), plot = p, width = 4, height = 8, dpi = 300)
  
  return(invisible(p)) # Return plot without printing it
}

# 3. ANALYSIS: Generate maps for each taxonomic group

# --- Terrestrial Vertebrates (TERR) ---
terr_params <- list(
  color = "brown",
  labels = c("1-17", "18-35", "36-53", "54-71", "72-86"),
  msom_path = "outputs/msom_outputs_Viorel covs/AmphRept_outputs_20241018/occupancy_bysite_msom_AmphRept_20241018.csv", # Placeholder, needs combined data
  occplus_path = "outputs/occPlus_outputs_Viorel covs/AmphRept_outputs_20241018_3factors/occupancybysite_occPlus_AmphRept_20241018_3factors.csv" # Placeholder
  # Note: The combined terrestrial model data requires a separate, more complex processing step not included in the simple function.
  # For now, we only plot the observed data as the logic is clear.
)
create_richness_map(glgs_obs_data_sf, "terr_rich", "Observed terrestrial\nspecies richness\nat ≥1/18 detections (PCRs)", terr_params$color, terr_params$labels, "GLGS_naive_TERR_richness_onlywild.jpg")
create_richness_map(glgs_obs_data_sf, "terr_rich3", "Observed terrestrial\nspecies richness\nat ≥3/18 detections (PCRs)", terr_params$color, terr_params$labels, "GLGS_naive_TERR_richness_onlywild_cutoff3of18.jpg")

# Modelled terrestrial maps would require a dedicated data prep step first.
# Load the individual model outputs for each terrestrial group
msom_amphrept <- read.csv(here("outputs/msom_outputs_Viorel covs/AmphRept_outputs_20241018/occupancy_bysite_msom_AmphRept_20241018.csv"))
msom_mammal <- read.csv(here("outputs/msom_outputs_Viorel covs/Mammal_outputs_20241018/occupancy_bysite_msom_Mammal_20241018.csv"))
msom_ave <- read.csv(here("outputs/msom_outputs_Viorel covs/Ave_outputs_20241018/occupancy_bysite_msom_Ave_20241018.csv"))
occplus_amphrept <- read.csv(here("outputs/occPlus_outputs_Viorel covs/AmphRept_outputs_20241018_3factors/occupancybysite_occPlus_AmphRept_20241018_3factors.csv"))
occplus_mammal <- read.csv(here("outputs/occPlus_outputs_Viorel covs/Mammal_outputs_20241018_4factors/occupancybysite_occPlus_Mammal_20241018_4factors.csv"))
occplus_ave <- read.csv(here("outputs/occPlus_outputs_Viorel covs/Ave_outputs_20241018_2factors/occupancybysite_occPlus_Ave_20241018_2factors.csv"))

# Process and combine MSOM data
msom_amphrept_sum <- msom_amphrept %>% 
  mutate(sum_amphrept = rowSums(across(starts_with("Mean.")))) %>%
  dplyr::select(site, longitude, latitude, sum_amphrept)

msom_mammal_sum <- msom_mammal %>% 
  dplyr::select(!contains(domestic_spp_model)) %>%
  mutate(sum_mammal = rowSums(across(starts_with("Mean.")))) %>%
  dplyr::select(site, sum_mammal)

msom_ave_sum <- msom_ave %>%
  dplyr::select(!contains(domestic_spp_model)) %>%
  mutate(sum_ave = rowSums(across(starts_with("Mean.")))) %>%
  dplyr::select(site, sum_ave)

# Join the summed probabilities and calculate total terrestrial richness
msom_terr_sf <- left_join(msom_amphrept_sum, msom_mammal_sum, by = "site") %>%
  left_join(msom_ave_sum, by = "site") %>%
  mutate(summed_probs = sum_amphrept + sum_mammal + sum_ave) %>%
  # Also need to get 'riverOrigin' for plotting alpha
  left_join(dplyr::select(glgs_sum2.1, site, riverOrigin), by = "site") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4236)

# Process and combine OccPlus data
occplus_amphrept_sum <- occplus_amphrept %>% 
  mutate(sum_amphrept = rowSums(across(starts_with("Mean.")))) %>%
  dplyr::select(site, sum_amphrept)

occplus_mammal_sum <- occplus_mammal %>%
  dplyr::select(!contains(domestic_spp_model)) %>%
  mutate(sum_mammal = rowSums(across(starts_with("Mean.")))) %>%
  dplyr::select(site, sum_mammal)

occplus_ave_sum <- occplus_ave %>%
  dplyr::select(!contains(domestic_spp_model)) %>%
  mutate(sum_ave = rowSums(across(starts_with("Mean.")))) %>%
  dplyr::select(site, sum_ave)

# Join, calculate total, and then join with coordinate/riverOrigin info from glgs_sum2.1
occplus_terr_sf <- left_join(occplus_amphrept_sum, occplus_mammal_sum, by = "site") %>%
  left_join(occplus_ave_sum, by = "site") %>%
  mutate(summed_probs = sum_amphrept + sum_mammal + sum_ave) %>%
  inner_join(glgs_sum2.1, by = "site") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Now, call the plotting function with the newly created data frames
create_richness_map(msom_terr_sf, "summed_probs", "Predicted terrestrial\nspecies richness (MSOM)\n[summed probabilities]", terr_params$color, terr_params$labels, "GLGS_msom_TERR_richness1018_nodomestic.jpg")
create_richness_map(occplus_terr_sf, "summed_probs", "Predicted terrestrial\nspecies richness (OccPlus)\n[summed probabilities]", terr_params$color, terr_params$labels, "GLGS_occPlus_TERR_richness1018_nodomestic.jpg")


# --- Fish (Total, Native, Invasive) ---
fish_params <- list(
  color = "blue",
  labels = c("0-10", "11-20", "21-30", "31-40", "41-48"),
  msom_path = "outputs/msom_outputs_Viorel covs/FISH_outputs_20241018/occupancy_bysite_msom_FISH_20241018.csv",
  occplus_path = "outputs/occPlus_outputs_Viorel covs/FISH_outputs_20241018_5factors/occupancybysite_occPlus_FISH_20241018_5factors.csv"
)
native_fish_params <- list(color = "royalblue", labels = c("0-4", "5-9", "10-14", "15-19", "≥20"))
invasive_fish_params <- list(color = "red", labels = c("0-4", "5-9", "10-14", "15-19", "≥20"))

fish_model_data <- process_model_data(fish_params$msom_path, fish_params$occplus_path, glgs_sum2.1, is_fish_data = TRUE)

create_richness_map(glgs_obs_data_sf, "fish_rich", "Observed fish\nspecies richness\nat ≥1/18 detections (PCRs)", fish_params$color, fish_params$labels, "GLGS_naive_FISH_richness.jpg")
create_richness_map(glgs_obs_data_sf, "fish_rich3", "Observed fish\nspecies richness\nat ≥3/18 detections (PCRs)", fish_params$color, fish_params$labels, "GLGS_naive_FISH_richness_cutoff3of18.jpg")
create_richness_map(fish_model_data$msom_sf, "summed_probs", "Predicted fish\nspecies richness (MSOM)\n[summed probabilities]", fish_params$color, fish_params$labels, "GLGS_msom_FISH_richness1018.jpg")
create_richness_map(fish_model_data$occplus_sf, "summed_probs", "Predicted fish\nspecies richness (OccPlus)\n[summed probabilities]", fish_params$color, fish_params$labels, "GLGS_occPlus_FISH_richness1018_5factors.jpg")
create_richness_map(glgs_obs_data_sf, "fish_native_rich", "Observed native fish\nspecies richness\nat ≥1/18 detections (PCRs)", native_fish_params$color, native_fish_params$labels, "GLGS_naive_nativeFISH_richness.jpg")
create_richness_map(glgs_obs_data_sf, "fish_native_rich3", "Observed native fish\nspecies richness\nat ≥3/18 detections (PCRs)", native_fish_params$color, native_fish_params$labels, "GLGS_naive_nativeFISH_richness_cutoff3of18.jpg")
create_richness_map(fish_model_data$msom_sf, "fish_native_rich", "Predicted Native Fish\nspecies richness (MSOM)\n[summed probabilities]", native_fish_params$color, native_fish_params$labels, "GLGS_msom_nativeFISH_richness1018.jpg")
create_richness_map(fish_model_data$occplus_sf, "fish_native_rich", "Predicted Native Fish\nspecies richness (OccPlus)\n[summed probabilities]", native_fish_params$color, native_fish_params$labels, "GLGS_occPlus_nativeFISH_richness1018.jpg")
create_richness_map(glgs_obs_data_sf, "fish_nonnative_rich", "Observed non-native fish\nspecies richness\nat ≥1/18 detections (PCRs)", invasive_fish_params$color, invasive_fish_params$labels, "GLGS_naive_invasiveFISH_richness.jpg")
create_richness_map(glgs_obs_data_sf, "fish_nonnative_rich3", "Observed non-native fish\nspecies richness\nat ≥3/18 detections (PCRs)", invasive_fish_params$color, invasive_fish_params$labels, "GLGS_naive_invasiveFISH_richness_cutoff3of18.jpg")
create_richness_map(fish_model_data$msom_sf, "fish_nonnative_rich", "Predicted Non-native Fish\nspecies richness (MSOM)\n[summed probabilities]", invasive_fish_params$color, invasive_fish_params$labels, "GLGS_msom_invasiveFISH_richness1018.jpg")
create_richness_map(fish_model_data$occplus_sf, "fish_nonnative_rich", "Predicted Non-native Fish\nspecies richness (OccPlus)\n[summed probabilities]", invasive_fish_params$color, invasive_fish_params$labels, "GLGS_occPlus_invasiveFISH_richness1018.jpg")


# --- Mammals ---
mammal_params <- list(
  color = "palevioletred4",
  labels = c("1-7", "8-14", "15-21", "22-28", "29-36"),
  msom_path = "outputs/msom_outputs_Viorel covs/Mammal_outputs_20241018/occupancy_bysite_msom_Mammal_20241018.csv",
  occplus_path = "outputs/occPlus_outputs_Viorel covs/Mammal_outputs_20241018_4factors/occupancybysite_occPlus_Mammal_20241018_4factors.csv"
)
mammal_model_data <- process_model_data(mammal_params$msom_path, mammal_params$occplus_path, glgs_sum2.1, domestic_spp_model)
create_richness_map(glgs_obs_data_sf, "mammal_rich", "Observed Mammalia\nspecies richness\nat ≥1/18 detections (PCRs)", mammal_params$color, mammal_params$labels, "GLGS_naive_mammal_richness_onlywild.jpg")
create_richness_map(glgs_obs_data_sf, "mammal_rich3", "Observed Mammalia\nspecies richness\nat ≥3/18 detections (PCRs)", mammal_params$color, mammal_params$labels, "GLGS_naive_mammal_richness_onlywild_cutoff3of18.jpg")
create_richness_map(mammal_model_data$msom_sf, "summed_probs", "Predicted Mammalia\nspecies richness (MSOM)\n[summed probabilities]", mammal_params$color, mammal_params$labels, "GLGS_msom_Mammal_richness1018_nodomestic.jpg")
create_richness_map(mammal_model_data$occplus_sf, "summed_probs", "Predicted Mammalia\nspecies richness (OccPlus)\n[summed probabilities]", mammal_params$color, mammal_params$labels, "GLGS_occPlus_Mammal_richness1018_nodomestic.jpg")


# --- Aves ---
aves_params <- list(
  color = "deepskyblue4",
  labels = c("0-6", "7-12", "13-18", "19-24", "25-31"),
  msom_path = "outputs/msom_outputs_Viorel covs/Ave_outputs_20241018/occupancy_bysite_msom_Ave_20241018.csv",
  occplus_path = "outputs/occPlus_outputs_Viorel covs/Ave_outputs_20241018_2factors/occupancybysite_occPlus_Ave_20241018_2factors.csv"
)
aves_model_data <- process_model_data(aves_params$msom_path, aves_params$occplus_path, glgs_sum2.1, domestic_spp_model)
create_richness_map(glgs_obs_data_sf, "ave_rich", "Observed Aves\nspecies richness\nat ≥1/18 detections (PCRs)", aves_params$color, aves_params$labels, "GLGS_naive_ave_richness_onlywild.jpg")
create_richness_map(glgs_obs_data_sf, "ave_rich3", "Observed Aves\nspecies richness\nat ≥3/18 detections (PCRs)", aves_params$color, aves_params$labels, "GLGS_naive_ave_richness_onlywild_cutoff3of18.jpg")
create_richness_map(aves_model_data$msom_sf, "summed_probs", "Predicted Aves\nspecies richness (MSOM)\n[summed probabilities]", aves_params$color, aves_params$labels, "GLGS_msom_Ave_richness1018_nodomestic.jpg")
create_richness_map(aves_model_data$occplus_sf, "summed_probs", "Predicted Aves\nspecies richness (OccPlus)\n[summed probabilities]", aves_params$color, aves_params$labels, "GLGS_occPlus_Ave_richness1018_nodomestic.jpg")


# --- Herptiles (Amphibians & Reptiles) ---
herp_params <- list(
  color = "blueviolet",
  labels = c("0-4", "5-8", "9-12", "13-16", "17-21"),
  msom_path = "outputs/msom_outputs_Viorel covs/AmphRept_outputs_20241018/occupancy_bysite_msom_AmphRept_20241018.csv",
  occplus_path = "outputs/occPlus_outputs_Viorel covs/AmphRept_outputs_20241018_3factors/occupancybysite_occPlus_AmphRept_20241018_3factors.csv"
)
herp_model_data <- process_model_data(herp_params$msom_path, herp_params$occplus_path, glgs_sum2.1)
create_richness_map(glgs_obs_data_sf, "amph_rich", "Observed herptile\nspecies richness\nat ≥1/18 detections (PCRs)", herp_params$color, herp_params$labels, "GLGS_naive_amphrept_richness_onlywild.jpg")
create_richness_map(glgs_obs_data_sf, "amph_rich3", "Observed herptile\nspecies richness\nat ≥3/18 detections (PCRs)", herp_params$color, herp_params$labels, "GLGS_naive_amphrept_richness_onlywild_cutoff3of18.jpg")
create_richness_map(herp_model_data$msom_sf, "summed_probs", "Predicted herptile\nspecies richness (MSOM)\n[summed probabilities]", herp_params$color, herp_params$labels, "GLGS_msom_AmphRept_richness1018.jpg")
create_richness_map(herp_model_data$occplus_sf, "summed_probs", "Predicted herptile\nspecies richness (OccPlus)\n[summed probabilities]", herp_params$color, herp_params$labels, "GLGS_occPlus_AmphRept_richness1018.jpg")



#######################################################################################################
# Plotting the maps of Species richness predictions for each order of mammals

# --- Data Preparation for Mammal Orders ---

# Load model outputs once
msom_mammal_raw <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Mammal_outputs_20241018", "occupancy_bysite_msom_Mammal_20241018.csv"))
occplus_mammal_raw <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Mammal_outputs_20241018_4factors", "occupancybysite_occPlus_Mammal_20241018_4factors.csv"))

# Clean domestic species
msom_mammal_clean <- msom_mammal_raw %>% 
  dplyr::select(!contains(domestic_spp_model))
occplus_mammal_clean <- inner_join(glgs_sum2.1, occplus_mammal_raw, by = "site") %>% 
  dplyr::select(!contains(domestic_spp_model))

# Extract Order information from the species list
# (Assuming 'spp_list' is the one from your previous script)
spp_list_mammal_orders <- spp_list %>%
  filter(class == "Mammalia") %>%
  mutate(ORDER = str_extract(otuID, "(?<=_)[^_]+(?=_)")) %>%
  mutate(spp_col_name = paste0("Mean.", gsub(" ", ".", species)))

# Get a unique list of orders to loop through
mammal_orders <- unique(spp_list_mammal_orders$ORDER)

# Loop to calculate richness for each order and store in the dataframes
for (order_name in mammal_orders) {
  
  # Get species column names for the current order
  spp_in_order <- spp_list_mammal_orders %>%
    filter(ORDER == order_name) %>%
    pull(spp_col_name)
  
  # Define the new richness column name
  richness_col <- paste0(order_name, "_rich")
  
  # Calculate richness and add as a new column
  # Use `any_of()` to avoid errors if a species is not in the model output
  msom_mammal_clean[[richness_col]] <- rowSums(msom_mammal_clean %>% dplyr::select(any_of(spp_in_order)))
  occplus_mammal_clean[[richness_col]] <- rowSums(occplus_mammal_clean %>% dplyr::select(any_of(spp_in_order)))
}

# Convert the final dataframes to sf objects
msom_mammal_orders_sf <- st_as_sf(msom_mammal_clean, coords = 6:7, crs = 4326)
occplus_mammal_orders_sf <- st_as_sf(occplus_mammal_clean, coords = c("longitude", "latitude"), crs = 4326)


# --- Define Plotting Parameters and Generate Maps ---

# Create a list that holds all the specific parameters for each order's map
order_plotting_params <- list(
  Primates = list(
    labels = c("0-0.35", "0.36-0.7", "0.71-1.05", "1.06-1.4", "1.41-1.75"),
    rich_col = "Primates_rich"
  ),
  Rodentia = list(
    labels = c("0-3.3", "3.4-6.6", "6.7-9.9", "10-13.2", "13.3-16.52"),
    rich_col = "Rodentia_rich"
  ),
  Chiroptera = list(
    labels = c("0-1.62", "1.63-3.24", "3.25-4.86", "4.87-6.48", "6.49-8.11"),
    rich_col = "Chiroptera_rich"
  ),
  Carnivora = list(
    labels = c("0-0.73", "0.74-1.46", "1.47-2.19", "2.2-2.92", "2.93-3.65"),
    rich_col = "Carnivora_rich"
  ),
  Artiodactyla = list(
    labels = c("0-0.81", "0.82-1.62", "1.63-2.43", "2.44-3.24", "3.25-4.06"),
    rich_col = "Artiodactyla_rich"
  ),
  Eulipotyphla = list(
    labels = c("0.23-1.18", "1.19-2.13", "2.14-3.08", "3.09-4.03", "4.04-4.98"),
    rich_col = "Eulipotyphla_rich"
  )
)

# Now, loop through the parameters and create the maps for each order
for (order_name in names(order_plotting_params)) {
  
  params <- order_plotting_params[[order_name]]
  
  # --- Create MSOM Map ---
  msom_title <- paste0("Predicted ", order_name, "\nspecies richness (MSOM)\n[summed probabilities]")
  msom_filename <- paste0("GLGS_msom_", order_name, "_richness1018.jpg")
  
  create_richness_map(
    data_sf = msom_mammal_orders_sf,
    size_var = params$rich_col,
    title = msom_title,
    point_fill = "palevioletred4", # Shared color for all mammal orders
    size_labels = params$labels,
    filename = msom_filename
  )
  
  # --- Create OccPlus Map ---
  occplus_title <- paste0("Predicted ", order_name, "\nspecies richness (OccPlus)\n[summed probabilities]")
  occplus_filename <- paste0("GLGS_occPlus_", order_name, "_richness1018.jpg")
  
  create_richness_map(
    data_sf = occplus_mammal_orders_sf,
    size_var = params$rich_col,
    title = occplus_title,
    point_fill = "palevioletred4",
    size_labels = params$labels,
    filename = occplus_filename
  )
}



#######################################################################################################

# MSOM predictions --------------------------------------------------------

msom_occ <- read.csv(here("outputs", "msom_outputs_Viorel covs", "AmphRept_outputs_20241018",
                          "occupancy_bysite_msom_AmphRept_20241018.csv"))
msom_occ <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Mammal_outputs_20241018",
                          "occupancy_bysite_msom_Mammal_20241018.csv"))
msom_occ <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Ave_outputs_20241018",
                          "occupancy_bysite_msom_Ave_20241018.csv"))
msom_occ <- read.csv(here("outputs", "msom_outputs_Viorel covs", "FISH_outputs_20241018",
                          "occupancy_bysite_msom_FISH_20241018.csv"))
str(msom_occ)
colnames(msom_occ)

# sum mean occupancies for all species for each site
# this provides a relative measure of species richness
msom_occ <- msom_occ |>
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))


# occPlus predictions -----------------------------------------------------

occPlus_occ <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors",
                               "occupancybysite_occPlus_AmphRept_20241018_3factors.csv"))

occPlus_occ <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Mammal_outputs_20241018_4factors",
                             "occupancybysite_occPlus_Mammal_20241018_4factors.csv"))

occPlus_occ <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Ave_outputs_20241018_2factors",
                               "occupancybysite_occPlus_Ave_20241018_2factors.csv"))

occPlus_occ <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors",
                             "occupancybysite_occPlus_FISH_20241018_5factors.csv"))

str(occPlus_occ)
colnames(occPlus_occ)

# sum mean occupancies for all species for each site
# this provides a relative measure of species richness
occPlus_occ <- occPlus_occ |>
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))

occPlus_occ <- inner_join(glgs_sum2.1[,1:21], occPlus_occ, by = "site")

str(occPlus_occ)

write.csv(occPlus_occ, here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors",
                            "occupancypredictions_occPlus_AmphRept_20241018_3factors.csv"))

write.csv(occPlus_occ, here("outputs", "occPlus_outputs_Viorel covs", "Mammal_outputs_20241018_4factors",
                            "occupancypredictions_occPlus_Mammal_20241018_4factors.csv"))

write.csv(occPlus_occ, here("outputs", "occPlus_outputs_Viorel covs", "Ave_outputs_20241018_2factors",
                            "occupancypredictions_occPlus_Ave_20241018_2factors.csv"))

write.csv(occPlus_occ, here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors",
                            "occupancypredictions_occPlus_FISH_20241018_5factors.csv"))


# --- Data Preparation for Individual Species Maps ---

# Prepare the sf object for observed data plots.
trial <- st_as_sf(glgs_sum2, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Centralize the creation of species lists for each taxonomic group.
get_species_list <- function(classes_of_interest) {
  spp_oldnames_df <- data.frame(spp_name = spp_names) %>% 
    filter(str_detect(spp_name, str_c("OTU_(", str_c(classes_of_interest, collapse = "|"), ")_")))
  aa <- str_split(spp_oldnames_df$spp_name, "_")
  spp_list <- sapply(aa, function(x) paste(x[5], x[6], sep = "."))
  spp_list <- intersect(spp_list, names(trial))
  return(sort(spp_list))
}

species_lists <- list(
  AmphRept = get_species_list(c("Amphibia", "Reptilia")),
  Mammal = get_species_list("Mammalia"),
  Aves = get_species_list("Aves"),
  FISH = get_species_list(c("Teleostei", "Actinopteri"))
)

# Define a generic size-mapping function for occupancy probabilities.
prob_size_mapping <- function(value) {
  cut(value, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = 1:5, right = TRUE)
}


# --- Create Plotting Functions for Individual Species ---

# Base map function to create shared layers (background, borders, etc.).
create_base_map <- function() {
  ggplot() +
    geom_sf(data = GLGS, color = "darkblue", fill = "lightblue", alpha = 0.8, linetype = "solid", size = 1.5) +
    geom_sf(data = PA, fill = "darkgreen", alpha = 0.8) +
    annotation_scale(location = "bl", width_hint = 0.2, style = "ticks", text_cex = 0.8, line_width = 0.5, bar_cols = c("black", "white")) +
    theme(axis.title = element_blank()) +
    geom_tile(data = data.frame(x = 97.25, y = 25.8), aes(x, y), fill = "darkgreen", width = 0.22, height = 0.14, alpha = 0.8) +
    annotate("text", x = 97.4, y = 25.8, label = "Protected Areas", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.25, y = 26), aes(x, y), fill = "lightblue", width = 0.22, height = 0.14, alpha = 0.6, color = "blue") +
    annotate("text", x = 97.4, y = 26, label = "Gaoligong Region", hjust = 0, size = 3)
}

# Function to plot observed detection maps.
plot_observed_map <- function(spp_col_name) {
  spp_label <- gsub("\\.", " ", spp_col_name)
  trial_detected <- trial %>% filter(!!sym(spp_col_name) > 0)
  trial_not_detected <- trial %>% filter(!!sym(spp_col_name) == 0)
  
  p <- create_base_map() +
    geom_sf(data = trial_not_detected, shape = 21, fill = NA, color = "brown", size = 4, alpha = 0.7) +
    geom_sf(data = trial_detected, aes(fill = !!sym(spp_col_name)), shape = 21, color = "brown", size = 4, alpha = 0.7) +
    scale_fill_gradient(low = "yellow", high = "red", name = paste("Detection of\n", spp_label), limits = c(1, 18)) +
    theme(legend.position = c(0.25, 0.72), legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2), legend.background = element_rect(fill = NA, color = NA)) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 6)) +
    coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE, default_crs = sf::st_crs(4326))
  return(p)
}

# Function to plot MSOM predicted occupancy maps.
plot_msom_map <- function(spp_col_name, msom_data_sf) {
  spp_label <- gsub("\\.", " ", spp_col_name)
  msom_spp_col <- paste0("Mean.", spp_col_name)
  msom_data_sf$size_mapped <- as.numeric(prob_size_mapping(msom_data_sf[[msom_spp_col]]))
  
  p <- create_base_map() +
    geom_sf(data = msom_data_sf, aes(size = size_mapped), shape = 21, color = "black", fill = 'purple', alpha = 0.6, show.legend = 'point') +
    labs(size = paste("MSOM\n", "Occupancy probablity\n", spp_label)) +
    theme(legend.position = c(0.32, 0.68), legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2), legend.background = element_rect(fill = NA, color = NA)) +
    scale_size_continuous(limits = c(1, 5), breaks = 1:5, labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
    coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE, default_crs = sf::st_crs(4326))
  return(p)
}

# Function to plot OccPlus predicted occupancy maps.
plot_occplus_map <- function(spp_col_name, occplus_data_sf) {
  spp_label <- gsub("\\.", " ", spp_col_name)
  occplus_spp_col <- paste0("Mean.", spp_col_name)
  occplus_data_sf$size_mapped <- as.numeric(prob_size_mapping(occplus_data_sf[[occplus_spp_col]]))
  
  p <- create_base_map() +
    geom_sf(data = occplus_data_sf, aes(size = size_mapped), shape = 21, color = "darkred", fill = 'orange', alpha = 0.6, show.legend = 'point') +
    labs(size = paste("OccPlus\n", "Occupancy probablity\n", spp_label)) +
    theme(legend.position = c(0.32, 0.68), legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2), legend.background = element_rect(fill = NA, color = NA)) +
    scale_size_continuous(limits = c(1, 5), breaks = 1:5, labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
    coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE, default_crs = sf::st_crs(4326))
  return(p)
}


# --- Loop Through Each Taxonomic Group and Generate PDFs ---
library(gridExtra)

# A list to hold all the necessary parameters for each plotting task.
plotting_tasks <- list(
  AmphRept = list(pdf_path = here("outputs","maps1018","AmphRept_DETvsMSOMvsOCC_maps_for_eachSP1018_3factors.pdf"), msom_path = here("outputs", "msom_outputs_Viorel covs", "AmphRept_outputs_20241018", "occupancy_bysite_msom_AmphRept_20241018.csv"), occplus_path = here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "occupancybysite_occPlus_AmphRept_20241018_3factors.csv"), spp_list = species_lists$AmphRept),
  Mammal = list(pdf_path = here("outputs","maps1018","Mammal_DETvsMSOMvsOCC_maps_for_eachSP1018_4factors.pdf"), msom_path = here("outputs", "msom_outputs_Viorel covs", "Mammal_outputs_20241018", "occupancy_bysite_msom_Mammal_20241018.csv"), occplus_path = here("outputs", "occPlus_outputs_Viorel covs", "Mammal_outputs_20241018_4factors", "occupancybysite_occPlus_Mammal_20241018_4factors.csv"), spp_list = species_lists$Mammal),
  Aves = list(pdf_path = here("outputs","maps1018","Ave_DETvsMSOMvsOCC_maps_for_eachSP1018_2factors.pdf"), msom_path = here("outputs", "msom_outputs_Viorel covs", "Ave_outputs_20241018", "occupancy_bysite_msom_Ave_20241018.csv"), occplus_path = here("outputs", "occPlus_outputs_Viorel covs", "Ave_outputs_20241018_2factors", "occupancybysite_occPlus_Ave_20241018_2factors.csv"), spp_list = species_lists$Aves),
  FISH = list(pdf_path = here("outputs","maps1018","FISH_DETvsMSOMvsOCC_maps_for_eachSP1018_5factors.pdf"), msom_path = here("outputs", "msom_outputs_Viorel covs", "FISH_outputs_20241018", "occupancy_bysite_msom_FISH_20241018.csv"), occplus_path = here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "occupancybysite_occPlus_FISH_20241018_5factors.csv"), spp_list = species_lists$FISH)
)

# The main loop that iterates through each task (AmphRept, Mammal, etc.).
for (task_name in names(plotting_tasks)) {
  task <- plotting_tasks[[task_name]]
  msom_occ_sf <- read.csv(task$msom_path) %>% st_as_sf(coords = 6:7, crs = 4326, remove = FALSE)
  occplus_occ_df <- read.csv(task$occplus_path) %>% inner_join(glgs_sum2.1, by = "site")
  occplus_occ_sf <- st_as_sf(occplus_occ_df, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  
  pdf(task$pdf_path, width = 9, height = 6)
  
  for (spp in task$spp_list) {
    p1 <- plot_observed_map(spp)
    p2 <- plot_msom_map(spp, msom_occ_sf)
    p3 <- plot_occplus_map(spp, occplus_occ_sf)
    grid.arrange(p1, p3, p2, ncol = 3)
    cat("Plotted maps for:", spp, "\n")
  }
  
  dev.off()
  cat("Finished generating PDF for:", task_name, "\n\n")
}



# Plot the observed and predicted data from MSOM and occPlus only for Budorcas taxicolor-------------------------------------------------------------------------------------------------------

# Load additional spatial layers needed for these specific maps
GLG_PA <- st_read(here("data", "glgs", "GLG_PA.shp")) %>% st_transform(4326)
TBG_PA <- st_read(here("data", "glgs", "TBG_PA.shp")) %>% st_transform(4326)
XHS_PA <- st_read(here("data", "glgs", "XHS_PA.shp")) %>% st_transform(4326)
GLGS_east <- st_read(here("data", "glgs", "hybas05_GLG_East.shp")) %>% st_transform(4326)
GLGS_SW <- st_read(here("data", "glgs", "hybas05_GLG_SW.shp")) %>% st_transform(4326)
GLGS_NW <- st_read(here("data", "glgs", "hybas05_GLG_NW.shp")) %>% st_transform(4326)
GLGS_dy <- st_read(here("data", "glgs", "hybas05_GLG_DY.shp")) %>% st_transform(4326)
GLGS_lc <- st_read(here("data", "glgs", "hybas05_GLG_LC.shp")) %>% st_transform(4326)

# Load and prepare model data
# We can reuse the msom_occ_sf and occplus_occ_sf from the Mammal section of the previous script
# If not already loaded, run this:
msom_occ_sf <- read.csv(here("outputs/msom_outputs_Viorel covs/Mammal_outputs_20241018/occupancy_bysite_msom_Mammal_20241018.csv")) %>% st_as_sf(coords = 6:7, crs = 4326)
occplus_occ_df <- read.csv(here("outputs/occPlus_outputs_Viorel covs/Mammal_outputs_20241018_4factors/occupancybysite_occPlus_Mammal_20241018_4factors.csv")) %>% inner_join(glgs_sum2.1, by = "site")
occplus_occ_sf <- st_as_sf(occplus_occ_df, coords = c("longitude", "latitude"), crs = 4326)

# Load Conditional Occupancy data
occplus_COPocc_df <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Mammal_outputs_20241018_4factors", "COPoccupancybysite_occPlus_Mammal_20241018_4factors.csv")) %>% inner_join(glgs_sum2.1, by = "site")
occplus_COPocc_sf <- st_as_sf(occplus_COPocc_df, coords = c("longitude", "latitude"), crs = 4326)

# Define the size mapping function for Budorcas taxicolor specifically
# This resolves the "'size_mapping' not found" error
budorcas_size_mapping <- function(value) {
  cut(value, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = 1:5, right = TRUE)
}

# Function for the customized "Observed" map base (with 3 PAs)
create_budorcas_observed_base <- function() {
  ggplot() +
    geom_sf(data = GLGS, color = "darkblue", fill = "lightblue", alpha = 0.8) +
    geom_sf(data = GLG_PA, color = "darkgreen", fill = "darkgreen", alpha = 0.8) + 
    geom_sf(data = TBG_PA, color = "deeppink", fill = "deeppink", alpha = 1) + 
    geom_sf(data = XHS_PA, color = "dodgerblue", fill = "dodgerblue", alpha = 1) +
    theme(
      legend.position = c(0.35, 0.65),
      legend.title = element_text(face = "bold", hjust = 0, size = 12, lineheight = 1.1),
      legend.spacing = unit(0, 'cm'),
      legend.background = element_rect(fill = "transparent", color = NA),
      axis.title = element_blank()
    ) +
    geom_tile(data = data.frame(x = 97.25, y = 25.7), aes(x, y), fill = "dodgerblue", width = 0.22, height = 0.14, alpha = 0.8) +
    annotate("text", x = 97.4, y = 25.7, label = "Xiaoheishan PA", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.25, y = 25.9), aes(x, y), fill = "deeppink", width = 0.22, height = 0.14, alpha = 0.8) +
    annotate("text", x = 97.4, y = 25.9, label = "Tongbiguan PA", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.25, y = 26.1), aes(x, y), fill = "darkgreen", width = 0.22, height = 0.14, alpha = 0.8) +
    annotate("text", x = 97.4, y = 26.1, label = "Gaoligongshan PA", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.25, y = 26.3), aes(x, y), fill = "lightblue", width = 0.22, height = 0.14, alpha = 0.6, color = "blue") +
    annotate("text", x = 97.4, y = 26.3, label = "Gaoligong Region", hjust = 0, size = 3)
}

# Function for the customized MSOM map base (with river basins)
create_budorcas_msom_base <- function() {
  ggplot() +
    geom_sf(data = GLGS_east, color = "dodgerblue4", fill = "palegreen", alpha = 0.8) + 
    geom_sf(data = GLGS_NW, color = "dodgerblue4", fill = "paleturquoise", alpha = 0.8) + 
    geom_sf(data = GLGS_dy, color = "dodgerblue4", fill = "rosybrown", alpha = 0.8) + 
    geom_sf(data = GLGS_lc, color = "dodgerblue4", fill = "moccasin", alpha = 0.8) +
    theme(legend.position = c(0.25, 0.72), 
          legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
          legend.background = element_rect(fill = "transparent", color = NA)) +
    geom_tile(data = data.frame(x = 97.15, y = 25.7), aes(x, y), fill = "rosybrown", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 25.7, label = "Daying River basin", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.15, y = 25.9), aes(x, y), fill = "moccasin", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 25.9, label = "Longchuan River basin", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.15, y = 26.1), aes(x, y), fill = "palegreen", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 26.1, label = "Nu River basin", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.15, y = 26.3), aes(x, y), fill = "paleturquoise", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 26.3, label = "Dulong River basin", hjust = 0, size = 3)
}

# Function for the customized OccPlus map base (with subregions)
create_budorcas_occplus_base <- function() {
  ggplot() +
    geom_sf(data = GLGS_east, color = "dodgerblue4", fill = "palegreen", alpha = 0.8) + 
    geom_sf(data = GLGS_NW, color = "dodgerblue4", fill = "paleturquoise", alpha = 0.8) + 
    geom_sf(data = GLGS_SW, color = "dodgerblue4", fill = "peachpuff", alpha = 0.8) +
    theme(legend.position = c(0.32, 0.68), 
          legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
          legend.background = element_rect(fill = "transparent", color = NA)) +
    geom_tile(data = data.frame(x = 97.15, y = 25.7), aes(x, y), fill = "peachpuff", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 25.7, label = "Southwest Subregion", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.15, y = 25.9), aes(x, y), fill = "palegreen", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 25.9, label = "East Subregion", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.15, y = 26.1), aes(x, y), fill = "paleturquoise", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 26.1, label = "Northwest Subregion", hjust = 0, size = 3)
}

# Plot 1: Observed Detections (reusing the function from the previous script is difficult due to the custom base map)
p1 <- create_budorcas_observed_base() +
  geom_sf(data = trial %>% filter(Budorcas.taxicolor == 0),
          shape = 21, fill = NA, color = "brown", size = 4, alpha = 0.7) +
  geom_sf(data = trial %>% filter(Budorcas.taxicolor > 0),
          aes(fill = Budorcas.taxicolor), shape = 21, color = "brown", size = 4, alpha = 0.7) +
  scale_fill_gradient(low = "yellow", high = "red", name = "Detections of\nBudorcas taxicolor", limits = c(1, 18)) +
  theme(legend.position = c(0.25, 0.72), legend.title = element_text(face = "bold", hjust = 0.5, size = 10),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 6)) +
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE, default_crs = st_crs(4326))

# Plot 2: MSOM Predictions
msom_occ_sf$size_mapped <- as.numeric(budorcas_size_mapping(msom_occ_sf[["Mean.Budorcas.taxicolor"]]))
p2 <- create_budorcas_msom_base() +
  geom_sf(data = msom_occ_sf, aes(size = size_mapped), shape = 21, color = "black", fill = 'purple', alpha = 0.6) +
  labs(size = "MSOM\nOccupancy probablity") +
  theme(legend.position = c(0.32, 0.72), legend.title = element_text(face = "bold", hjust = 0.5, size = 10),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  scale_size_continuous(limits = c(1, 5), breaks = 1:5, labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE, default_crs = st_crs(4326))

# Plot 3: OccPlus Predictions
occplus_occ_sf$size_mapped <- as.numeric(budorcas_size_mapping(occplus_occ_sf[["Mean.Budorcas.taxicolor"]]))
p3 <- create_budorcas_occplus_base() +
  geom_sf(data = occplus_occ_sf, aes(size = size_mapped), shape = 21, color = "darkred", fill = 'orange', alpha = 0.6) +
  labs(size = "OccPlus\nOccupancy probablity") +
  theme(legend.position = c(0.32, 0.7), legend.title = element_text(face = "bold", hjust = 0.5, size = 10),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  scale_size_continuous(limits = c(1, 5), breaks = 1:5, labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE, default_crs = st_crs(4326))

# Plot 4: OccPlus Conditional Predictions
occplus_COPocc_sf$size_mapped <- as.numeric(budorcas_size_mapping(occplus_COPocc_sf[["Mean.Budorcas.taxicolor"]]))
p4 <- create_budorcas_occplus_base() +
  geom_sf(data = occplus_COPocc_sf, aes(size = size_mapped), shape = 21, color = "darkred", fill = 'orange', alpha = 0.6) +
  labs(size = "OccPlus\nConditional posterior probability") +
  theme(legend.position = c(0.32, 0.7), legend.title = element_text(face = "bold", hjust = 0.5, size = 10),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  scale_size_continuous(limits = c(1, 5), breaks = 1:5, labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE, default_crs = st_crs(4326))

# Arrange and save the two PDF files
pdf(here("outputs", "maps1018", "Budorcas_taxicolor_DETvsMSOMvsOCC_maps.pdf"), width = 9, height = 6)
grid.arrange(p1, p3, p2, ncol = 3)
dev.off()

pdf(here("outputs", "maps1018", "Budorcas_taxicolor_Conditional.pdf"), width = 9, height = 6)
grid.arrange(p1, p4, p3, ncol = 3)
dev.off()



# Plot the observed and predicted data for Rana catesbeiana -------------------------------------------------------------------------------------------------------
# --- Data Preparation for Bullfrog Maps ---

# Load spatial data (if not already loaded)
# GLGS <- st_read(here("data", "glgs", "GLG_Region_FULL_UTM.shp")) %>% st_transform(4326)
# PA <- st_read(here("data", "glgs", "GLGS_PAs.shp")) %>% st_transform(4326)

# Prepare observed data specifically for Rana catesbeiana
# This uses the raw 'glgs' dataframe and summarizes detections for the specific OTU.
binarise1 <- function(x, na.rm = FALSE) vegan::decostand(x, method = "pa") 
glgs_bullfrog_obs <- glgs %>% 
  mutate(across(starts_with("Amphibia"), binarise1)) %>% 
  summarise(
    across(Filter_vol:DistFromPAedge, first),
    across(starts_with("Amphibia"), sum), 
    .by = c(site)
  )  %>% 
  filter(status == "sample")

# Create the sf object for observed data
trial_bullfrog <- st_as_sf(glgs_bullfrog_obs, coords = c("longitude", "latitude"), crs = 4326)

# Load and prepare model prediction data
msom_bullfrog_sf <- read.csv(here("outputs", "output_bullfrog", "occupancy_preds_MSOM_Bullfrogs.csv")) %>% 
  left_join(glgs_sum2.1[, c("site", "longitude", "latitude")], by = c("obs.covs1.site" = "site")) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

occplus_bullfrog_sf <- read.csv(here("outputs","output_bullfrog","occupancybysite_occPlus_bullfrog_20250217.csv")) %>% 
  left_join(glgs_sum2.1[, c("site", "longitude", "latitude")], by = "site") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Define the specific OTU column name for observed data
bullfrog_otu_col <- "Amphibia_Anura_Ranidae_Rana_catesbeiana_OTU414_31973"


# --- Reuse and Adapt Plotting Functions ---

# We can reuse the generic 'prob_size_mapping' function from the previous section.
# We will also reuse the 'create_base_map' function logic inside our new plot calls.

# For this one-off case, we can create the plots directly without new functions,
# but we will follow the clean logic established before.

# Plot 1: Observed Detections for Bullfrog
p1_bullfrog <- create_base_map() +
  geom_sf(data = trial_bullfrog %>% filter(!!sym(bullfrog_otu_col) == 0),
          shape = 21, fill = NA, color = "brown", size = 4, alpha = 0.7) +
  geom_sf(data = trial_bullfrog %>% filter(!!sym(bullfrog_otu_col) > 0),
          aes(fill = !!sym(bullfrog_otu_col)),
          shape = 21, color = "brown", size = 4, alpha = 0.7) +
  scale_fill_gradient(low = "yellow", high = "red", 
                      name = "Detection of\nRana catesbeiana", 
                      limits = c(1, 18)) +
  theme(legend.position = c(0.25, 0.72), 
        legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 6)) +
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE, default_crs = st_crs(4326))


# Plot 2: Single-Species Occupancy Model Predictions
msom_bullfrog_sf$size_mapped <- as.numeric(prob_size_mapping(msom_bullfrog_sf$psi.mean))
p2_bullfrog <- create_base_map() +
  geom_sf(data = msom_bullfrog_sf, aes(size = size_mapped),
          shape = 21, color = "black", fill = 'purple', alpha = 0.6) +
  labs(size = "Single-Species Occupancy Model\nOccupancy probablity\nRana catesbeiana") +
  theme(legend.position = c(0.32, 0.68), 
        legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  scale_size_continuous(limits = c(1, 5), breaks = 1:5, 
                        labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE, default_crs = st_crs(4326))


# Plot 3: OccPlus Predictions
occplus_bullfrog_sf$size_mapped <- as.numeric(prob_size_mapping(occplus_bullfrog_sf$Mean.Rana.catesbeiana))
p3_bullfrog <- create_base_map() +
  geom_sf(data = occplus_bullfrog_sf, aes(size = size_mapped),
          shape = 21, color = "darkred", fill = 'orange', alpha = 0.6) +
  labs(size = "OccPlus\nPredictive posterior probability\nRana catesbeiana") +
  theme(legend.position = c(0.32, 0.68), 
        legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  scale_size_continuous(limits = c(1, 5), breaks = 1:5, 
                        labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE, default_crs = st_crs(4326))


# --- Arrange and Save the Final PDF ---
library(gridExtra)

pdf(here("outputs","output_bullfrog","Bullfrog_DETvsSSPOCC_maps20250217.pdf"), width = 9, height = 6)
grid.arrange(p1_bullfrog, p3_bullfrog, p2_bullfrog, ncol = 3)
dev.off()

