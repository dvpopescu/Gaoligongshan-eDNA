library(here)
library(tidyverse)


# READ DATA ---------------------------------------------------------------



filename <- "OTUtable_12S_toSP_GLG23_20240528.txt"

data <- read_tsv(file = here("data", filename))
str(data)

new_covs <- read_tsv(here("data", "23GLG_covs_info_20240927.txt"))
zones <- read_tsv(here("data", "23GLG_zones_info_20240903.txt"))
climate_covs <- read_csv(here("data","GLGS_climate.csv"))
climate_covs <- climate_covs %>% select(site, annual_precipitation, mean_temperature)

new_covs <- left_join(new_covs, zones, by = "site")
new_covs <- left_join(new_covs, climate_covs, by = "site")

str(new_covs)

binarise1 <- function(x, na.rm = FALSE) vegan::decostand(x, method = "pa") 


# SPECIES TO BE REMOVED FROM DATASET --------------------------------------

# RETAIN ONLY SPECIES THAT WERE DETECTED AT >1 SITE

# AmphRept species
spptoremove_AmphRept <- data %>% 
  as.data.frame() %>%
  #dplyr::select(!contains(c("Spikein","Actinopteri", "Aves", "Amphibia", "Reptilia", "Teleostei"))) %>% # extract just a taxon data
  dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Mammalia", "Aves"))) %>% 
  dplyr::select(SampleID, PCR, site, contains("OTU")) %>% 
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

# Mammal species
spptoremove_Mammal <- data %>% 
  as.data.frame() %>%
  dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Aves", "Amphibia", "Reptilia"))) %>% 
  dplyr::select(SampleID, PCR, site, contains("OTU")) %>% 
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

# Ave species
spptoremove_Ave <- data %>% 
  as.data.frame() %>%
  dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Mammalia", "Amphibia", "Reptilia"))) %>% 
  dplyr::select(SampleID, PCR, site, contains("OTU")) %>% 
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

# Fish species
spptoremove_fish <- data %>% 
  as.data.frame() %>%
  dplyr::select(!contains(c("Spikein","Mammalia","Amphibia", "Reptilia", "Aves"))) %>%
  dplyr::select(!contains(c("Cololabis_saira", "Sardina_pilchardus",
                            "Sardinella_lemuru", "Engraulis_ringens",
                            "Micromesistius_poutassou", "Scophthalmus_maximus",
                            "Scomberomorus_niphonius", "Trichiurus_haumela"))) %>%
  dplyr::select(SampleID, PCR, site, contains("OTU")) %>% 
  #dplyr::rename_with(.cols = starts_with(c("Actinopteri", "Teleostei")), 
  #                   function(x){paste0("OTU_", x)}) |>
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

# create OTU table
OTUtable_AmphRept <- data %>%
  as.data.frame() %>%
  # FOR AmphRept SPECIES DATA
  dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Mammalia", "Aves"))) %>% 
  dplyr::select(SampleID, PCR, site, contains("OTU")) %>% 
  group_by(SampleID) %>% 
  summarise(across(contains("OTU"), function(x){
    sum(x > 0)
  })) %>% 
  arrange(SampleID) %>% 
  dplyr::select(contains("OTU")) %>% 
  dplyr::select(-one_of(spptoremove_AmphRept)) %>%
  dplyr::filter(!row_number() %in% c(304:321)) %>%
  dplyr::select(!contains(c("Rana_catesbeiana"))) %>%  #（remove bullfrog）
  as.matrix

OTUtable_mammal <- data %>% 
  as.data.frame() %>%
  # FOR AmphRept SPECIES DATA
  dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Aves", "Amphibia", "Reptilia"))) %>% 
#  dplyr::select(!contains(c("Bos_gaurus", "Bos_grunniens", "Bos_taurus", "Bubalus_bubalis", 
#                     "Capra_hircus", "Ovis_aries", "Sus_scrofa", "Canis_lupus", 
#                     "Felis_catus", "Oryctolagus_cuniculus", "Equus_asinus", "Equus_caballus"))) %>% 
  dplyr::select(SampleID, PCR, site, contains("OTU")) %>% 
  group_by(SampleID) %>% 
  summarise(across(contains("OTU"), function(x){
    sum(x > 0)
  })) %>% 
  arrange(SampleID) %>% 
  dplyr::select(contains("OTU")) %>% 
  dplyr::select(-one_of(spptoremove_Mammal)) %>%
  dplyr::filter(!row_number() %in% c(304:321)) %>%
  as.matrix

OTUtable_Ave <- data %>% 
  as.data.frame() %>%
  # FOR Ave SPECIES DATA
  dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Mammalia", "Amphibia", "Reptilia"))) %>% 
  dplyr::select(SampleID, PCR, site, contains("OTU")) %>% 
  group_by(SampleID) %>% 
  summarise(across(contains("OTU"), function(x){
    sum(x > 0)
  })) %>% 
  arrange(SampleID) %>% 
  dplyr::select(contains("OTU")) %>% 
  dplyr::select(-one_of(spptoremove_Ave)) %>%
  dplyr::filter(!row_number() %in% c(304:321)) %>%
  as.matrix

OTUtable_fish <- data %>% 
  as.data.frame() %>%
  # FOR FISH DATA (remove all spp other than fish + seafish)
  dplyr::select(!contains(c("Spikein","Mammalia","Amphibia", "Reptilia", "Aves"))) %>%
  dplyr::select(!contains(c("Cololabis_saira", "Sardina_pilchardus",
                            "Sardinella_lemuru", "Engraulis_ringens",
                            "Micromesistius_poutassou", "Scophthalmus_maximus",
                            "Scomberomorus_niphonius", "Trichiurus_haumela"))) %>%
  dplyr::select(SampleID, PCR, site, contains("OTU")) %>% 
  group_by(SampleID) %>% 
  summarise(across(contains("OTU"), function(x){
    sum(x > 0)
  })) %>% 
  arrange(SampleID) %>% 
  dplyr::select(contains("OTU")) %>% 
  dplyr::select(-one_of(spptoremove_fish)) %>%
  dplyr::filter(!row_number() %in% c(304:321)) %>%
  as.matrix


# create object with environmental and detection covariates
data_infos <- data %>% 
  arrange(SampleID) %>% 
  group_by(SampleID) %>%
  dplyr::filter(row_number() == 1) %>% 
  mutate(streamOrigin = case_when(
    is.na(DistFromPAedge) ~ "outsideGLGS",
    TRUE ~ "insideGLGS"
  )) %>% 
  rename(Site = site) %>% 
  dplyr::select(Site, latitude, longitude, river, tributary, team,
                streamOrigin, DistFromPAedge, altitude, TEMP, river_width, water_velocity, 
                Filter_vol, PrevDayRain) %>%
  as.data.frame %>%
  dplyr::filter(!row_number() %in% c(304:321))


land.covs <- new_covs |>
  rename(Site = site) |>
  as.data.frame() 

data_infos1 <- left_join(data_infos, land.covs, by = "Site") |>
  as.data.frame() 

summary(data_infos1)

#ggplot(data_infos, aes(x = longitude, y = latitude)) + geom_point()

# combine all data in a list for analysis and save it
data <- list(
  "info" = data_infos1,
  "OTU" = OTUtable_AmphRept,
  #"OTU" = OTUtable_mammal,
  #"OTU" = OTUtable_Ave,
  #"OTU" = OTUtable_fish,
  "K" = 6) # number of PCR replicates per water sample

save(data, file = here("data", "gaoligongshan_AmphRept_data_20241014.rda"))
save(data, file = here("data", "gaoligongshan_Mammal_data_20240912.rda"))
save(data, file = here("data", "gaoligongshan_Ave_data_20240912.rda"))
save(data, file = here("data", "gaoligongshan_FISH_data_20240912.rda"))
