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


# create OTU table
OTUtable_bullfrog <- data %>%
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
  dplyr::select(contains(c("Rana_catesbeiana"))) %>% #（only keep bullfrog）
  as.vector()


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


# combine all data in a list for analysis and save it
data <- list(
  "info" = data_infos1,
  "OTU" = OTUtable_bullfrog,
  "K" = 6) # number of PCR replicates per water sample

save(data, file = here("data", "gaoligongshan_bullfrog_data_20250127.rda"))
