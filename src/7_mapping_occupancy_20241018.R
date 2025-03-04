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
spp_list <- read_tsv(here("data","SP_INFO_20240528.txt"))
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

# Preparing the required objects
GLGS <- st_read(here("data", "glgs", "GLG_Region_FULL_UTM.shp"))
PA <- st_read(here("data", "glgs", "GLGS_PAs.shp"))
raster_map <- raster(here("data", "glgs", "forest_1km.tif"))
raster_df <- as.data.frame(raster_map, xy = TRUE)
GLGS <- st_transform(GLGS, 4326)
PA <- st_transform(PA, 4326)


# for TERR
msom_occ_AmphRept <- read.csv(here("outputs", "msom_outputs_Viorel covs", "AmphRept_outputs_20241018",
                                   "occupancy_bysite_msom_AmphRept_20241018.csv"))
occPlus_occ_AmphRept <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors",
                                      "occupancybysite_occPlus_AmphRept_20241018_3factors.csv"))
msom_occ_Mammal <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Mammal_outputs_20241018",
                                 "occupancy_bysite_msom_Mammal_20241018.csv"))
occPlus_occ_Mammal <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Mammal_outputs_20241018_4factors",
                                    "occupancybysite_occPlus_Mammal_20241018_4factors.csv"))
msom_occ_Ave <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Ave_outputs_20241018",
                              "occupancy_bysite_msom_Ave_20241018.csv"))
occPlus_occ_Ave <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Ave_outputs_20241018_2factors",
                                 "occupancybysite_occPlus_Ave_20241018_2factors.csv"))
msom_occ <- left_join(msom_occ_AmphRept, msom_occ_Mammal, by = "site") %>% 
  left_join(msom_occ_Ave, by = "site") %>% 
  dplyr::select(!contains(c("Bos.gaurus","Bos.grunniens","Bos.taurus","Bubalus.bubalis","Capra.hircus","Ovis.aries","Cervus.nippon","Sus.scrofa","Canis.lupus","Felis.catus","Oryctolagus.cuniculu","Equus.asinus","Equus.caballus","Anas.platyrhynchos","Gallus.gallus","Melopsittacus.undulatus","Pelodiscus.sinensis","Trachemys.scripta","Rana.catesbeiana"))) %>% 
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))
occPlus_occ <- left_join(occPlus_occ_AmphRept, occPlus_occ_Mammal, by = "site") %>% 
  left_join(occPlus_occ_Ave, by = "site")
occPlus_occ <- inner_join(glgs_sum2.1[,1:21], occPlus_occ, by = "site")
occPlus_occ <- occPlus_occ %>% 
  dplyr::select(!contains(c("Bos.gaurus","Bos.grunniens","Bos.taurus","Bubalus.bubalis","Capra.hircus","Ovis.aries","Cervus.nippon","Sus.scrofa","Canis.lupus","Felis.catus","Oryctolagus.cuniculu","Equus.asinus","Equus.caballus","Anas.platyrhynchos","Gallus.gallus","Melopsittacus.undulatus","Pelodiscus.sinensis","Trachemys.scripta","Rana.catesbeiana"))) %>% 
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))
msom_occ_sf <- st_as_sf(msom_occ,coords = 6:7) %>% st_set_crs(4326)
occPlus_occ_sf <- st_as_sf(occPlus_occ,coords = 3:4) %>% st_set_crs(4326)
range(glgs_sum2.1_sf$terr_rich) #[1]  7 73
range(glgs_sum2.1_sf$terr_rich3) #[1]  1 49
range(msom_occ_sf$summed_probs) #[1] 26.73645 78.80499
range(occPlus_occ$summed_probs) #[1]  2.734734 51.739717
size_mapping <- function(value) {
  if (value <= 17) {
    return(1)
  } else if (value <= 35) {
    return(2)
  } else if (value <= 53) {
    return(3)
  } else if (value <= 71) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("1-17", "18-35", "36-53", "54-71", "72-86"),
)

# for FISH
msom_occ <- read.csv(here("outputs", "msom_outputs_Viorel covs", "FISH_outputs_20241018",
                          "occupancy_bysite_msom_FISH_20241018.csv"))
occPlus_occ <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors",
                             "occupancybysite_occPlus_FISH_20241018_5factors.csv"))
msom_occ <- msom_occ %>% 
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))
occPlus_occ <- inner_join(glgs_sum2.1[,1:21], occPlus_occ, by = "site")
occPlus_occ <- occPlus_occ %>% 
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))
msom_occ_sf <- st_as_sf(msom_occ,coords = 6:7) %>% st_set_crs(4326)
occPlus_occ_sf <- st_as_sf(occPlus_occ,coords = 3:4) %>% st_set_crs(4326)
range(glgs_sum2.1_sf$fish_rich) #[1]  1 47
range(glgs_sum2.1_sf$fish_rich3) #[1]  0 35
range(msom_occ_sf$summed_probs) #[1]  9.071901 47.706529
range(occPlus_occ$summed_probs) #[1]  1.272384 31.982179
size_mapping <- function(value) {
  if (value <= 10) {
    return(1)
  } else if (value <= 20) {
    return(2)
  } else if (value <= 30) {
    return(3)
  } else if (value <= 40) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("0-10", "11-20", "21-30", "31-40", "41-48"),
)

# for native and invasive fish
msom_occ <- read.csv(here("outputs", "msom_outputs_Viorel covs", "FISH_outputs_20241018",
                          "occupancy_bysite_msom_FISH_20241018.csv"))
occPlus_occ <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors",
                             "occupancybysite_occPlus_FISH_20241018_5factors.csv"))
spp_list_native <- spp_list_native %>% 
  mutate(spp = gsub(" ", ".", species)) %>% 
  mutate(spp = paste0("Mean.", spp))
spp_list_invasive <- spp_list_invasive %>% 
  mutate(spp = gsub(" ", ".", species)) %>% 
  mutate(spp = paste0("Mean.", spp))
col_names <- colnames(msom_occ)
fish_native <- which(col_names %in% spp_list_native$spp)
msom_occ$fish_native_rich <- rowSums(msom_occ[fish_native])
fish_nonnative <- which(col_names %in% spp_list_invasive$spp)
msom_occ$fish_nonnative_rich <- rowSums(msom_occ[fish_nonnative])
col_names <- colnames(occPlus_occ)
fish_native <- which(col_names %in% spp_list_native$spp)
occPlus_occ$fish_native_rich <- rowSums(occPlus_occ[fish_native])
fish_nonnative <- which(col_names %in% spp_list_invasive$spp)
occPlus_occ$fish_nonnative_rich <- rowSums(occPlus_occ[fish_nonnative])
occPlus_occ <- inner_join(glgs_sum2.1[,1:21], occPlus_occ, by = "site")
msom_occ_sf <- st_as_sf(msom_occ,coords = 6:7) %>% st_set_crs(4326)
occPlus_occ_sf <- st_as_sf(occPlus_occ,coords = 3:4) %>% st_set_crs(4326)

range(glgs_sum2.1$fish_native_rich) # [1]  0 23
range(glgs_sum2.1$fish_native_rich3) # [1]  0 17
range(msom_occ$fish_native_rich) #[1]  4.030184 20.176333
range(occPlus_occ$fish_native_rich) #[1]  0.7027038 12.6810212
range(glgs_sum2.1$fish_nonnative_rich) # [1]  0 28
range(glgs_sum2.1$fish_nonnative_rich3) # [1]  0 21
range(msom_occ$fish_nonnative_rich) #[1]  4.384914 27.144421
range(occPlus_occ$fish_nonnative_rich) #[1]  0.3975384 19.8265722
size_mapping <- function(value) {
  if (value <= 4) {
    return(1)
  } else if (value <= 9) {
    return(2)
  } else if (value <= 14) {
    return(3)
  } else if (value <= 19) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("0-4", "5-9", "10-14", "15-19", "≥20"),
)

# for mammal
msom_occ <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Mammal_outputs_20241018",
                          "occupancy_bysite_msom_Mammal_20241018.csv"))
occPlus_occ <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Mammal_outputs_20241018_4factors",
                             "occupancybysite_occPlus_Mammal_20241018_4factors.csv"))
msom_occ <- msom_occ %>% 
  dplyr::select(!contains(c("Bos.gaurus","Bos.grunniens","Bos.taurus","Bubalus.bubalis","Capra.hircus","Ovis.aries","Cervus.nippon","Sus.scrofa","Canis.lupus","Felis.catus","Oryctolagus.cuniculu","Equus.asinus","Equus.caballus","Anas.platyrhynchos","Gallus.gallus","Melopsittacus.undulatus","Pelodiscus.sinensis","Trachemys.scripta","Rana.catesbeiana"))) %>%
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))
occPlus_occ <- inner_join(glgs_sum2.1[,1:21], occPlus_occ, by = "site")
occPlus_occ <- occPlus_occ %>% 
  dplyr::select(!contains(c("Bos.gaurus","Bos.grunniens","Bos.taurus","Bubalus.bubalis","Capra.hircus","Ovis.aries","Cervus.nippon","Sus.scrofa","Canis.lupus","Felis.catus","Oryctolagus.cuniculu","Equus.asinus","Equus.caballus","Anas.platyrhynchos","Gallus.gallus","Melopsittacus.undulatus","Pelodiscus.sinensis","Trachemys.scripta","Rana.catesbeiana"))) %>%
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))
msom_occ_sf <- st_as_sf(msom_occ,coords = 6:7) %>% st_set_crs(4326)
occPlus_occ_sf <- st_as_sf(occPlus_occ,coords = 3:4) %>% st_set_crs(4326)
range(glgs_sum2.1$mammal_rich) # [1]  3 36
range(glgs_sum2.1$mammal_rich3) # [1]  0 27
range(msom_occ$summed_probs) #[1] 12.78687 34.67871
range(occPlus_occ$summed_probs) #[1]  1.021445 26.971971
size_mapping <- function(value) { # for mammal
  if (value <= 7) {
    return(1)
  } else if (value <= 14) {
    return(2)
  } else if (value <= 21) {
    return(3)
  } else if (value <= 28) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("1-7", "8-14", "15-21", "22-28", "29-36"),
)

# for Ave
msom_occ <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Ave_outputs_20241018",
                          "occupancy_bysite_msom_Ave_20241018.csv"))
occPlus_occ <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Ave_outputs_20241018_2factors",
                             "occupancybysite_occPlus_Ave_20241018_2factors.csv"))
msom_occ <- msom_occ %>% 
  dplyr::select(!contains(c("Bos.gaurus","Bos.grunniens","Bos.taurus","Bubalus.bubalis","Capra.hircus","Ovis.aries","Cervus.nippon","Sus.scrofa","Canis.lupus","Felis.catus","Oryctolagus.cuniculu","Equus.asinus","Equus.caballus","Anas.platyrhynchos","Gallus.gallus","Melopsittacus.undulatus","Pelodiscus.sinensis","Trachemys.scripta","Rana.catesbeiana"))) %>%
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))
occPlus_occ <- inner_join(glgs_sum2.1[,1:21], occPlus_occ, by = "site")
occPlus_occ <- occPlus_occ %>% 
  dplyr::select(!contains(c("Bos.gaurus","Bos.grunniens","Bos.taurus","Bubalus.bubalis","Capra.hircus","Ovis.aries","Cervus.nippon","Sus.scrofa","Canis.lupus","Felis.catus","Oryctolagus.cuniculu","Equus.asinus","Equus.caballus","Anas.platyrhynchos","Gallus.gallus","Melopsittacus.undulatus","Pelodiscus.sinensis","Trachemys.scripta","Rana.catesbeiana"))) %>%
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))
msom_occ_sf <- st_as_sf(msom_occ,coords = 6:7) %>% st_set_crs(4326)
occPlus_occ_sf <- st_as_sf(occPlus_occ,coords = 3:4) %>% st_set_crs(4326)
range(glgs_sum2.1$ave_rich) # [1]  2 27
range(glgs_sum2.1$ave_rich3) # [1]  0 16
range(msom_occ$summed_probs) #[1] 7.274817 30.401589
range(occPlus_occ$summed_probs) #[1]  0.3006476 14.9181595
size_mapping <- function(value) { # for bird
  if (value <= 6) {
    return(1)
  } else if (value <= 12) {
    return(2)
  } else if (value <= 18) {
    return(3)
  } else if (value <= 24) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("0-6", "7-12", "13-18", "19-24", "25-31"),
)

# for amphrept
msom_occ <- read.csv(here("outputs", "msom_outputs_Viorel covs", "AmphRept_outputs_20241018",
                          "occupancy_bysite_msom_AmphRept_20241018.csv"))
occPlus_occ <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors",
                             "occupancybysite_occPlus_AmphRept_20241018_3factors.csv"))
msom_occ <- msom_occ %>% 
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))
occPlus_occ <- inner_join(glgs_sum2.1[,1:21], occPlus_occ, by = "site")
occPlus_occ <- occPlus_occ %>% 
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))
msom_occ_sf <- st_as_sf(msom_occ,coords = 6:7) %>% st_set_crs(4326)
occPlus_occ_sf <- st_as_sf(occPlus_occ,coords = 3:4) %>% st_set_crs(4326)
range(glgs_sum2.1$amph_rich) # [1]  0 21
range(glgs_sum2.1$amph_rich3) # [1]  0 14
range(msom_occ$summed_probs) #[1] 1.671985 15.991755
range(occPlus_occ$summed_probs) #[1]  0.8010411 13.6949536
size_mapping <- function(value) { # for amphrept
  if (value <= 4) {
    return(1)
  } else if (value <= 8) {
    return(2)
  } else if (value <= 12) {
    return(3)
  } else if (value <= 16) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("0-4", "5-8", "9-12", "13-16", "17-21"),
)


glgs_sum2.1_sf <- st_as_sf(glgs_sum2.1,coords = 3:4) %>% st_set_crs(4326)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$terr_rich, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$terr_rich3, size_mapping)
msom_occ_sf$size_mapped <- sapply(msom_occ_sf$summed_probs, size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf$summed_probs, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$fish_rich, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$fish_rich3, size_mapping)
msom_occ_sf$size_mapped <- sapply(msom_occ_sf$summed_probs, size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf$summed_probs, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$fish_native_rich, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$fish_native_rich3, size_mapping)
msom_occ_sf$size_mapped <- sapply(msom_occ_sf$fish_native_rich, size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf$fish_native_rich, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$fish_nonnative_rich, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$fish_nonnative_rich3, size_mapping)
msom_occ_sf$size_mapped <- sapply(msom_occ_sf$fish_nonnative_rich, size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf$fish_nonnative_rich, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$mammal_rich, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$ave_rich, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$amph_rich, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$mammal_rich3, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$ave_rich3, size_mapping)
glgs_sum2.1_sf$size_mapped <- sapply(glgs_sum2.1_sf$amph_rich3, size_mapping)
msom_occ_sf$size_mapped <- sapply(msom_occ_sf$summed_probs, size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf$summed_probs, size_mapping)

ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = forest_1km)) +
  scale_fill_gradient(low = "white", high = "grey70", na.value = NA) +
  guides(fill = "none") +
  geom_sf(data = GLGS, color = "darkblue", fill = NA, linetype = "solid", size = 8) +
  geom_sf(data = PA, color = "green", fill = "green", alpha = 0.8) +
  geom_sf(data = glgs_sum2.1_sf, 
  #geom_sf(data = msom_occ_sf, 
  #geom_sf(data = occPlus_occ_sf, 
          aes(size = size_mapped, alpha = riverOrigin), 
          #shape = 21, color = "black", fill = "brown", # terr
          #shape = 21, color = "black", fill = "blue", # fish
          #shape = 21, color = "black", fill = 'royalblue', # for native fish
          #shape = 21, color = "brown", fill = 'red', # for non-native fish
          #shape = 21, color = "black", fill = "palevioletred4", # mam
          #shape = 21, color = "black", fill = "deepskyblue4", # ave
          shape = 21, color = "black", fill = "blueviolet", # amphrept
          show.legend = 'point', inherit.aes = FALSE) +
  #labs(size = paste("Observed terrestrial\n", "species richness\n", "at ≥1/18 detections (PCRs)")) +
  #labs(size = paste("Observed terrestrial\n", "species richness\n", "at ≥3/18 detections (PCRs)")) +
  #labs(size = paste("Predicted terrestrial\n", "species richness (MSOM)\n","[summed probabilities]")) +
  #labs(size = paste("Predicted terrestrial\n", "species richness (OccPlus)\n","[summed probabilities]")) +
  #labs(size = paste("Observed fish\n", "species richness\n", "at ≥1/18 detections (PCRs)")) +
  #labs(size = paste("Observed fish\n", "species richness\n", "at ≥3/18 detections (PCRs)")) +
  #labs(size = paste("Predicted fish\n", "species richness (MSOM)\n","[summed probabilities]")) +
  #labs(size = paste("Predicted fish\n", "species richness (OccPlus)\n","[summed probabilities]")) +
  #labs(size = paste("Observed native fish\n", "species richness\n", "at ≥1/18 detections (PCRs)")) +
  #labs(size = paste("Observed native fish\n", "species richness\n", "at ≥3/18 detections (PCRs)")) +
  #labs(size = paste("Predicted Native Fish\n", "species richness (MSOM)\n"," [summed probabilities]")) +
  #labs(size = paste("Predicted Native Fish\n", "species richness (OccPlus)\n","[summed probabilities]")) +
  #labs(size = paste("Observed non-native fish\n", "species richness\n", "at ≥1/18 detections (PCRs)")) +
  #labs(size = paste("Observed non-native fish\n", "species richness\n", "at ≥3/18 detections (PCRs)")) +
  #labs(size = paste("Predicted Non-native Fish\n", "species richness (MSOM)\n", "[summed probabilities]")) +
  #labs(size = paste("Predicted Non-native Fish\n", "species richness (OccPlus)\n", "[summed probabilities]")) +
  #labs(size = paste("Observed Mammalia\n", "species richness\n", "at ≥1/18 detections (PCRs)")) +
  #labs(size = paste("Observed Mammalia\n", "species richness\n", "at ≥3/18 detections (PCRs)")) +
  #labs(size = paste("Predicted Mammalia\n", "species richness (MSOM)\n","[summed probabilities]")) +
  #labs(size = paste("Predicted Mammalia\n", "species richness (OccPlus)\n","[summed probabilities]")) +
  #labs(size = paste("Observed Aves\n", "species richness\n", "at ≥1/18 detections (PCRs)")) +
  #labs(size = paste("Observed Aves\n", "species richness\n", "at ≥3/18 detections (PCRs)")) +
  #labs(size = paste("Predicted Aves\n", "species richness (MSOM)\n","[summed probabilities]")) +
  #labs(size = paste("Predicted Aves\n", "species richness (OccPlus)\n","[summed probabilities]")) +
  #labs(size = paste("Observed herptile\n", "species richness\n", "at ≥1/18 detections (PCRs)")) +
  labs(size = paste("Observed herptile\n", "species richness\n", "at ≥3/18 detections (PCRs)")) +
  #labs(size = paste("Predicted herptile\n", "species richness (MSOM)\n", "[summed probabilities]")) +
  #labs(size = paste("Predicted herptile\n", "species richness (OccPlus)\n", "[summed probabilities]")) +
  theme(legend.position = c(0.35, 0.65), 
        legend.title = element_text(face = "bold", hjust = 0, size = 12, lineheight = 1),
        legend.spacing = unit(0, 'cm'),
        legend.background = element_rect(fill = NA, color = NA)) +
  scale_alpha_manual(values = c("insideGLGS" = 1, "outsideGLGS" = 0.5)) +
  guides(alpha = guide_legend(override.aes = list(size = 4), keywidth = 1, keyheight = 1),
         size = guide_legend(order = 1)) +
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
  custom_scale_size +
  annotation_scale(location = "bl", 
                   width_hint = 0.2, 
                   style = "ticks", 
                   text_cex = 0.8, 
                   line_width = 0.5, 
                   bar_cols = c("black", "white")) +
  geom_tile(data = data.frame(x = 97.1, y = 25.7), aes(x, y), 
            fill = "green", width = 0.22, height = 0.14, alpha = 0.8) +
  annotate("text", x = 97.25, y = 25.7, label = "Protected Areas", hjust = 0, size = 4) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave(here("outputs", "maps1018", "GLGS_naive_TERR_richness_onlywild.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_TERR_richness_onlywild_cutoff3of18.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_TERR_richness1018_nodomestic.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_TERR_richness1018_nodomestic.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_FISH_richness.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_FISH_richness_cutoff3of18.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_FISH_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_FISH_richness1018_5factors.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_nativeFISH_richness.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_invasiveFISH_richness.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_nativeFISH_richness_cutoff3of18.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_invasiveFISH_richness_cutoff3of18.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_nativeFISH_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_nativeFISH_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_invasiveFISH_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_invasiveFISH_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_mammal_richness_onlywild.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_mammal_richness_onlywild_cutoff3of18.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_Mammal_richness1018_nodomestic.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_Mammal_richness1018_nodomestic.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_ave_richness_onlywild.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_ave_richness_onlywild_cutoff3of18.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_Ave_richness1018_nodomestic.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_Ave_richness1018_nodomestic.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_amphrept_richness_onlywild.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_naive_amphrept_richness_onlywild_cutoff3of18.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_AmphRept_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_AmphRept_richness1018.jpg"), width = 4, height = 8)


#######################################################################################################
# Plotting the maps of Species richness predictions for each order of mammals

# Preparing the required objects
GLGS <- st_read(here("data", "glgs", "GLG_Region_FULL_UTM.shp"))
PA <- st_read(here("data", "glgs", "GLGS_PAs.shp"))
raster_map <- raster(here("data", "glgs", "forest_1km.tif"))
raster_df <- as.data.frame(raster_map, xy = TRUE)
GLGS <- st_transform(GLGS, 4326)
PA <- st_transform(PA, 4326)

msom_occ <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Mammal_outputs_20241018",
                          "occupancy_bysite_msom_Mammal_20241018.csv"))
occPlus_occ <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Mammal_outputs_20241018_4factors",
                             "occupancybysite_occPlus_Mammal_20241018_4factors.csv"))
msom_occ <- msom_occ %>% 
  dplyr::select(!contains(c("Bos.gaurus","Bos.grunniens","Bos.taurus","Bubalus.bubalis","Capra.hircus","Ovis.aries","Cervus.nippon","Sus.scrofa","Canis.lupus","Felis.catus","Oryctolagus.cuniculu","Equus.asinus","Equus.caballus","Anas.platyrhynchos","Gallus.gallus","Melopsittacus.undulatus","Pelodiscus.sinensis","Trachemys.scripta","Rana.catesbeiana"))) 
occPlus_occ <- inner_join(glgs_sum2.1[,1:21], occPlus_occ, by = "site")
occPlus_occ <- occPlus_occ %>% 
  dplyr::select(!contains(c("Bos.gaurus","Bos.grunniens","Bos.taurus","Bubalus.bubalis","Capra.hircus","Ovis.aries","Cervus.nippon","Sus.scrofa","Canis.lupus","Felis.catus","Oryctolagus.cuniculu","Equus.asinus","Equus.caballus","Anas.platyrhynchos","Gallus.gallus","Melopsittacus.undulatus","Pelodiscus.sinensis","Trachemys.scripta","Rana.catesbeiana")))

spp_list <- spp_list %>% mutate(ORDER = str_extract(otuID, "(?<=_)[^_]+(?=_)"))
spp_list_Primates <- spp_list %>% 
  filter(ORDER == "Primates") %>% 
  mutate(spp = gsub(" ", ".", species)) %>% 
  mutate(spp = paste0("Mean.", spp))
spp_list_Rodentia <- spp_list %>% 
  filter(ORDER == "Rodentia") %>% 
  mutate(spp = gsub(" ", ".", species)) %>% 
  mutate(spp = paste0("Mean.", spp))
spp_list_Chiroptera <- spp_list %>% 
  filter(ORDER == "Chiroptera") %>% 
  mutate(spp = gsub(" ", ".", species)) %>% 
  mutate(spp = paste0("Mean.", spp))
spp_list_Carnivora <- spp_list %>% 
  filter(ORDER == "Carnivora") %>% 
  mutate(spp = gsub(" ", ".", species)) %>% 
  mutate(spp = paste0("Mean.", spp))
spp_list_Artiodactyla <- spp_list %>% 
  filter(ORDER == "Artiodactyla") %>% 
  mutate(spp = gsub(" ", ".", species)) %>% 
  mutate(spp = paste0("Mean.", spp))
spp_list_Eulipotyphla <- spp_list %>% 
  filter(ORDER == "Eulipotyphla") %>% 
  mutate(spp = gsub(" ", ".", species)) %>% 
  mutate(spp = paste0("Mean.", spp))

col_names <- colnames(msom_occ)
Primates <- which(col_names %in% spp_list_Primates$spp)
msom_occ$Primates_rich <- rowSums(msom_occ[Primates])
Rodentia <- which(col_names %in% spp_list_Rodentia$spp)
msom_occ$Rodentia_rich <- rowSums(msom_occ[Rodentia])
Chiroptera <- which(col_names %in% spp_list_Chiroptera$spp)
msom_occ$Chiroptera_rich <- rowSums(msom_occ[Chiroptera])
Carnivora <- which(col_names %in% spp_list_Carnivora$spp)
msom_occ$Carnivora_rich <- rowSums(msom_occ[Carnivora])
Artiodactyla <- which(col_names %in% spp_list_Artiodactyla$spp)
msom_occ$Artiodactyla_rich <- rowSums(msom_occ[Artiodactyla])
Eulipotyphla <- which(col_names %in% spp_list_Eulipotyphla$spp)
msom_occ$Eulipotyphla_rich <- rowSums(msom_occ[Eulipotyphla])
col_names <- colnames(occPlus_occ)
Primates <- which(col_names %in% spp_list_Primates$spp)
occPlus_occ$Primates_rich <- rowSums(occPlus_occ[Primates])
Rodentia <- which(col_names %in% spp_list_Rodentia$spp)
occPlus_occ$Rodentia_rich <- rowSums(occPlus_occ[Rodentia])
Chiroptera <- which(col_names %in% spp_list_Chiroptera$spp)
occPlus_occ$Chiroptera_rich <- rowSums(occPlus_occ[Chiroptera])
Carnivora <- which(col_names %in% spp_list_Carnivora$spp)
occPlus_occ$Carnivora_rich <- rowSums(occPlus_occ[Carnivora])
Artiodactyla <- which(col_names %in% spp_list_Artiodactyla$spp)
occPlus_occ$Artiodactyla_rich <- rowSums(occPlus_occ[Artiodactyla])
Eulipotyphla <- which(col_names %in% spp_list_Eulipotyphla$spp)
occPlus_occ$Eulipotyphla_rich <- rowSums(occPlus_occ[Eulipotyphla])

msom_occ_sf <- st_as_sf(msom_occ,coords = 6:7) %>% st_set_crs(4326)
occPlus_occ_sf <- st_as_sf(occPlus_occ,coords = 3:4) %>% st_set_crs(4326)

# for Primates
range(msom_occ_sf$Primates_rich) #[1]  0.05667595 1.72017203
range(occPlus_occ_sf$Primates_rich) #[1]  0.004582053 1.747433823
size_mapping <- function(value) {
  if (value <= 0.35) {
    return(1)
  } else if (value <= 0.7) {
    return(2)
  } else if (value <= 1.05) {
    return(3)
  } else if (value <= 1.4) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("0-0.35", "0.36-0.7", "0.71-1.05", "1.06-1.4", "1.41-1.75"),
)
msom_occ_sf$size_mapped <- sapply(msom_occ_sf$Primates_rich, size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf$Primates_rich, size_mapping)

# for Rodentia
range(msom_occ_sf$Rodentia_rich) #[1]  6.37969 16.51746
range(occPlus_occ_sf$Rodentia_rich) #[1]  0.3592355 13.9475996
size_mapping <- function(value) {
  if (value <= 3.3) {
    return(1)
  } else if (value <= 6.6) {
    return(2)
  } else if (value <= 9.9) {
    return(3)
  } else if (value <= 13.2) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("0-3.3", "3.4-6.6", "6.7-9.9", "10-13.2", "13.3-16.52"),
)
msom_occ_sf$size_mapped <- sapply(msom_occ_sf$Rodentia_rich, size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf$Rodentia_rich, size_mapping)

# for Chiroptera
range(msom_occ_sf$Chiroptera_rich) #[1]  0.4927859 8.1028818
range(occPlus_occ_sf$Chiroptera_rich) #[1]  0.02112693 4.34793740
size_mapping <- function(value) { # for amphrept
  if (value <= 1.62) {
    return(1)
  } else if (value <= 3.24) {
    return(2)
  } else if (value <= 4.86) {
    return(3)
  } else if (value <= 6.48) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("0-1.62", "1.63-3.24", "3.25-4.86", "4.87-6.48", "6.49-8.11"),
)
msom_occ_sf$size_mapped <- sapply(msom_occ_sf$Chiroptera_rich, size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf$Chiroptera_rich, size_mapping)

# for Carnivora
range(msom_occ_sf$Carnivora_rich) #[1]  1.745233 3.643870
range(occPlus_occ_sf$Carnivora_rich) #[1]  0.02304934 2.38214883
size_mapping <- function(value) { # for amphrept
  if (value <= 0.73) {
    return(1)
  } else if (value <= 1.46) {
    return(2)
  } else if (value <= 2.19) {
    return(3)
  } else if (value <= 2.92) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("0-0.73", "0.74-1.46", "1.47-2.19", "2.2-2.92", "2.93-3.65"),
)
msom_occ_sf$size_mapped <- sapply(msom_occ_sf$Carnivora_rich, size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf$Carnivora_rich, size_mapping)

# for Artiodactyla
range(msom_occ_sf$Artiodactyla_rich) #[1]  0.4338355 4.0581934
range(occPlus_occ_sf$Artiodactyla_rich) #[1]  0.01733508 3.78040043
size_mapping <- function(value) { # for amphrept
  if (value <= 0.81) {
    return(1)
  } else if (value <= 1.62) {
    return(2)
  } else if (value <= 2.43) {
    return(3)
  } else if (value <= 3.24) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("0-0.81", "0.82-1.62", "1.63-2.43", "2.44-3.24", "3.25-4.06"),
)
msom_occ_sf$size_mapped <- sapply(msom_occ_sf$Artiodactyla_rich, size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf$Artiodactyla_rich, size_mapping)

# for Eulipotyphla
range(msom_occ_sf$Eulipotyphla_rich) #[1]  2.202987 4.976948
range(occPlus_occ_sf$Eulipotyphla_rich) #[1]  0.2303644 4.3788727
size_mapping <- function(value) { # for amphrept
  if (value <= 1.18) {
    return(1)
  } else if (value <= 2.13) {
    return(2)
  } else if (value <= 3.08) {
    return(3)
  } else if (value <= 4.03) {
    return(4)
  } else {
    return(5)
  }
}
custom_scale_size <- scale_size_continuous(
  limits = c(1, 5), 
  breaks = 1:5, 
  labels = c("0.23-1.18", "1.19-2.13", "2.14-3.08", "3.09-4.03", "4.04-4.98"),
)
msom_occ_sf$size_mapped <- sapply(msom_occ_sf$Eulipotyphla_rich, size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf$Eulipotyphla_rich, size_mapping)

ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = forest_1km)) +
  scale_fill_gradient(low = "white", high = "grey70", na.value = NA) +
  guides(fill = "none") +
  geom_sf(data = GLGS, color = "darkblue", fill = NA, linetype = "solid", size = 8) +
  geom_sf(data = PA, color = "green", fill = "green", alpha = 0.8) +
  geom_sf(data = msom_occ_sf,
  #geom_sf(data = occPlus_occ_sf,
          aes(size = size_mapped, alpha = riverOrigin), 
          shape = 21, color = "black", fill = 'palevioletred4',
          show.legend = 'point', inherit.aes = FALSE) +
  #labs(size = paste("Predicted Primates\n", "species richness (MSOM)\n"," [summed probabilities]")) +
  #labs(size = paste("Predicted Primates\n", "species richness (OccPlus)\n","[summed probabilities]")) +
  #labs(size = paste("Predicted Rodentia\n", "species richness (MSOM)\n"," [summed probabilities]")) +
  #labs(size = paste("Predicted Rodentia\n", "species richness (OccPlus)\n","[summed probabilities]")) +
  #labs(size = paste("Predicted Chiroptera\n", "species richness (MSOM)\n"," [summed probabilities]")) +
  #labs(size = paste("Predicted Chiroptera\n", "species richness (OccPlus)\n","[summed probabilities]")) +
  #labs(size = paste("Predicted Carnivora\n", "species richness (MSOM)\n"," [summed probabilities]")) +
  #labs(size = paste("Predicted Carnivora\n", "species richness (OccPlus)\n","[summed probabilities]")) +
  #labs(size = paste("Predicted Artiodactyla\n", "species richness (MSOM)\n"," [summed probabilities]")) +
  #labs(size = paste("Predicted Artiodactyla\n", "species richness (OccPlus)\n","[summed probabilities]")) +
  labs(size = paste("Predicted Eulipotyphla\n", "species richness (MSOM)\n"," [summed probabilities]")) +
  #labs(size = paste("Predicted Eulipotyphla\n", "species richness (OccPlus)\n","[summed probabilities]")) +
  theme(legend.position = c(0.35, 0.65), 
        legend.title = element_text(face = "bold", hjust = 0, size = 12, lineheight = 1),
        legend.spacing = unit(0, 'cm'),
        legend.background = element_rect(fill = NA, color = NA)) +
  scale_alpha_manual(values = c("insideGLGS" = 1, "outsideGLGS" = 0.5)) +
  guides(alpha = guide_legend(override.aes = list(size = 4), keywidth = 1, keyheight = 1),
         size = guide_legend(order = 1)) +
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
  custom_scale_size +
  annotation_scale(location = "bl", 
                   width_hint = 0.2, 
                   style = "ticks", 
                   text_cex = 0.8, 
                   line_width = 0.5, 
                   bar_cols = c("black", "white")) +
  geom_tile(data = data.frame(x = 97.1, y = 25.7), aes(x, y), 
            fill = "green", width = 0.22, height = 0.14, alpha = 0.8) +
  annotate("text", x = 97.25, y = 25.7, label = "Protected Areas", hjust = 0, size = 4) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave(here("outputs", "maps1018", "GLGS_msom_Primates_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_Primates_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_Rodentia_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_Rodentia_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_Chiroptera_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_Chiroptera_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_Carnivora_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_Carnivora_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_Artiodactyla_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_Artiodactyla_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_msom_Eulipotyphla_richness1018.jpg"), width = 4, height = 8)
ggsave(here("outputs", "maps1018", "GLGS_occPlus_Eulipotyphla_richness1018.jpg"), width = 4, height = 8)


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

# Preparing spp lists for plotting maps of INDIVIDUALS SPECIES -----------------------------------------------

trial <- st_as_sf(glgs_sum2,coords = 3:4) %>% st_set_crs(4326)
msom_occ_sf <- st_as_sf(msom_occ,coords = 6:7) %>% st_set_crs(4326)
occPlus_occ_sf <- st_as_sf(occPlus_occ,coords = 3:4) %>% st_set_crs(4326)
GLGS <- st_read(here("data", "glgs", "GLG_Region_FULL_UTM.shp"))
PA <- st_read(here("data", "glgs", "GLGS_PAs.shp"))
GLGS <- st_transform(GLGS, 4326)
PA <- st_transform(PA, 4326)

classes_of_interest <- c("Amphibia", "Reptilia")
AmphRept_spp_oldnames <- data.frame(spp_name = spp_names) %>% 
  filter(str_detect(spp_name, str_c("OTU_(", str_c(classes_of_interest, collapse = "|"), ")_"))) %>% 
  arrange(spp_name)
aa <- str_split(AmphRept_spp_oldnames$spp_name, "_")
AmphRept_spp <- rep(NA, length(aa))
for(i in 1:length(aa)) {
  tmp <- aa[[i]][5:6]
  AmphRept_spp[i] <- paste(tmp[1],".",tmp[2], sep = "")
}

classes_of_interest <- c("Mammalia")
Mammal_spp_oldnames <- data.frame(spp_name = spp_names) %>% 
  filter(str_detect(spp_name, str_c("OTU_(", str_c(classes_of_interest, collapse = "|"), ")_"))) %>% 
  arrange(spp_name)
aa <- str_split(Mammal_spp_oldnames$spp_name, "_")
Mammal_spp <- rep(NA, length(aa))
for(i in 1:length(aa)) {
  tmp <- aa[[i]][5:6]
  Mammal_spp[i] <- paste(tmp[1],".",tmp[2], sep = "")
}

classes_of_interest <- c("Aves")
Ave_spp_oldnames <- data.frame(spp_name = spp_names) %>% 
  filter(str_detect(spp_name, str_c("OTU_(", str_c(classes_of_interest, collapse = "|"), ")_"))) %>% 
  arrange(spp_name)
aa <- str_split(Ave_spp_oldnames$spp_name, "_")
Ave_spp <- rep(NA, length(aa))
for(i in 1:length(aa)) {
  tmp <- aa[[i]][5:6]
  Ave_spp[i] <- paste(tmp[1],".",tmp[2], sep = "")
}

classes_of_interest <- c("Teleostei", "Actinopteri")
FISH_spp_oldnames <- data.frame(spp_name = spp_names) %>% 
  filter(str_detect(spp_name, str_c("OTU_(", str_c(classes_of_interest, collapse = "|"), ")_"))) %>% 
  arrange(spp_name)
aa <- str_split(FISH_spp_oldnames$spp_name, "_")
FISH_spp <- rep(NA, length(aa))
for(i in 1:length(aa)) {
  tmp <- aa[[i]][5:6]
  FISH_spp[i] <- paste(tmp[1],".",tmp[2], sep = "")
}

# Custom Size Mapping Functions

size_mapping <- function(value) {
  if (value <= 0.2) {
    return(1)
  } else if (value <= 0.4) {
    return(2)
  } else if (value <= 0.6) {
    return(3)
  } else if (value <= 0.8) {
    return(4)
  } else {
    return(5)
  }
}

# Plot the observed and predicted data from MSOM and occPlus for each sp on the same page-------------------------------------------------------------------------------------------------------
library(gridExtra)

pdf(here("outputs","maps1018","AmphRept_DETvsMSOMvsOCC_maps_for_eachSP1018_3factors.pdf"), width = 9, height = 6)
pdf(here("outputs","maps1018","Mammal_DETvsMSOMvsOCC_maps_for_eachSP1018_4factors.pdf"), width = 9, height = 6)
pdf(here("outputs","maps1018","Ave_DETvsMSOMvsOCC_maps_for_eachSP1018_2factors.pdf"), width = 9, height = 6)
pdf(here("outputs","maps1018","FISH_DETvsMSOMvsOCC_maps_for_eachSP1018_5factors.pdf"), width = 9, height = 6)

#for (spp in AmphRept_spp) {
#for (spp in Mammal_spp) {
#for (spp in Ave_spp) {
for (spp in FISH_spp) {
  
  spp_label <- gsub("\\.", " ", spp)
  trial_0 <- trial[trial[[spp]] == 0, ]
  trial_filtered <- trial[trial[[spp]] != 0, ]
  
  p1 <- ggplot() +
    geom_sf(data = GLGS, color = "darkblue", fill = "lightblue", alpha = 0.8, linetype = "solid", size = 1.5) + 
    geom_sf(data = PA, fill = "darkgreen", alpha = 0.8) + 
    
    geom_sf(data = trial_0, 
            shape = 21,
            fill = NA, color = "brown", size = 4, alpha = 0.7,
            inherit.aes = FALSE, show.legend = FALSE) +
    
    geom_sf(data = trial_filtered, 
            aes(fill = !!sym(spp)), shape = 21,
            color = "brown", size = 4, alpha = 0.7,
            inherit.aes = FALSE) +
    
    scale_shape_manual(values = c(1, 21)) +
    
    # geographic limits of map:
    coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
    scale_fill_gradient(low = "yellow", high = "red", 
                        name = paste("Detection of\n", spp_label), 
                        limits = c(1, 18), 
                        na.value = "transparent") +
    theme(legend.position = c(0.25, 0.72), 
          legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
          legend.background = element_rect(fill = NA, color = NA)) +
    
    geom_tile(data = data.frame(x = 97.25, y = 25.8), aes(x, y), 
              fill = "darkgreen", width = 0.22, height = 0.14, alpha = 0.8) +
    annotate("text", x = 97.4, y = 25.8, label = "Protected Areas", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.25, y = 26), aes(x, y), 
              fill = "lightblue", width = 0.22, height = 0.14, alpha = 0.6, color = "blue") +
    annotate("text", x = 97.4, y = 26, label = "Gaoligong Region", hjust = 0, size = 3) +
    # Adding a scale
    annotation_scale(location = "bl", 
                     width_hint = 0.2, 
                     style = "ticks", 
                     text_cex = 0.8, 
                     line_width = 0.5, 
                     bar_cols = c("black", "white")) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 6)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  msom_occ_sf$size_mapped <- sapply(msom_occ_sf[[paste0("Mean.", spp)]], size_mapping)

  p2 <- ggplot() +
    geom_sf(data = GLGS, color = "darkblue", fill = "lightblue", alpha = 0.8, linetype = "solid", size = 1.5) + 
    geom_sf(data = PA, fill = "darkgreen", alpha = 0.8) + 
    geom_sf(data = msom_occ_sf, 
            aes(size = size_mapped), 
            shape = 21, color = "black", fill = 'purple', alpha = 0.6,
            show.legend = 'point', inherit.aes = FALSE) +
    # geographic limits of map:
    coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
    labs(size = paste("MSOM\n", "Occupancy probablity\n", spp_label), alpha = 0.6) +
    theme(legend.position = c(0.32, 0.68), 
          legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
          legend.background = element_rect(fill = NA, color = NA)) +
    geom_tile(data = data.frame(x = 97.25, y = 25.7), aes(x, y), 
              fill = "darkgreen", width = 0.22, height = 0.14, alpha = 0.8) +
    annotate("text", x = 97.4, y = 25.7, label = "Protected Areas", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.25, y = 25.9), aes(x, y), 
              fill = "lightblue", width = 0.22, height = 0.14, alpha = 0.6, color = "blue") +
    annotate("text", x = 97.4, y = 25.9, label = "Gaoligong Region", hjust = 0, size = 3) +
    scale_size_continuous(limits = c(1, 5), breaks = 1:5, 
                          #labels = c("0-0.1", "0.1-0.3", "0.3-0.5", "0.5-0.7", "0.7-1")) +
                          labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
    # Adding a scale
    annotation_scale(location = "bl", 
                     width_hint = 0.2, 
                     style = "ticks", 
                     text_cex = 0.8, 
                     line_width = 0.5, 
                     bar_cols = c("black", "white")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf[[paste0("Mean.", spp)]], size_mapping)
  
  p3 <- ggplot() +
    geom_sf(data = GLGS, color = "darkblue", fill = "lightblue", alpha = 0.8, linetype = "solid", size = 1.5) + 
    geom_sf(data = PA, fill = "darkgreen", alpha = 0.8) + 
    geom_sf(data = occPlus_occ_sf, 
            aes(size = size_mapped), 
            shape = 21, color = "darkred", fill = 'orange', alpha = 0.6,
            show.legend = 'point', inherit.aes = FALSE) +
    # geographic limits of map:
    coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
    labs(size = paste("OccPlus\n", "Occupancy probablity\n", spp_label), alpha = 0.6) +
    theme(legend.position = c(0.32, 0.68), 
          legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
          legend.background = element_rect(fill = NA, color = NA)) +
    geom_tile(data = data.frame(x = 97.25, y = 25.7), aes(x, y), 
              fill = "darkgreen", width = 0.22, height = 0.14, alpha = 0.8) +
    annotate("text", x = 97.4, y = 25.7, label = "Protected Areas", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.25, y = 25.9), aes(x, y), 
              fill = "lightblue", width = 0.22, height = 0.14, alpha = 0.6, color = "blue") +
    annotate("text", x = 97.4, y = 25.9, label = "Gaoligong Region", hjust = 0, size = 3) +
    scale_size_continuous(limits = c(1, 5), breaks = 1:5, 
                          #labels = c("0-0.1", "0.1-0.3", "0.3-0.5", "0.5-0.7", "0.7-1")) +
                          labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
  # Adding a scale
    annotation_scale(location = "bl", 
                     width_hint = 0.2, 
                     style = "ticks", 
                     text_cex = 0.8, 
                     line_width = 0.5, 
                     bar_cols = c("black", "white")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  grid.arrange(p1, p3, p2, ncol = 3)
}
dev.off()


# Plot the observed and predicted data from MSOM and occPlus only for Budorcas taxicolor-------------------------------------------------------------------------------------------------------
GLGS <- st_read(here("data", "glgs", "GLG_Region_FULL_UTM.shp"))
GLG_PA <- st_read(here("data", "glgs", "GLG_PA.shp"))
TBG_PA <- st_read(here("data", "glgs", "TBG_PA.shp"))
XHS_PA <- st_read(here("data", "glgs", "XHS_PA.shp"))
GLGS <- st_transform(GLGS, 4326)
GLG_PA <- st_transform(GLG_PA, 4326)
TBG_PA <- st_transform(TBG_PA, 4326)
XHS_PA <- st_transform(XHS_PA, 4326)
GLGS_east <- st_read(here("data", "glgs", "hybas05_GLG_East.shp"))
GLGS_SW <- st_read(here("data", "glgs", "hybas05_GLG_SW.shp"))
GLGS_NW <- st_read(here("data", "glgs", "hybas05_GLG_NW.shp"))
GLGS_dy <- st_read(here("data", "glgs", "hybas05_GLG_DY.shp"))
GLGS_lc <- st_read(here("data", "glgs", "hybas05_GLG_LC.shp"))
GLGS_east <- st_transform(GLGS_east, 4326)
GLGS_SW <- st_transform(GLGS_SW, 4326)
GLGS_NW <- st_transform(GLGS_NW, 4326)
GLGS_dy <- st_transform(GLGS_dy, 4326)
GLGS_lc <- st_transform(GLGS_lc, 4326)

  trial_0 <- trial[trial[["Budorcas.taxicolor"]] == 0, ]
  trial_filtered <- trial[trial[["Budorcas.taxicolor"]] != 0, ]
  
  p1 <- ggplot() +
    geom_sf(data = GLGS, color = "darkblue", fill = "lightblue", alpha = 0.8, linetype = "solid", size = 1.5) + 
    geom_sf(data = GLG_PA, color = "darkgreen", fill = "darkgreen", alpha = 0.8) + 
    geom_sf(data = TBG_PA, color = "deeppink", fill = "deeppink", alpha = 1) + 
    geom_sf(data = XHS_PA, color = "dodgerblue", fill = "dodgerblue", alpha = 1) + 
    
    geom_sf(data = trial_0, 
            shape = 21,
            fill = NA, color = "brown", size = 4, alpha = 0.7,
            inherit.aes = FALSE, show.legend = FALSE) +
    
    geom_sf(data = trial_filtered, 
            aes(fill = Budorcas.taxicolor), shape = 21,
            color = "brown", size = 4, alpha = 0.7,
            inherit.aes = FALSE) +
    
    scale_shape_manual(values = c(1, 21)) +
    
    # geographic limits of map:
    coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
    scale_fill_gradient(low = "yellow", high = "red", 
                        name = paste("Detection of\n", "Budorcas taxicolor"), 
                        limits = c(1, 18), 
                        na.value = "transparent") +
    theme(legend.position = c(0.25, 0.72), 
          legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
          legend.background = element_rect(fill = NA, color = NA)) +
    
    geom_tile(data = data.frame(x = 97.25, y = 25.7), aes(x, y), 
              fill = "dodgerblue", width = 0.22, height = 0.14, alpha = 0.8) +
    annotate("text", x = 97.4, y = 25.7, label = "Xiaoheishan PA", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.25, y = 25.9), aes(x, y), 
              fill = "deeppink", width = 0.22, height = 0.14, alpha = 0.8) +
    annotate("text", x = 97.4, y = 25.9, label = "Tongbiguan PA", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.25, y = 26.1), aes(x, y), 
              fill = "darkgreen", width = 0.22, height = 0.14, alpha = 0.8) +
    annotate("text", x = 97.4, y = 26.1, label = "Gaoligongshan PA", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.25, y = 26.3), aes(x, y), 
              fill = "lightblue", width = 0.22, height = 0.14, alpha = 0.6, color = "blue") +
    annotate("text", x = 97.4, y = 26.3, label = "Gaoligong Region", hjust = 0, size = 3) +
    # Adding a scale
    annotation_scale(location = "bl", 
                     width_hint = 0.2, 
                     style = "ticks", 
                     text_cex = 0.8, 
                     line_width = 0.5, 
                     bar_cols = c("black", "white")) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 6)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  msom_occ_sf$size_mapped <- sapply(msom_occ_sf[["Mean.Budorcas.taxicolor"]], size_mapping)
  
  p2 <- ggplot() +
    geom_sf(data = GLGS_east, color = "dodgerblue4", fill = "palegreen", alpha = 0.8, linetype = "solid", size = 1.5) + 
    geom_sf(data = GLGS_NW, color = "dodgerblue4", fill = "paleturquoise", alpha = 0.8, linetype = "solid", size = 1.5) + 
    geom_sf(data = GLGS_dy, color = "dodgerblue4", fill = "rosybrown", alpha = 0.8, linetype = "solid", size = 1.5) + 
    geom_sf(data = GLGS_lc, color = "dodgerblue4", fill = "moccasin", alpha = 0.8, linetype = "solid", size = 1.5) + 
    geom_sf(data = msom_occ_sf, 
            aes(size = size_mapped), 
            shape = 21, color = "black", fill = 'purple', alpha = 0.6,
            show.legend = 'point', inherit.aes = FALSE) +
    # geographic limits of map:
    coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
    labs(size = paste("MSOM\n", "Occupancy probablity"), alpha = 0.6) +
    theme(legend.position = c(0.32, 0.72), 
          legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
          legend.background = element_rect(fill = NA, color = NA)) +
    geom_tile(data = data.frame(x = 97.15, y = 25.7), aes(x, y), 
              fill = "rosybrown", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 25.7, label = "Daying River basin", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.15, y = 25.9), aes(x, y), 
              fill = "moccasin", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 25.9, label = "Longchuan River basin", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.15, y = 26.1), aes(x, y), 
              fill = "palegreen", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 26.1, label = "Nu River basin", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.15, y = 26.3), aes(x, y), 
              fill = "paleturquoise", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 26.3, label = "Dulong River basin", hjust = 0, size = 3) +
    scale_size_continuous(limits = c(1, 5), breaks = 1:5, 
                          labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
    #guides(size = guide_legend(ncol = 2)) +
    # Adding a scale
    annotation_scale(location = "bl", 
                     width_hint = 0.2, 
                     style = "ticks", 
                     text_cex = 0.8, 
                     line_width = 0.5, 
                     bar_cols = c("black", "white")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf[["Mean.Budorcas.taxicolor"]], size_mapping)
  
  p3 <- ggplot() +
    geom_sf(data = GLGS_east, color = "dodgerblue4", fill = "palegreen", alpha = 0.8, linetype = "solid", size = 1.5) + 
    geom_sf(data = GLGS_NW, color = "dodgerblue4", fill = "paleturquoise", alpha = 0.8, linetype = "solid", size = 1.5) + 
    geom_sf(data = GLGS_SW, color = "dodgerblue4", fill = "peachpuff", alpha = 0.8, linetype = "solid", size = 1.5) + 
    geom_sf(data = occPlus_occ_sf, 
            aes(size = size_mapped), 
            shape = 21, color = "darkred", fill = 'orange', alpha = 0.6,
            show.legend = 'point', inherit.aes = FALSE) +
    # geographic limits of map:
    coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
    labs(size = paste("OccPlus\n", "Occupancy probablity"), alpha = 0.6) +
    theme(legend.position = c(0.32, 0.7), 
          legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
          legend.background = element_rect(fill = NA, color = NA)) +
    geom_tile(data = data.frame(x = 97.15, y = 25.7), aes(x, y), 
              fill = "peachpuff", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 25.7, label = "Southwest Subregion", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.15, y = 25.9), aes(x, y), 
              fill = "palegreen", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 25.9, label = "East Subregion", hjust = 0, size = 3) +
    geom_tile(data = data.frame(x = 97.15, y = 26.1), aes(x, y), 
              fill = "paleturquoise", width = 0.22, height = 0.14, alpha = 0.8, color = "dodgerblue4") +
    annotate("text", x = 97.3, y = 26.1, label = "Northwest Subregion", hjust = 0, size = 3) +
    scale_size_continuous(limits = c(1, 5), breaks = 1:5, 
                          #labels = c("0-0.1", "0.1-0.3", "0.3-0.5", "0.5-0.7", "0.7-1")) +
                          labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
    # Adding a scale
    annotation_scale(location = "bl", 
                     width_hint = 0.2, 
                     style = "ticks", 
                     text_cex = 0.8, 
                     line_width = 0.5, 
                     bar_cols = c("black", "white")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  pdf(here("outputs","maps1018","Budorcas_taxicolor_DETvsMSOMvsOCC_maps.pdf"), width = 9, height = 6)
  
  grid.arrange(p1, p3, p2, ncol = 3)

dev.off()


# Plot the observed and predicted data for Rana catesbeiana -------------------------------------------------------------------------------------------------------
GLGS <- st_read(here("data", "glgs", "GLG_Region_FULL_UTM.shp"))
GLGS <- st_transform(GLGS, 4326)
PA <- st_read(here("data", "glgs", "GLGS_PAs.shp"))
PA <- st_transform(PA, 4326)
glgs_tmp <- glgs %>% 
  mutate(across(starts_with("Amphibia"), binarise1)) %>% 
  summarise(
    across(Filter_vol:DistFromPAedge, \(x) first(x)),
    across(starts_with("Amphibia"), \(x) sum(x)), 
    .by = c(site, Rep_of_sample)
  )  %>% 
  filter(status == "sample") %>%  # remove negative controls
  summarise(
    across(Filter_vol:DistFromPAedge, \(x) first(x)),
    across(starts_with("Amphibia"), \(x) sum(x)), 
    .by = c(site)
  ) 
trial <- st_as_sf(glgs_tmp,coords = 3:4) %>% st_set_crs(4326)
msom_occ <- read.csv(here("outputs", "output_bullfrog", "occupancy_preds_MSOM_Bullfrogs.csv"))
msom_occ <- msom_occ %>% 
  left_join(glgs_sum2[,1:4],by = c("obs.covs1.site"="site"))
occPlus_occ <- read.csv(here("outputs","output_bullfrog","occupancybysite_occPlus_bullfrog_20250217.csv"))
occPlus_occ <- occPlus_occ %>% 
  left_join(glgs_sum2[,1:4],by = "site")
msom_occ_sf <- st_as_sf(msom_occ,coords = 8:9) %>% st_set_crs(4326)
occPlus_occ_sf <- st_as_sf(occPlus_occ,coords = 9:10) %>% st_set_crs(4326)
size_mapping <- function(value) {
  if (value <= 0.2) {
    return(1)
  } else if (value <= 0.4) {
    return(2)
  } else if (value <= 0.6) {
    return(3)
  } else if (value <= 0.8) {
    return(4)
  } else {
    return(5)
  }
}
msom_occ_sf$size_mapped <- sapply(msom_occ_sf[["psi.mean"]], size_mapping)
occPlus_occ_sf$size_mapped <- sapply(occPlus_occ_sf[["Mean.Rana.catesbeiana"]], size_mapping)
trial_0 <- trial[trial[["Amphibia_Anura_Ranidae_Rana_catesbeiana_OTU414_31973"]] == 0, ]
trial_filtered <- trial[trial[["Amphibia_Anura_Ranidae_Rana_catesbeiana_OTU414_31973"]] != 0, ]

p1 <- ggplot() +
  geom_sf(data = GLGS, color = "darkblue", fill = "lightblue", alpha = 0.8, linetype = "solid", size = 1.5) + 
  geom_sf(data = PA, fill = "darkgreen", alpha = 0.8) + 
  
  geom_sf(data = trial_0, 
          shape = 21,
          fill = NA, color = "brown", size = 4, alpha = 0.7,
          inherit.aes = FALSE, show.legend = FALSE) +
  
  geom_sf(data = trial_filtered, 
          aes(fill = Amphibia_Anura_Ranidae_Rana_catesbeiana_OTU414_31973), shape = 21,
          color = "brown", size = 4, alpha = 0.7,
          inherit.aes = FALSE) +
  
  scale_shape_manual(values = c(1, 21)) +
  
  # geographic limits of map:
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
  scale_fill_gradient(low = "yellow", high = "red", 
                      name = paste("Detection of\n", "Rana catesbeiana"), 
                      limits = c(1, 18), 
                      na.value = "transparent") +
  theme(legend.position = c(0.25, 0.72), 
        legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
        legend.background = element_rect(fill = NA, color = NA)) +
  
  geom_tile(data = data.frame(x = 97.25, y = 25.8), aes(x, y), 
            fill = "darkgreen", width = 0.22, height = 0.14, alpha = 0.8) +
  annotate("text", x = 97.4, y = 25.8, label = "Protected Areas", hjust = 0, size = 3) +
  geom_tile(data = data.frame(x = 97.25, y = 26), aes(x, y), 
            fill = "lightblue", width = 0.22, height = 0.14, alpha = 0.6, color = "blue") +
  annotate("text", x = 97.4, y = 26, label = "Gaoligong Region", hjust = 0, size = 3) +
  # Adding a scale
  annotation_scale(location = "bl", 
                   width_hint = 0.2, 
                   style = "ticks", 
                   text_cex = 0.8, 
                   line_width = 0.5, 
                   bar_cols = c("black", "white")) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 6)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p2 <- ggplot() +
  geom_sf(data = GLGS, color = "darkblue", fill = "lightblue", alpha = 0.8, linetype = "solid", size = 1.5) + 
  geom_sf(data = PA, fill = "darkgreen", alpha = 0.8) + 
  geom_sf(data = msom_occ_sf, 
          aes(size = size_mapped), 
          shape = 21, color = "black", fill = 'purple', alpha = 0.6,
          show.legend = 'point', inherit.aes = FALSE) +
  # geographic limits of map:
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
  labs(size = paste("MSOM\n", "Single-Species Occupancy Model\n", "Occupancy probablity\n", "Rana catesbeiana"), alpha = 0.6) +
  theme(legend.position = c(0.32, 0.68), 
        legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
        legend.background = element_rect(fill = NA, color = NA)) +
  geom_tile(data = data.frame(x = 97.25, y = 25.7), aes(x, y), 
            fill = "darkgreen", width = 0.22, height = 0.14, alpha = 0.8) +
  annotate("text", x = 97.4, y = 25.7, label = "Protected Areas", hjust = 0, size = 3) +
  geom_tile(data = data.frame(x = 97.25, y = 25.9), aes(x, y), 
            fill = "lightblue", width = 0.22, height = 0.14, alpha = 0.6, color = "blue") +
  annotate("text", x = 97.4, y = 25.9, label = "Gaoligong Region", hjust = 0, size = 3) +
  scale_size_continuous(limits = c(1, 5), breaks = 1:5, 
                        #labels = c("0-0.1", "0.1-0.3", "0.3-0.5", "0.5-0.7", "0.7-1")) +
                        labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
  # Adding a scale
  annotation_scale(location = "bl", 
                   width_hint = 0.2, 
                   style = "ticks", 
                   text_cex = 0.8, 
                   line_width = 0.5, 
                   bar_cols = c("black", "white")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p3 <- ggplot() +
  geom_sf(data = GLGS, color = "darkblue", fill = "lightblue", alpha = 0.8, linetype = "solid", size = 1.5) + 
  geom_sf(data = PA, fill = "darkgreen", alpha = 0.8) + 
  geom_sf(data = occPlus_occ_sf, 
          aes(size = size_mapped), 
          shape = 21, color = "darkred", fill = 'orange', alpha = 0.6,
          show.legend = 'point', inherit.aes = FALSE) +
  # geographic limits of map:
  coord_sf(xlim = c(97, 99.3), ylim = c(23.6, 28.4), expand = FALSE) +
  labs(size = paste("OccPlus\n", "Occupancy probablity\n", "Rana catesbeiana"), alpha = 0.6) +
  theme(legend.position = c(0.32, 0.68), 
        legend.title = element_text(face = "bold", hjust = 0.5, size = 10, lineheight = 1.2),
        legend.background = element_rect(fill = NA, color = NA)) +
  geom_tile(data = data.frame(x = 97.25, y = 25.7), aes(x, y), 
            fill = "darkgreen", width = 0.22, height = 0.14, alpha = 0.8) +
  annotate("text", x = 97.4, y = 25.7, label = "Protected Areas", hjust = 0, size = 3) +
  geom_tile(data = data.frame(x = 97.25, y = 25.9), aes(x, y), 
            fill = "lightblue", width = 0.22, height = 0.14, alpha = 0.6, color = "blue") +
  annotate("text", x = 97.4, y = 25.9, label = "Gaoligong Region", hjust = 0, size = 3) +
  scale_size_continuous(limits = c(1, 5), breaks = 1:5, 
                        #labels = c("0-0.1", "0.1-0.3", "0.3-0.5", "0.5-0.7", "0.7-1")) +
                        labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")) +
  # Adding a scale
  annotation_scale(location = "bl", 
                   width_hint = 0.2, 
                   style = "ticks", 
                   text_cex = 0.8, 
                   line_width = 0.5, 
                   bar_cols = c("black", "white")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())


library(gridExtra)
pdf(here("outputs","output_bullfrog","Bullfrog_DETvsSSPOCC_maps20250217.pdf"), width = 9, height = 6)
grid.arrange(p1, p3, p2, ncol = 3)
dev.off()



#######################################################################################################
# Correlation plot of species richness among different groups

# observed
observed_speciesrichness <- glgs_sum2.1 %>% dplyr::select(amph_rich, mammal_rich, ave_rich, fish_rich, fish_native_rich, fish_nonnative_rich)
library(corrplot)
cor_matrix <- cor(observed_speciesrichness, use = "complete.obs")
corrplot(cor_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45)

col_names <- colnames(occPlus_occ)
occPlus_occ$Primates_rich <- rowSums(occPlus_occ[Primates])
occPlus_occ$Rodentia_rich <- rowSums(occPlus_occ[Rodentia])
occPlus_occ$Chiroptera_rich <- rowSums(occPlus_occ[Chiroptera])
occPlus_occ$Carnivora_rich <- rowSums(occPlus_occ[Carnivora])
occPlus_occ$Artiodactyla_rich <- rowSums(occPlus_occ[Artiodactyla])
occPlus_occ$Eulipotyphla_rich <- rowSums(occPlus_occ[Eulipotyphla])


# occPlus
occPlus_orderrich_mam <- occPlus_occ %>% dplyr::select(Primates_rich:Eulipotyphla_rich)
cor_matrix <- cor(occPlus_orderrich_mam, use = "complete.obs")
corrplot(cor_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45)
model <- lm(Rodentia_rich ~ Carnivora_rich, data = occPlus_orderrich_mam)
model <- lm(Eulipotyphla_rich ~ Carnivora_rich, data = occPlus_orderrich_mam)
summary(model)

# estimate Pearson r correltion coefficient
cor_test <- cor.test(occPlus_orderrich_mam$Carnivora_rich, occPlus_orderrich_mam$Rodentia_rich, method = "pearson")
cor_test$estimate
cor_test$p.value
cor_test <- cor.test(occPlus_orderrich_mam$Carnivora_rich, occPlus_orderrich_mam$Chiroptera_rich, method = "pearson")
cor_test <- cor.test(occPlus_orderrich_mam$Carnivora_rich, occPlus_orderrich_mam$Eulipotyphla_rich, method = "pearson")


