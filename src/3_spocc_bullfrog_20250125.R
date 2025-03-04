
# script-specific packages
suppressPackageStartupMessages({
  # add script-specific packages here
  library(mapview)
  library(sf)
})

# general-use packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(fs)
  library(glue)
  library(readxl)
  library(cowplot)
  library(patchwork)
  library(lubridate)
  library(broom)
  library(ggeffects)
  library(viridis)
  library(arsenal) # for tableby()
  library(waldo) # for compare()
  library(sjmisc) # for rotate_df()
  library(envDocument)
  library(conflicted)
  library(knitr)
  library(beepr)
  library(pivottabler)
  library(furrr)
  library(scales)
  library(janitor)
  library(tictoc)
  library(inspectdf) # summary figures and tables of dataframes
  library(visdat)
  library(equatiomatic) 
  library(report) # generate text summaries of model outputs
  library(RColorBrewer)
  library(ggrepel)
  library(spOccupancy)
  library(coda)
  library(stars)
  library(ggplot2)
  library(cowplot)
  library(ggspatial)
  library(shapefiles)
  library(gstat)
  library(sp)
  library(terra)
})

# Sometimes, two or more packages use the same function names. The {conflicted} package lets you set which package gets precedence. For example, the next line enforces that filter() refers to the {dplyr} package. If you want to use the command filter() from a different package, you just need to precede it with the desired package name like this: stats::filter.
conflict_prefer("mutate", "dplyr", quiet = TRUE)
conflict_prefer("select", "dplyr", quiet = TRUE)
conflict_prefer("summarise", "dplyr", quiet = TRUE)
conflict_prefer("filter", "dplyr", quiet = TRUE)
conflict_prefer("first", "dplyr", quiet = TRUE)
conflict_prefer("here", "here", quiet = TRUE)
conflict_prefer("separate", "tidyr", quiet = TRUE)
conflict_prefer("unite", "tidyr", quiet = TRUE)
conflict_prefer("trim", "sjmisc", quiet=TRUE)
conflict_prefer("rescale", "scales", quiet=TRUE)
conflict_prefer("intersect", "dplyr", quiet = TRUE)
conflict_prefer("setdiff", "dplyr", quiet = TRUE) # w/out this, R crashes
conflict_prefer("to_factor", "sjmisc", quiet = TRUE)
conflict_prefer("trim", "glue", quiet = TRUE)
conflict_prefer("discard", "purrr", quiet = TRUE)
conflict_prefer("extract", "tidyr", quiet = TRUE)
conflict_prefer("na_if", "dplyr", quiet = TRUE)
conflict_prefer("expand", "tidyr", quiet = TRUE)
conflict_prefer("lmer", "lme4", quiet = TRUE)
conflict_prefer("col_factor", "scales", quiet = TRUE)

# Print real numbers, not scientific notation.
options(scipen = 999)

# R version
R.version.string

# This command reports the R environment that you are using, especially the package and R versions. 
# It is a good idea to include this information at the end of every report you generate, 
# to document the details of your computing environment. 
envDocument::env_doc("table", git = FALSE)
# info <- envDocument::env_doc("return", git = FALSE)
# sessionInfo() # base R method
# sessioninfo::session_info() # package sessioninfo method

# Print real numbers, not scientific notation.
options(scipen = 999)



# READ DATA ---------------------------------------------------------------


glgs <- read_tsv(here("data","OTUtable_12S_toSP_GLG23_20240528.txt"))

# add covariate:  riverOrigin (insideGLGS, outsideGLGS). outsideGLGS for all NA values of DistFromPAedge
glgs <- glgs |> 
  mutate(
    riverOrigin = case_when(
      is.na(DistFromPAedge) ~ "outsideGLGS",
      TRUE ~ "insideGLGS"
    )
  ) |> 
  relocate(riverOrigin, .before = DistFromPAedge)

new_covs <- read_tsv(here("data","23GLG_covs_info_20240927.txt"))

zones <- read_tsv(here("data","23GLG_zones_info_20240903.txt"))

new_covs <- left_join(new_covs, zones, by = "site")

binarise1 <- function(x, na.rm = FALSE) vegan::decostand(x, method = "pa") 

climate_covs <- read_csv(here("data","GLGS_climate.csv"))

new_covs <- left_join(new_covs, climate_covs, by = "site") %>% select(site:zone, annual_precipitation, mean_temperature)


# IDENTIFY SPECIES DETECTED AT A SINGLE SITE --------------------------------

# Amphibia and Reptilia
spptoremove <- glgs %>% 
  as.data.frame() %>%
  #dplyr::select(!contains(c("Spikein","Actinopteri", "Aves", "Amphibia", "Reptilia", "Teleostei"))) %>% # extract just a taxon data
  dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Aves", "Mammalia"))) %>% # remove fish
  dplyr::select(SampleID, PCR, site, contains("OTU")) %>% 
  rename_with(.cols = starts_with(c("Amphibia", "Reptilia")), 
              function(x){paste0("OTU_", x)}) %>% 
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


# SUMMARIZE DATA AND PREPARE FOR MSOM -------------------------------------

# binarise read number and sum by site and sample replicate. For instance, if a species is detected in all 3 sample replicates and all 6 PCRs per sample replicate, the output will be 6,6,6 for that species and site.
glgs_sum1 <- glgs %>% 
  
  # remove fish, Aves and Mammalia
  select(!starts_with(c("Spikein", "Actinopteri", "Teleostei", "Aves", "Mammalia"))) |>
  
  # add "OTU_" to the beginning of all OTU columns
  rename_with(.cols = starts_with(c("Amphibia", "Reptilia")), 
              function(x){paste0("OTU_", x)}) |>  
  # remove bullfrog
  #dplyr::select(!contains(c("Rana_catesbeiana"))) |>
  # retain bullfrog
  #dplyr::select(contains(c("Rana_catesbeiana"))) |>
  
  # change numbers in OTU columns to 0/1
  mutate(across(starts_with("OTU_"), binarise1)) |> 
  
  # sum detections by site and Rep_of_sample (stage 1 replicates)
  summarise(
    across(Filter_vol:DistFromPAedge, \(x) first(x)),
    across(starts_with("OTU_"), \(x) sum(x)), 
    .by = c(site, Rep_of_sample)
  ) |>
  filter(status == "sample") # remove negative controls


# sum by site, so if a species is detected in all 3 sample replicates and all 
#6 PCRs per sample replicate, the output will be 18 for that species and site. 
glgs_sum2 <- glgs_sum1 |> 
  # sum detections within site
  summarise(
    across(Filter_vol:DistFromPAedge, \(x) first(x)),
    across(starts_with("OTU_"), \(x) sum(x)), 
    .by = c(site)
  ) 
str(glgs_sum2)


# prepare data for evaluating naive occupancy
glgs_sum2.1 <- glgs_sum2 |>
  mutate(across(starts_with("OTU_"), binarise1))

str(glgs_sum2.1)

# naive detections
glgs_sum3 <- glgs_sum2 |> 
  # binarise
  mutate(across(starts_with("OTU_"), binarise1)) |> 
  dplyr::select(-one_of(spptoremove)) %>%
  # sum all detections
  summarise(
    across(starts_with("OTU_"), \(x) sum(x)) 
  ) |> 
  # rotate
  sjmisc::rotate_df(rn = "OTU") |> 
  rename(site_detections = V1)


# COMPILE OBJECT FOR MSOM ANALYSIS ----------------------------------------

# create list for site covariates
site.covs <- glgs_sum2 |>
  select(site, altitude, river, tributary, longitude, latitude, riverOrigin, river_width, TEMP)

# add land cover vars
site.covs1 <- inner_join(site.covs, new_covs, by = "site") 

# tributary is character; it needs to be numeric for random effect
site.covs1$tributary <- as.factor(site.covs1$tributary)
site.covs1$tributary <- as.numeric(site.covs1$tributary)

# river, riverOrigin and Regional_Division are characters; it needs to be factor for model
site.covs1$river <- as.factor(site.covs1$river)
site.covs1$riverOrigin <- as.factor(site.covs1$riverOrigin)
site.covs1$Regional_Division <- as.factor(site.covs1$Regional_Division)

# create list for observation covariates
obs.covs <- glgs_sum2 |>
  select(site, Filter_vol, team, PrevDayRain, water_velocity, TEMP)

obs.covs1 <- with(obs.covs, {
  list(site = site, Filter_vol = Filter_vol, team = team, TEMP = TEMP, 
       PrevDayRain = PrevDayRain, water_velocity = water_velocity)})

# create a 3D array for species data (3 occasions)
species_for_y <- glgs_sum1 |>
  select(starts_with("OTU_")) |>
  mutate(across(starts_with("OTU_"), binarise1)) |>
  dplyr::select(-one_of(spptoremove)) %>% # remove species with a single site detection
  dplyr::select(contains(c("Rana_catesbeiana"))) 
data.frame()

str(species_for_y)

species_for_y.1 <- matrix(as.matrix(species_for_y), nrow = 101, ncol = 3, byrow = T)

spp_names <- colnames(species_for_y)
str(spp_names)
spp_names

# combine observations, occupancy covariates and detection covariates
# coordinates can be provided for spatial model 
bullfrog <- list(y = species_for_y.1, 
                  occ.covs = site.covs1,
                  det.covs = obs.covs1)
#coords = coords)
str(bullfrog)


# SET UP AND RUN THE MODEL ------------------------------------------------

# occupancy formula
occ.formula <- ~ scale(annual_precipitation) + scale(mean_temperature) + 
  scale(Percent_of_forest_cover_1km) + scale(Human_modification_1km) + 
  scale(Patch_density_1km) + scale(Edge_density_1km) + 
  riverOrigin + Regional_Division + river +
  (1|tributary)
# delete TEMP, because it is strongly correlated with human_modification and mean_temperature
# keep mean_temperature and delete altitude AND latitude

# detection formula
det.formula <- ~ scale(Filter_vol) + scale(PrevDayRain) + (1|team)

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))
# Initial values
inits.list <- list(alpha = 0, beta = 0,
                   z = apply(bullfrog$y, 1, max, na.rm = TRUE))


# RUN MODEL!
out <- PGOcc(occ.formula = occ.formula, 
                  det.formula = det.formula, 
                  data = bullfrog, 
                  inits = inits.list, 
                  n.samples = 50000, 
                  priors = prior.list, 
                  n.omp.threads = 3, 
                  verbose = TRUE, 
                  n.report = 2000, 
                  n.burn = 20000,
                  n.thin = 20, 
                  n.chains = 3)


summary(out)

# predictive check
ppc.out <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out)

# traceplots to examine convergence
#pdf(here("outputs","msom_outputs_Viorel covs","AmphRept_outputs_20241018","msom_AmphRept_beta_traceplot_20241018.pdf"))
plot(out, 'beta', density = TRUE) # Occupancy parameters.
#dev.off()
#pdf(here("outputs","msom_outputs_Viorel covs","AmphRept_outputs_20241018","msom_AmphRept_alpha_traceplot_20241018.pdf"))
plot(out, 'alpha', density = TRUE) # Detection parameters.
#dev.off()

# EXTRACT RELEVANT INFORMATION FROM MODEL FOR FURTHER ANALYSES --------

# extract species names from OTUs
aa <- str_split(spp_names, "_")
spp <- rep(NA, length(aa))

for(i in 1:length(aa)) {
  tmp <- aa[[i]][5:6]
  spp[i] <- paste(tmp[1], tmp[2])
}

spp

# extract classes names
classesNames <- sapply(spp_names, function(x){
  strsplit(x, split = "_")[[1]][2]
})

unique(classesNames)


# OCCUPANCY PREDICTIONS ---------------------------------------------------


## extract psi values
str(out$psi.samples)

# mean occupancy across all posterior draws by site
psi.mean <- apply(out$psi.samples, c(2), mean)
str(psi.mean)
cbind(psi.mean, spp)

# 95% credible interval acorss all posterior draws
psi.CI <- apply(out$psi.samples, c(2), function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% 
  as.data.frame %>%
  mutate(species = spp) %>%
  rename(lowerCI = `2.5%`, upperCI = `97.5%`)

str(psi.CI)
as.data.frame(psi.CI)


# combine occupancy mean, CI and site name
psi_preds <- cbind(psi.CI, psi.mean, obs.covs1$site)
psi_preds
str(psi_preds)


# plot occupancy probabilities 
occ_plot <- ggplot(psi_preds, aes(x = obs.covs1$site, y = psi.mean)) +  
  geom_point() + 
  geom_errorbar(data = psi_preds, 
                aes(x = obs.covs1$site,
                    ymin = lowerCI,
                    ymax = upperCI)) +
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0,
                             size = 8)) +
  coord_flip() + ylab("Occupancy Probability") + xlab("Site")

occ_plot

write.csv(psi_preds, here("outputs", "output_bullfrog", "occupancy_preds_MSOM_Bullfrogs.csv"))


# DETECTION PREDICTIONS ---------------------------------------------------

# extract detection probability estimates
# fitted values for detection probability
fit <- fitted(out)
str(fit)

str(fit$p.samples)

# extract the detection 95% credible interval for each species from all posterior draws
params_CI_p <- apply(fit$p.samples[,,], 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
})  %>% t %>% 
  as.data.frame %>%
  mutate(species = spp) %>%
  rename(lowerCI = `2.5%`, upperCI = `97.5%`)

str(params_CI_p)

# extract mean detection probability for each species from all posterior draws
params_mean_p <- apply(fit$p.samples[,,], 2, mean)

# combine detection mean, CI and site name
p_preds <- cbind(params_CI_p, params_mean_p, obs.covs1$site)

write.csv(p_preds, here("outputs", "output_bullfrog", "detection_preds_MSOM_Bullfrogs.csv"))


# plot detection probability predictions 
det_plot <- ggplot() +  
  geom_errorbar(data = p_preds, 
                aes(x = obs.covs1$site,
                    ymin = lowerCI,
                    ymax = upperCI)) +
  geom_point(data = p_preds, 
             aes(x = obs.covs1$site,
                 y = params_mean_p)) +
  coord_flip() + ylab("Detection Probability") + xlab("Site")

det_plot

