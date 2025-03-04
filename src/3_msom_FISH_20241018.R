
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

# fish species
spptoremove_forfish <- glgs %>% 
  as.data.frame() %>%
  dplyr::select(!contains(c("Spikein", "Mammalia", 
                            "Amphibia", "Reptilia", "Aves"))) %>%
  dplyr::select(!contains(c("Cololabis_saira", "Sardina_pilchardus",
                            "Sardinella_lemuru", "Engraulis_ringens",
                            "Micromesistius_poutassou", "Scophthalmus_maximus",
                            "Scomberomorus_niphonius", "Trichiurus_haumela"))) %>%
  dplyr::select(SampleID, PCR, site, contains("OTU")) %>% 
  rename_with(.cols = contains("OTU"), 
              function(x){paste0("OTU_", x)}) %>% 
  #dplyr::rename_with(.cols = starts_with(c("Teleostei")), 
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



# SUMMARIZE DATA AND PREPARE FOR MSOM -------------------------------------

# binarise read number and sum by site and sample replicate. For instance, if a species is detected in all 3 sample replicates and all 6 PCRs per sample replicate, the output will be 6,6,6 for that species and site.
glgs_sum1 <- glgs %>% 
  
  # retain only fish
  select(!contains(c("Spikein", "Mammalia", "Amphibia", "Reptilia", "Aves"))) |>
  
  # remove all seafish
  select(!contains(c("Cololabis_saira", "Sardina_pilchardus",
                     "Sardinella_lemuru", "Engraulis_ringens",
                     "Micromesistius_poutassou", "Scophthalmus_maximus",
                     "Scomberomorus_niphonius", "Trichiurus_haumela"))) |>
  
  # add "OTU_" to the beginning of all OTU columns
  rename_with(.cols = starts_with(c("Teleostei", "Actinopteri")), 
              function(x){paste0("OTU_", x)}) |> 
  
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
# 6 PCRs per sample replicate, the output will be 18 for that species and site. 
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

# naive detections
glgs_sum3 <- glgs_sum2 |> 
  # binarise
  mutate(across(starts_with("OTU_"), binarise1)) |> 
  dplyr::select(-one_of(spptoremove_forfish)) %>%
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
  select(site, altitude, river, tributary, longitude, latitude, riverOrigin, river_width, TEMP, water_velocity)

# add land cover vars
site.covs1 <- inner_join(site.covs, new_covs, by = "site")

# tributary is character; it needs to be numeric for random effect
site.covs1$tributary <- as.factor(site.covs1$tributary)
site.covs1$tributary <- as.numeric(site.covs1$tributary)

# river, riverOrigin and Reservoir_downstream are characters; it needs to be factor for model
site.covs1$river <- as.factor(site.covs1$river)
site.covs1$riverOrigin <- as.factor(site.covs1$riverOrigin)
site.covs1$Reservoir_downstream <- as.factor(site.covs1$Reservoir_downstream)

site.covs1$water_velocity[is.na(site.covs1$water_velocity)] <- mean(site.covs1$water_velocity, na.rm = T)

# create list for observation covariates
obs.covs <- glgs_sum2 |>
  select(Filter_vol, team, PrevDayRain)

obs.covs1 <- with(obs.covs, {
  list(Filter_vol = Filter_vol, team = team, PrevDayRain = PrevDayRain)})

# create a 3D array for species data (3 occasions)
species_for_y <- glgs_sum1 |>
  select(starts_with("OTU_")) |>
  mutate(across(starts_with("OTU_"), binarise1)) |>
  dplyr::select(-one_of(spptoremove_forfish)) %>% # remove species with a single site detection
  data.frame()

str(species_for_y)

spp_names <- colnames(species_for_y)
str(spp_names)

# Create a 3D array by site, species and visit, in this order (reqd by spOccupancy) 

#  number of species
N <- length(spp_names)
# number of replicates 
K <- 3
# Number of sites
J <- 101

y <- array(dim = c(N, J, K))

for(i in 1:N) {
  matrix1 <- matrix(species_for_y[,i], nrow = J, ncol = K, byrow = T)
  rownames(matrix1) <- glgs_sum2$site
  print(spp_names[i])
  y[i,,] <- matrix1
}

str(y)

# combine observations, occupancy covariates and detection covariates
# coordinates can be provided for spatial model 
multi_spp <- list(y = y, 
                  occ.covs = site.covs1,
                  det.covs = obs.covs1)
                  #coords = coords)
str(multi_spp)



# SET UP AND RUN THE MODEL ------------------------------------------------

# occupancy formula
occ.ms.formula <- ~ scale(TEMP) + scale(annual_precipitation) + 
  scale(Percent_of_forest_cover_1km) + scale(Human_modification_1km) + 
  scale(Patch_density_1km) + scale(Edge_density_1km) + 
  Reservoir_downstream + scale(river_width) + scale(water_velocity) + 
  riverOrigin + river +
  (1|tributary)

# detection formula
det.ms.formula <- ~ scale(Filter_vol) + scale(PrevDayRain) + (1|team) 

# model inits (initial values for parameters)
ms.inits <- list(alpha.comm = 0, 
                 beta.comm = 0, 
                 beta = 0, 
                 alpha = 0,
                 tau.sq.beta = 1, 
                 tau.sq.alpha = 1, 
                 z = apply(multi_spp$y, c(1, 2), max, na.rm = F))

# model priors 
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72), 
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

# RUN MODEL!
out.ms <- msPGOcc(occ.formula = occ.ms.formula, 
                  det.formula = det.ms.formula, 
                  data = multi_spp, 
                  inits = ms.inits, 
                  n.samples = 50000, 
                  priors = ms.priors, 
                  n.omp.threads = 3, 
                  verbose = TRUE, 
                  n.report = 2000, 
                  n.burn = 20000,
                  n.thin = 50, 
                  n.chains = 3)


summary(out.ms)
summary(out.ms, level = "community")

# predictive check
ppc.out <- ppcOcc(out.ms, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out)

# examine WAIC for compariosn with other models (if needed)
waicOcc(out.ms, by.sp = TRUE)

# traceplots to examine convergence
pdf(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","msom_fish_beta_traceplot_20241018.pdf"))
plot(out.ms, 'beta', density = TRUE) # Occupancy parameters.
dev.off()
pdf(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","msom_fish_alpha_traceplot_20241018.pdf"))
plot(out.ms, 'alpha', density = TRUE) # Detection parameters.
dev.off()


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
str(out.ms$psi.samples)

# mean occupancy across all posterior draws BY SPECIES
psi.mean.ms <- apply(out.ms$psi.samples, c(2), mean)
str(psi.mean.ms)
cbind(psi.mean.ms, spp)

# 95% credible interval acorss all posterior draws
psi.CI.ms <- apply(out.ms$psi.samples, c(2), function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% 
  as.data.frame %>%
  mutate(species = spp) %>%
  rename(lowerCI = `2.5%`, upperCI = `97.5%`)

str(psi.CI.ms)
as.data.frame(psi.CI.ms)

psi_preds <- cbind(psi.CI.ms, psi.mean.ms)

# order species from high-to-low MEAN occupancy probability
psi_preds$phylum <- classesNames
orderSpecies <- order(psi_preds$phylum, psi_preds$psi.mean.ms)

psi_preds$species <- factor(psi_preds$species,
                      levels = psi_preds$species[orderSpecies])

# plot occupancy probabilities for each species
occ_plot <- ggplot(psi_preds, aes(x = species, y = psi.mean.ms, col=phylum)) +  
  geom_point() + 
  geom_errorbar(data = psi_preds, 
                aes(x = species,
                    ymin = lowerCI,
                    ymax = upperCI,
                    color = phylum)) +
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0,
                             size = 8)) +
  coord_flip() + ylab("Occupancy Probability") + xlab("Species")

occ_plot

# save occupancy plot
ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","occupancy_MSOM_FISH_20241018.jpeg"), 
       width = 8,
       height = 10,
       units = "in",
       dpi = 300)

write.csv(psi_preds, here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","occupancy_preds_msom_FISH_20241018.csv"))


# mean occupancy across all posterior draws BY SPECIES AND SITE
psi.mean.bysite <- apply(out.ms$psi.samples, c(3,2), mean) |>
  as.data.frame() 

# 95% credible interval acorss all posterior draws
psi.lower95CI <- apply(out.ms$psi.samples, c(3, 2), function(x) quantile(x, 0.025)) |>
  as.data.frame()
psi.upper95CI <- apply(out.ms$psi.samples, c(3, 2), function(x) quantile(x, 0.975)) |>
  as.data.frame()


# extract standard deviation
psi.SD.bysite <- apply(out.ms$psi.samples, c(3,2), sd) |>
  as.data.frame()

psi.bysite <- cbind(
  psi.mean.bysite,
  psi.SD.bysite,
  psi.lower95CI,
  psi.upper95CI
) |>
  as.data.frame()

oldnames <- as.character(names(psi.bysite))
newnames <- as.character(c(paste0("Mean ", spp), paste0("SD ", spp), paste0("lower95CI ", spp), paste0("upper95CI ", spp)))

names(psi.bysite) <- c(newnames)
psi.bysite$site <- new_covs$site

# combine speciesxsite occupancy probabilities with site covariates and save the file
psi.bysite <- inner_join(site.covs1,psi.bysite, by = "site")


# extract predicted occupied status (z.samples)

#occ.ms.formula <- ~ scale(TEMP) + scale(annual_precipitation) + 
#  scale(Percent_of_forest_cover_1km) + scale(Human_modification_1km) + 
#  scale(Patch_density_1km) + scale(Edge_density_1km) + 
#  Reservoir_downstream + scale(river_width) + scale(water_velocity) + 
#  riverOrigin + river +
#  (1|tributary)

TEMP <- scale(multi_spp$occ.covs$TEMP)
annual_precipitation <- scale(multi_spp$occ.covs$annual_precipitation)
Percent_of_forest_cover_1km <- scale(multi_spp$occ.covs$Percent_of_forest_cover_1km)
Human_modification_1km <- scale(multi_spp$occ.covs$Human_modification_1km)
Patch_density_1km <- scale(multi_spp$occ.covs$Patch_density_1km)
Edge_density_1km <- scale(multi_spp$occ.covs$Edge_density_1km)
Reservoir_downstream <- multi_spp$occ.covs$Reservoir_downstream
river_width <- scale(multi_spp$occ.covs$river_width)
water_velocity <- scale(multi_spp$occ.covs$water_velocity)
riverOrigin <- as.factor(multi_spp$occ.covs$riverOrigin)
riverOrigin <- as.integer(riverOrigin)
riverOrigin[riverOrigin == 1] <- 0
riverOrigin[riverOrigin == 2] <- 1
river <-  multi_spp$occ.covs$riverOrigin
tributary <- multi_spp$occ.covs$tributary

X.0 <- cbind(1, TEMP, annual_precipitation, Percent_of_forest_cover_1km, Human_modification_1km, 
             Patch_density_1km, Edge_density_1km, Reservoir_downstream, river_width,
             water_velocity, riverOrigin, river, river, river,
             tributary)
str(X.0)

out.pred <- predict(out.ms, X.0, ignore.RE = F)

rich.samples <- apply(out.pred$z.0.samples, c(1, 3), sum)
rich.mean <- round(apply(rich.samples, 2, mean))
rich.sd <- apply(rich.samples, 2, sd)

pred.richness <- cbind(rich.mean, rich.sd)

psi.bysite <- cbind(psi.bysite, pred.richness)
colnames(psi.bysite)

write.csv(psi.bysite, here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","occupancy_bysite_msom_FISH_20241018.csv"))


# DETECTION PREDICTIONS ---------------------------------------------------

# extract detection probability estimates
# fitted values for detection probability
fit <- fitted(out.ms)
str(fit)

str(fit$p.samples)

# extract the detection 95% credible interval for each species from all posterior draws
params_CI_p <- apply(fit$p.samples[,,,], 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
})

str(params_CI_p)
params_CI_p1 <- as.data.frame(t(params_CI_p))

# extract mean detection probability for each species from all posterior draws
params_mean_p <- apply(fit$p.samples[,,,], 2, mean)

p_preds <- cbind(params_CI_p1, params_mean_p)

# create new variable "phylum" for plotting purposes
p_preds$species <- spp
p_preds$phylum <- classesNames

colnames(p_preds) <- c("lower95CI", "upper95CI", "mean",
                       "species", "phylum")

# order species from high-to-low MEAN detection probability
orderSpecies <- order(p_preds$phylum, p_preds$mean)

p_preds$species <- factor(p_preds$species,
                               levels = p_preds$species[orderSpecies])

write.csv(p_preds, here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","detection_preds_msom_FISH_20241018.csv"))

# plot detection probability predictions 
ggplot() +  
  geom_errorbar(data = p_preds, 
                aes(x = species,
                    ymin = lower95CI,
                    ymax = upper95CI,
                    color = phylum)) +
  geom_point(data = p_preds, 
             aes(x = species,
                 y = params_mean_p,
                 color = phylum)) +
  coord_flip() + ylab("Detection Probability") + xlab("Species")

# save detection probability figure
ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","detection_preds_msom_FISH_20241018.jpeg"), width = 8,
       height = 9,
       units = "in",
       dpi = 300)



# PREDICTIONS - EFFECTS OF VARIABLES ON OCCUPANCY ----------------------

# for FISH, we are interested in the distribution of invasive vs native
# spp, so we will extract that data first and specify colors for plotting

# Invasiveness...
spp_list <- read_tsv(here("data","SP_INFO_20240528.txt"))
spp_list <- as.data.frame(spp_list)
spp_list$species <- as.factor(spp_list$species)

spp1 <- as.data.frame(spp)
colnames(spp1) <- "species"
spp1$species <- as.factor(spp1$species)

spp_list1 <- left_join(spp1, spp_list, by = "species")
inv <- spp_list1$status

colorr <- function(x){
  case_when(
    inv == "invasive" ~ 'red',
    inv != "invasive" ~ 'blue'
  )
}

cols <- colorr(inv)


# first, let's remind ourselves about the model structure
#occ.ms.formula <- ~ scale(TEMP) + scale(annual_precipitation) + 
#  scale(Percent_of_forest_cover_1km) + scale(Human_modification_1km) + 
#  scale(Patch_density_1km) + scale(Edge_density_1km) + 
#  Reservoir_downstream + scale(river_width) + scale(water_velocity) + 
#  riverOrigin + river +
#  (1|tributary)

# 11 variables
#X.0 <- cbind(1, TEMP, annual_precipitation, Percent_of_forest_cover_1km, Human_modification_1km, 
#             Patch_density_1km, Edge_density_1km, Reservoir_downstream, river_width,
#             water_velocity, riverOrigin, river, river, river,
#             tributary)

# TEMP
# create prediction dataframe for variable of interest
TEMP.pred <- seq(min(scale(multi_spp$occ.covs$TEMP)), 
                  max(scale(multi_spp$occ.covs$TEMP)), 0.1)
range(TEMP.pred) # -2.11358  2.38642
X.1 <- cbind(1, TEMP.pred, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
out.pred1 <- predict(out.ms, X.1, ignore.RE = T)
str(out.pred1)

# plots all species 
b <- ggplot() +
  xlim(-2.12,2.8)

for (i in 1:length(spp)) {
  
  plot.dat1 <- data.frame(TEMP.pred,
                          mean.psi = apply(out.pred1$psi.0.samples[,i,], 2, mean), 
                          upper.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.975)),
                          lower.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.025)),
                          stringsAsFactors = FALSE)
  
  col <- cols[i]
  
  b <- b + geom_line(data=plot.dat1, aes(x=TEMP.pred, y=mean.psi), 
                     color = col) +
    labs(title = "Fish species occupancy probability",
         x =paste0("low", strrep(" ", 25), "water temperature (scaled)", strrep(" ", 25), "high"), y = "Occupancy probability")
}

b <- b +
  theme(plot.margin = unit(c(1,3,1,1), "lines")) +
  annotate(geom="text", x=2.6, y=0.7, label="invasive",
           color="red") +
  annotate(geom="text", x=2.6, y=0.65, label="native",
           color="blue")

print(b)

ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","Occ_FISH_TEMP_20241018.jpeg"), width = 6,
       height = 5,
       units = "in",
       dpi = 300)


# annual_precipitation
# create prediction dataframe for variable of interest
precipitation.pred <- seq(min(scale(multi_spp$occ.covs$annual_precipitation)), 
                          max(scale(multi_spp$occ.covs$annual_precipitation)), 0.1)
range(precipitation.pred) # -2.13619  2.66381
X.1 <- cbind(1, 0, precipitation.pred, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
out.pred1 <- predict(out.ms, X.1, ignore.RE = T)
str(out.pred1)

# plots all species 
b <- ggplot() +
  xlim(-2.2,3)

for (i in 1:length(spp)) {
  
  plot.dat1 <- data.frame(precipitation.pred,
                          mean.psi = apply(out.pred1$psi.0.samples[,i,], 2, mean), 
                          upper.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.975)),
                          lower.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.025)),
                          stringsAsFactors = FALSE)
  
  col <- cols[i]
  
  b <- b + geom_line(data=plot.dat1, aes(x=precipitation.pred, y=mean.psi), 
                     color = col) +
    labs(title = "Fish species occupancy probability",
         x =paste0("low", strrep(" ", 25), "Annual precipitation (scaled)", strrep(" ", 25), "high"), y = "Occupancy probability")
  
}

b <- b +
  theme(plot.margin = unit(c(1,3,1,1), "lines")) +
  annotate(geom="text", x=2.8, y=0.7, label="invasive",
           color="red") +
  annotate(geom="text", x=2.8, y=0.65, label="native",
           color="blue")

print(b)

ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","Occ_FISH_precipitation_20241018.jpeg"), width = 5,
       height = 4,
       units = "in",
       dpi = 300)


# FOREST COVER
# create prediction dataframe for variable of interest
forestcov.pred <- seq(min(scale(multi_spp$occ.covs$Percent_of_forest_cover_1km)), 
                      max(scale(multi_spp$occ.covs$Percent_of_forest_cover_1km)), 0.1)
range(forestcov.pred) # -2.0031944  0.9968056
X.1 <- cbind(1, 0, 0, forestcov.pred, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
out.pred1 <- predict(out.ms, X.1, ignore.RE = T)
str(out.pred1)

# plots all species 
b <- ggplot() +
  xlim(-2.01,1.5)

for (i in 1:length(spp)) {
  
  plot.dat1 <- data.frame(forestcov.pred,
                          mean.psi = apply(out.pred1$psi.0.samples[,i,], 2, mean), 
                          upper.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.975)),
                          lower.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.025)),
                          stringsAsFactors = FALSE)
  
  col <- cols[i]
  
  b <- b + geom_line(data=plot.dat1, aes(x=forestcov.pred, y=mean.psi), 
                     color = col) +
    labs(title = "Fish species occupancy probability",
         x =paste0("low", strrep(" ", 25), "Forest cover (scaled)", strrep(" ", 25), "high"), y = "Occupancy probability")
  
}

b <- b +
  annotate(geom="text", x=1.3, y=0.7, label="invasive",
           color="red") +
  annotate(geom="text", x=1.3, y=0.65, label="native",
           color="blue")


print(b)

ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","Occ_FISH_forestcover_20241018.jpeg"), width = 6,
       height = 5,
       units = "in",
       dpi = 300)


# HUMAN MODIFICATION
# create prediction dataframe for variable of interest
humanimp.pred <- seq(min(scale(multi_spp$occ.covs$Human_modification_1km)), 
                     max(scale(multi_spp$occ.covs$Human_modification_1km)), 0.1)
range(humanimp.pred) # -1.448473  2.251527
X.1 <- cbind(1, 0, 0, 0, humanimp.pred, 0, 0, 0, 0, 0, 0, 0, 0, 0)
out.pred1 <- predict(out.ms, X.1, ignore.RE = T)
str(out.pred1)

# plots all species 
b <- ggplot() +
  xlim(-1.5,2.5)

for (i in 1:length(spp)) {
  
  plot.dat1 <- data.frame(humanimp.pred,
                          mean.psi = apply(out.pred1$psi.0.samples[,i,], 2, mean), 
                          upper.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.975)),
                          lower.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.025)),
                          stringsAsFactors = FALSE)
  
  col <- cols[i]
  
  b <- b + geom_line(data=plot.dat1, aes(x=humanimp.pred, y=mean.psi), 
                     color = col) +
    #xlim(-1.5,3) +
    labs(title = "Fish species occupancy probability",
         x =paste0("low", strrep(" ", 25), "Human impact (scaled)", strrep(" ", 25), "high"), y = "Occupancy probability")
  
}

b <- b +
  theme(plot.margin = unit(c(1,1,1,1), "lines")) +
  annotate(geom="text", x=2.4, y=0.7, label="invasive",
           color="red") +
  annotate(geom="text", x=2.4, y=0.65, label="native",
           color="blue")


print(b)

ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","Occ_FISH_humanimpact_20241018.jpeg"), width = 6,
       height = 5,
       units = "in",
       dpi = 300)


# PATCH DENSITY
# create prediction dataframe for variable of interest
patchdens.pred <- seq(min(scale(multi_spp$occ.covs$Patch_density_1km)), 
                      max(scale(multi_spp$occ.covs$Patch_density_1km)), 0.1)
range(patchdens.pred) # -1.621791  2.578209
X.1 <- cbind(1, 0, 0, 0, 0, patchdens.pred, 0, 0, 0, 0, 0, 0, 0, 0)
out.pred1 <- predict(out.ms, X.1, ignore.RE = T)
str(out.pred1)

# plots all species 
b <- ggplot() +
  xlim(-1.8,3)

for (i in 1:length(spp)) {
  
  plot.dat1 <- data.frame(patchdens.pred,
                          mean.psi = apply(out.pred1$psi.0.samples[,i,], 2, mean), 
                          upper.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.975)),
                          lower.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.025)),
                          stringsAsFactors = FALSE)
  
  col <- cols[i]
  
  b <- b + geom_line(data=plot.dat1, aes(x=patchdens.pred, y=mean.psi), 
                     color = col) +
    labs(title = "Fish species occupancy probability",
         x =paste0("low", strrep(" ", 25), "Patch density (scaled)", strrep(" ", 25), "high"), y = "Occupancy probability")
  
}

b <- b +
  theme(plot.margin = unit(c(1,1,1,1), "lines")) +
  annotate(geom="text", x=2.8, y=0.7, label="invasive",
           color="red") +
  annotate(geom="text", x=2.8, y=0.65, label="native",
           color="blue")

print(b)

ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","Occ_FISH_patchdens_20241018.jpeg"), width = 5,
       height = 4,
       units = "in",
       dpi = 300)



# EDGE DENSITY
# create prediction dataframe for variable of interest
edgedens.pred <- seq(min(scale(multi_spp$occ.covs$Edge_density_1km)), 
                     max(scale(multi_spp$occ.covs$Edge_density_1km)), 0.1)
range(edgedens.pred) # -1.485431  2.614569
X.1 <- cbind(1, 0, 0, 0, 0, 0, edgedens.pred, 0, 0, 0, 0, 0, 0, 0)
out.pred1 <- predict(out.ms, X.1, ignore.RE = T)
str(out.pred1)

# plots all species 
b <- ggplot() +
  xlim(-1.5,3)

for (i in 1:length(spp)) {
  
  plot.dat1 <- data.frame(edgedens.pred,
                          mean.psi = apply(out.pred1$psi.0.samples[,i,], 2, mean), 
                          upper.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.975)),
                          lower.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.025)),
                          stringsAsFactors = FALSE)
  
  col <- cols[i]
  
  b <- b + geom_line(data=plot.dat1, aes(x=edgedens.pred, y=mean.psi), 
                     color = col) +
    labs(title = "Fish species occupancy probability",
         x =paste0("low", strrep(" ", 25), "Edge density (scaled)", strrep(" ", 25), "high"), y = "Occupancy probability")
  
}

b <- b +
  theme(plot.margin = unit(c(1,1,1,1), "lines")) +
  annotate(geom="text", x=2.85, y=0.7, label="invasive",
           color="red") +
  annotate(geom="text", x=2.85, y=0.65, label="native",
           color="blue")

print(b)

ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","Occ_FISH_edgedens_20241018.jpeg"), width = 5,
       height = 4,
       units = "in",
       dpi = 300)


# Reservoir downstream
# create prediction dataframe for variable of interest
reservdown.pred <- c(0,1)

X.1 <- cbind(1, 0, 0, 0, 0, 0, 0, reservdown.pred, 0, 0, 0, 0, 0, 0)
out.pred1 <- predict(out.ms, X.1, ignore.RE = T)
str(out.pred1)

# plots all species 
b <- ggplot() +
  scale_x_discrete(limits = c("0", "1"))
for (i in 1:length(spp)) {
  
  plot.dat1 <- data.frame(reservdown.pred,
                          mean.psi = apply(out.pred1$psi.0.samples[,i,], 2, mean), 
                          #upper.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                          #                  function(x) quantile(x, probs = 0.975)),
                          #lower.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                          #                  function(x) quantile(x, probs = 0.025)),
                          stringsAsFactors = FALSE)
  
  plot.dat1$reservdown.pred <- as.factor(plot.dat1$reservdown.pred)
  print(plot.dat1$mean.psi)
  
  newdata <- as.numeric(t(plot.dat1)[2,])
  
  newdata <- matrix(newdata, nrow=1) |>
    as.data.frame() |>
    dplyr::rename("0" = "V1", "1" = "V2")
  
  col <- cols[i]
  
  b <- b + geom_point(data=plot.dat1, aes(x=reservdown.pred, y=mean.psi), 
                      color = col) +
    geom_segment(data=newdata, aes(x = "0", xend = "1", 
                                   y = `0`, yend = `1`), 
                 color = col) +
    labs(title = "Fish species occupancy probability",
         x ="Reservoir downstream (0 = no; 1 = yes)", y = "Occupancy probability")
  
}

b <- b +
  annotate(geom="text", x=2.2, y=0.5, label="invasive",
           color="red") +
  annotate(geom="text", x=2.2, y=0.45, label="native",
           color="blue")

print(b)

ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","Occ_FISH_reservoir_downstream_20241018.jpeg"), width = 6,
       height = 5,
       units = "in",
       dpi = 300)


# river_width
# create prediction dataframe for variable of interest
rivwid.pred <- seq(min(scale(multi_spp$occ.covs$river_width)), 
                   max(scale(multi_spp$occ.covs$river_width)), 0.1)
range(rivwid.pred) # -0.5367418  6.1632582
X.1 <- cbind(1, 0, 0, 0, 0, 0, 0, 0, rivwid.pred, 0, 0, 0, 0, 0)
out.pred1 <- predict(out.ms, X.1, ignore.RE = T)
str(out.pred1)

# plots all species 
b <- ggplot() +
  xlim(-0.55,6.4)

for (i in 1:length(spp)) {
  
  plot.dat1 <- data.frame(rivwid.pred,
                          mean.psi = apply(out.pred1$psi.0.samples[,i,], 2, mean), 
                          upper.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.975)),
                          lower.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.025)),
                          stringsAsFactors = FALSE)
  
  col <- cols[i]
  
  b <- b + geom_line(data=plot.dat1, aes(x=rivwid.pred, y=mean.psi), 
                     color = col) +
    labs(title = "Fish species occupancy probability",
         x =paste0("low", strrep(" ", 25), "River width (scaled)", strrep(" ", 25), "high"), y = "Occupancy probability")
  
}

b <- b +
  annotate(geom="text", x=6.3, y=0.74, label="invasive",
           color="red") +
  annotate(geom="text", x=6.3, y=0.69, label="native",
           color="blue")


print(b)

ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","Occ_FISH_river_width_20241018.jpeg"), width = 6,
       height = 5,
       units = "in",
       dpi = 300)


# water_velocity
# create prediction dataframe for variable of interest
velocity.pred <- seq(min(scale(multi_spp$occ.covs$water_velocity)), 
                   max(scale(multi_spp$occ.covs$water_velocity)), 0.1)
range(velocity.pred) # -1.741327  3.558673
X.1 <- cbind(1, 0, 0, 0, 0, 0, 0, 0, 0, velocity.pred, 0, 0, 0, 0)
out.pred1 <- predict(out.ms, X.1, ignore.RE = T)
str(out.pred1)

# plots all species 
b <- ggplot() +
  xlim(-1.75,4)

for (i in 1:length(spp)) {
  
  plot.dat1 <- data.frame(velocity.pred,
                          mean.psi = apply(out.pred1$psi.0.samples[,i,], 2, mean), 
                          upper.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.975)),
                          lower.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                                            function(x) quantile(x, probs = 0.025)),
                          stringsAsFactors = FALSE)
  
  col <- cols[i]
  
  b <- b + geom_line(data=plot.dat1, aes(x=velocity.pred, y=mean.psi), 
                     color = col) +
    labs(title = "Fish species occupancy probability",
         x =paste0("low", strrep(" ", 25), "Water velocity (scaled)", strrep(" ", 25), "high"), y = "Occupancy probability")
  
}

b <- b +
  annotate(geom="text", x=3.8, y=0.74, label="invasive",
           color="red") +
  annotate(geom="text", x=3.8, y=0.69, label="native",
           color="blue")


print(b)

ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","Occ_FISH_velocity_20241018.jpeg"), width = 6,
       height = 5,
       units = "in",
       dpi = 300)


# RIVER ORIGIN
# create prediction dataframe for variable of interest
riverOri.pred <- c(0,1)

X.1 <- cbind(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, riverOri.pred, 0, 0, 0)
out.pred1 <- predict(out.ms, X.1, ignore.RE = T)
str(out.pred1)

# plots all species 
b <- ggplot()

for (i in 1:length(spp)) {
  
  plot.dat1 <- data.frame(riverOri.pred,
                          mean.psi = apply(out.pred1$psi.0.samples[,i,], 2, mean), 
                          #upper.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                          #                  function(x) quantile(x, probs = 0.975)),
                          #lower.psi = apply(out.pred1$psi.0.samples[,i,], 2, 
                          #                  function(x) quantile(x, probs = 0.025)),
                          stringsAsFactors = FALSE)
  
  plot.dat1$riverOri.pred <- as.factor(plot.dat1$riverOri.pred)
  print(plot.dat1$mean.psi)
  
  newdata <- as.numeric(t(plot.dat1)[2,])
  
  newdata <- matrix(newdata, nrow=1) |>
    as.data.frame() |>
    dplyr::rename("0" = "V1", "1" = "V2")
  
  col <- cols[i]
  
  b <- b + geom_point(data=plot.dat1, aes(x=riverOri.pred, y=mean.psi), 
                      color = col) +
    geom_segment(data=newdata, aes(x = "0", xend = "1", 
                                   y = `0`, yend = `1`), 
                 color = col) +
    labs(title = "Fish species occupancy probability",
         x ="River origin (0 = inside; 1 = outside)", y = "Occupancy probability")
  
}

b <- b +
  annotate(geom="text", x=2.2, y=0.5, label="invasive",
           color="red") +
  annotate(geom="text", x=2.2, y=0.45, label="native",
           color="blue")

print(b)

ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","Occ_FISH_riverOrigin_20241018.jpeg"), width = 6,
       height = 5,
       units = "in",
       dpi = 300)



# Plot the occupancy covariate effects ---------------------------------------------------------------
beta_comm <- as.data.frame(out.ms$beta.comm.samples)
effects_summary <- beta_comm %>%
  summarise(
    TEMP_mean = mean(beta_comm$`scale(TEMP)`),
    TEMP_lci = quantile(beta_comm$`scale(TEMP)`, 0.025),
    TEMP_uci = quantile(beta_comm$`scale(TEMP)`, 0.975),
    annual_precipitation_mean = mean(beta_comm$`scale(annual_precipitation)`), 
    annual_precipitation_lci = quantile(beta_comm$`scale(annual_precipitation)`, 0.025),
    annual_precipitation_uci = quantile(beta_comm$`scale(annual_precipitation)`, 0.975),
    Percent_of_forest_cover_1km_mean = mean(beta_comm$`scale(Percent_of_forest_cover_1km)`),
    Percent_of_forest_cover_1km_lci = quantile(beta_comm$`scale(Percent_of_forest_cover_1km)`, 0.025),
    Percent_of_forest_cover_1km_uci = quantile(beta_comm$`scale(Percent_of_forest_cover_1km)`, 0.975),
    Human_modification_1km_mean = mean(beta_comm$`scale(Human_modification_1km)`),
    Human_modification_1km_lci = quantile(beta_comm$`scale(Human_modification_1km)`, 0.025),
    Human_modification_1km_uci = quantile(beta_comm$`scale(Human_modification_1km)`, 0.975),
    Patch_density_1km_mean = mean(beta_comm$`scale(Patch_density_1km)`),
    Patch_density_1km_lci = quantile(beta_comm$`scale(Patch_density_1km)`, 0.025),
    Patch_density_1km_uci = quantile(beta_comm$`scale(Patch_density_1km)`, 0.975),
    Edge_density_1km_mean = mean(beta_comm$`scale(Edge_density_1km)`),
    Edge_density_1km_lci = quantile(beta_comm$`scale(Edge_density_1km)`, 0.025),
    Edge_density_1km_uci = quantile(beta_comm$`scale(Edge_density_1km)`, 0.975), 
    Reservoir_downstreamyes_mean = mean(Reservoir_downstreamyes),
    Reservoir_downstreamyes_lci = quantile(Reservoir_downstreamyes, 0.025),
    Reservoir_downstreamyes_uci = quantile(Reservoir_downstreamyes, 0.975),
    river_width_mean = mean(beta_comm$`scale(river_width)`), 
    river_width_lci = quantile(beta_comm$`scale(river_width)`, 0.025),
    river_width_uci = quantile(beta_comm$`scale(river_width)`, 0.975),
    water_velocity_mean = mean(beta_comm$`scale(water_velocity)`), 
    water_velocity_lci = quantile(beta_comm$`scale(water_velocity)`, 0.025),
    water_velocity_uci = quantile(beta_comm$`scale(water_velocity)`, 0.975),
    streamOriginoutsideGLGS_mean = mean(riverOriginoutsideGLGS),
    streamOriginoutsideGLGS_lci = quantile(riverOriginoutsideGLGS, 0.025),
    streamOriginoutsideGLGS_uci = quantile(riverOriginoutsideGLGS, 0.975),
    riverdy_mean = mean(riverdy),
    riverdy_lci = quantile(riverdy, 0.025),
    riverdy_uci = quantile(riverdy, 0.975),
    riverlc_mean = mean(riverlc),
    riverlc_lci = quantile(riverlc, 0.025),
    riverlc_uci = quantile(riverlc, 0.975),
    rivernu_mean = mean(rivernu),
    rivernu_lci = quantile(rivernu, 0.025),
    rivernu_uci = quantile(rivernu, 0.975),
  )
effects_long <- effects_summary %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_(mean|lci|uci)"
  )
ordered_variables <- c(
  "water_velocity", "river_width", "Reservoir_downstreamyes", "annual_precipitation", 
  "Edge_density_1km", "Human_modification_1km", "TEMP", 
  "Patch_density_1km", "Percent_of_forest_cover_1km", "riverdy", "riverlc", 
  "rivernu", "streamOriginoutsideGLGS"
)
effects_long$variable <- factor(effects_long$variable, levels = ordered_variables)

ggplot(effects_long, aes(x = variable, ymin = lci, ymax = uci)) +
  geom_errorbar() +
  geom_hline(aes(yintercept = 0), color = "red") +
  theme_bw() +
  coord_flip() +
  labs(title = "MSOM - Posterior Distribution of Environmental Variable Effects",
       y = "Effect Size",
       x = "Covariate") +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  ) +
  scale_y_continuous(breaks = seq(-15, 2, 5), limits = c(-15, 4))
ggsave(here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","MSOM_Covariates_FISH_20241018.jpeg"), width = 5,
       height = 2,
       units = "in",
       dpi = 300)

beta <- as.data.frame(out.ms$beta.samples)
beta_long <- beta %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "effect") %>%
  separate(variable, into = c("covariate", "species"), sep = "-") %>%
  mutate(species = gsub("\\(Intercept\\)", "Intercept", species)) %>%
  group_by(species, covariate) %>%
  summarize(
    mean_effect = mean(effect, na.rm = TRUE),
    lower_95 = quantile(effect, 0.025, na.rm = TRUE),
    upper_95 = quantile(effect, 0.975, na.rm = TRUE)
  ) %>% 
  filter(covariate != '(Intercept)')
spp_df <- data.frame(
  species = paste0("sp", 1:71),
  latin_name = spp
)
beta_long <- left_join(beta_long, spp_df, by = "species") %>% 
  mutate(species = latin_name) %>%
  select(-latin_name)
write.csv(beta_long, here("outputs","msom_outputs_Viorel covs","FISH_outputs_20241018","envCov_perSP_table_msom_FISH.csv"))





save(list = ls(), file = here("data","Rdata files","msom_FISH_20241018_Vioreldata.Rdata"))
