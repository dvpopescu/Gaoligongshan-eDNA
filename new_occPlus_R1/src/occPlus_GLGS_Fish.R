
# LOAD DATA ------

#devtools::install_github("AlexDiana/occPlus", build_vignettes = TRUE,
#                         force = TRUE)

library(occPlus)

#vignette("vignette", package = "occPlus")

library(tidyverse)

filename <- "OTUtable_12S_toSP_GLG23_20240528.txt"

data <- read_tsv(file = filename)
str(data)
View(data)


new_covs <- read_tsv("23GLG_covs_info_20240927.txt")
zones <- read_tsv("23GLG_zones_info_20240903.txt")
climate_covs <- read_csv("GLGS_climate.csv")
climate_covs <- climate_covs %>% select(site, annual_precipitation, mean_temperature)

new_covs <- left_join(new_covs, zones, by = "site")
new_covs <- left_join(new_covs, climate_covs, by = "site")

str(new_covs)


binarise1 <- function(x, na.rm = FALSE) vegan::decostand(x, method = "pa")

# species to remove
spptoremove <- data %>%
  as.data.frame() %>%
  # Mammalia
  #dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Aves", "Amphibia", "Reptilia"))) %>%
  # Aves
  #dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Mammalia", "Amphibia", "Reptilia"))) %>%
  # Amphia + Reptilia
  #dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Aves", "Mammalia"))) %>%
  # Fish
  dplyr::select(!contains(c("Spikein", "Mammalia", "Aves", "Amphibia", "Reptilia"))) %>%
  dplyr::select(!contains(c("Cololabis_saira", "Sardina_pilchardus",
                            "Sardinella_lemuru", "Engraulis_ringens",
                            "Micromesistius_poutassou", "Scophthalmus_maximus",
                            "Scomberomorus_niphonius", "Trichiurus_haumela"))) %>%
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

str(spptoremove)

# all species OTU table (minus species that were only detected at 1 site)
OTUtable <- data %>%
  as.data.frame() %>%
  # Mammalia
  #dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Aves", "Amphibia", "Reptilia"))) %>%
  # Aves
  #dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Mammalia", "Amphibia", "Reptilia"))) %>%
  # Amphia + Reptilia
  #dplyr::select(!contains(c("Spikein", "Actinopteri", "Teleostei", "Aves", "Mammalia"))) %>%
  # Fish
  dplyr::select(!contains(c("Spikein", "Mammalia", "Aves", "Amphibia", "Reptilia"))) %>%
  dplyr::select(!contains(c("Cololabis_saira", "Sardina_pilchardus",
                            "Sardinella_lemuru", "Engraulis_ringens",
                            "Micromesistius_poutassou", "Scophthalmus_maximus",
                            "Scomberomorus_niphonius", "Trichiurus_haumela"))) %>%
  dplyr::filter(status == "sample") %>%
  dplyr::select(contains("OTU")) %>% 
  dplyr::select(-one_of(spptoremove)) %>%
  as.matrix

str(OTUtable)
View(OTUtable)


# covariates for each site
data_infos <- data %>% 
  mutate(primer = "12S") |>
  arrange(SampleID) %>% 
  group_by(SampleID) %>%
  dplyr::filter(row_number() == 1) %>% 
  mutate(streamOrigin = case_when(
    is.na(DistFromPAedge) ~ "outsideGLGS",
    TRUE ~ "insideGLGS"
  )) %>% 
  rename(Site = site,
         Sample = SampleID,
         Primer = primer) %>%
  dplyr::select(Site, Primer, latitude, longitude, river, tributary, team,
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


data1 <- list(info = data_infos1,
             OTU = OTUtable)


# RUN MODEL -----

# Fish 
fitmodel2  <- runOccPlus(data1,
                         d = 2,
                         occCovariates = c(),
                         ordCovariates = c("TEMP", "annual_precipitation",
                                           "Human_modification_1km",
                                           "Percent_of_forest_cover_1km",
                                           "Patch_density_1km",
                                           "Edge_density_1km",
                                           "Reservoir_downstream",
                                           "river_width",
                                           "water_velocity",
                                           "streamOrigin",
                                           "river"),
                         detCovariates = c("Filter_vol","PrevDayRain"),
                         numSamples = 1000)



fitmodel3  <- runOccPlus(data1,
                         d = 3,
                         occCovariates = c(),
                         ordCovariates = c("TEMP", "annual_precipitation",
                                           "Human_modification_1km",
                                           "Percent_of_forest_cover_1km",
                                           "Patch_density_1km",
                                           "Edge_density_1km",
                                           "Reservoir_downstream",
                                           "river_width",
                                           "water_velocity",
                                           "streamOrigin",
                                           "river"),
                         detCovariates = c("Filter_vol","PrevDayRain"),
                         numSamples = 1000)


fitmodel4  <- runOccPlus(data1,
                         d = 4,
                         occCovariates = c(),
                         ordCovariates = c("TEMP", "annual_precipitation",
                                           "Human_modification_1km",
                                           "Percent_of_forest_cover_1km",
                                           "Patch_density_1km",
                                           "Edge_density_1km",
                                           "Reservoir_downstream",
                                           "river_width",
                                           "water_velocity",
                                           "streamOrigin",
                                           "river"),
                         detCovariates = c("Filter_vol","PrevDayRain"),
                         numSamples = 1000)


fitmodel5  <- runOccPlus(data1,
                         d = 5,
                         occCovariates = c(),
                         ordCovariates = c("TEMP", "annual_precipitation",
                                           "Human_modification_1km",
                                           "Percent_of_forest_cover_1km",
                                           "Patch_density_1km",
                                           "Edge_density_1km",
                                           "Reservoir_downstream",
                                           "river_width",
                                           "water_velocity",
                                           "streamOrigin",
                                           "river"),
                         detCovariates = c("Filter_vol","PrevDayRain"),
                         numSamples = 1000)

### model comparison
library(loo)

log_lik2 <- extract_log_lik(fitmodel2$vb_fit, parameter_name = "log_lik", merge_chains = FALSE)
log_lik3 <- extract_log_lik(fitmodel3$vb_fit, parameter_name = "log_lik", merge_chains = FALSE)
log_lik4 <- extract_log_lik(fitmodel4$vb_fit, parameter_name = "log_lik", merge_chains = FALSE)
log_lik5 <- extract_log_lik(fitmodel5$vb_fit, parameter_name = "log_lik", merge_chains = FALSE)


loo2 <- loo::loo(log_lik2)
loo3 <- loo::loo(log_lik3)
loo4 <- loo::loo(log_lik4)
loo5 <- loo::loo(log_lik5)


comparison <- loo_compare(loo2, loo3, loo4, loo5)
print(comparison)


saveRDS(fitmodel5, "occuPlus_Fish_d5_9June2025.rda")

# Investigate Model Outputs -----------------------------------------------

fitmodel <- readRDS("occuPlus_Fish_d5_9June2025.rda")

# OUTPUT  ---------

colnames(fitmodel$X_theta)
colnames(fitmodel$X_ord)
colnames(fitmodel$X_psi)

# plotOccupancyCovariates(fitmodel,
#                         covName = "")

#plot effects of covariates on detection 
plotDetectionCovariates(fitmodel,
                        covName = "Filter_vol")

plotDetectionCovariates(fitmodel,
                        covName = "PrevDayRain")


#plot effects of covariates on community ordination
plotOrdinationCovariates(fitmodel,
                         covName = "Human_modification_1km")

plotOrdinationCovariates(fitmodel,
                         covName = "streamOriginoutsideGLGS")

plotOrdinationCovariates(fitmodel,
                         covName = "annual_precipitation")

plotOrdinationCovariates(fitmodel,
                         covName = "mean_temperature")

plotOrdinationCovariates(fitmodel,
                         covName = "Patch_density_1km")

plotOrdinationCovariates(fitmodel,
                         covName = "Percent_of_forest_cover_1km")

plotOrdinationCovariates(fitmodel,
                         covName = "Edge_density_1km")

plotOrdinationCovariates(fitmodel,
                         covName = "Regional_DivisionNorthwest")

plotOrdinationCovariates(fitmodel,
                         covName = "Regional_DivisionSouthwest")

# plot occupancy by species (95% Credible Interval)
plotOccupancyRates(fitmodel)
# can use plotOccupancyRates(fitmodel, idx_species = 1:20) to plot only the first 20 species

# plot Stage 1 detection (95% Credible Interval)
plotCollectionRates(fitmodel)

# plot Stage 2 detection (95% Credible Interval)
plotDetectionRates(fitmodel)

# plot Stage 1 (water sampling) false positive rates (95% Credible Interval)
plotStage1FPRates(fitmodel)

# plot Stage 1 (water sampling) false positive rates (95% Credible Interval)
plotStage2FPRates(fitmodel)

# plot read intensity for true and false positives
plotReadIntensity(fitmodel) 

# plot covariance matrix
plotCovarianceMatrix(fitmodel)

plotSigElementsCovMatrix(fitmodel)


###################################################
### EXPLORE OUTPUTS AND EXTRACT RELEVANT METRICS

### you need to create or re-create the "data1" object from above to get 
### the OTU names... getting the OTU's can probably be done in a different way as well

matrix_of_draws <- fitmodel$matrix_of_draws

nsites <- length(unique(data$info$Site))
nspecies <- ncol(data$OTU)
niter <- nrow(matrix_of_draws)
ncov_psi <- ncol(fitmodel$X_psi)
ncov_ord <- ncol(fitmodel$X_ord)
ncov_theta <- ncol(fitmodel$X_theta)

X_psi <- fitmodel$X_psi
X_ord <- fitmodel$X_ord
X_theta <- fitmodel$X_theta
d <- 5

U_output0 <- 
  matrix_of_draws[,grepl("U\\[", colnames(matrix_of_draws))]
L_output0 <- 
  matrix_of_draws[,grepl("L\\[", colnames(matrix_of_draws))]
beta_psi_output0 <- 
  matrix_of_draws[,grepl("beta_psi\\[", colnames(matrix_of_draws))]
beta_ord_output0 <- 
  matrix_of_draws[,grepl("beta_ord\\[", colnames(matrix_of_draws))]
beta_theta_output0 <- 
  matrix_of_draws[,grepl("beta_theta\\[", colnames(matrix_of_draws))]
theta0_output <- 
  matrix_of_draws[,grepl("theta0\\[", colnames(matrix_of_draws))]
p_output <- 
  matrix_of_draws[,grepl("p\\[", colnames(matrix_of_draws))]
q_output <- 
  matrix_of_draws[,grepl("q\\[", colnames(matrix_of_draws))]

U_output <- array(NA, dim = c(niter, nsites, d))
for(iter in 1:niter){
  U_output[iter,,] <- matrix(U_output0[iter,], nsites, d, byrow = F)
}

L_output <- array(NA, dim = c(niter, d, nspecies))
for(iter in 1:niter){
  L_output[iter,,] <- matrix(L_output0[iter,], d, nspecies, byrow = F)
}

beta_psi_output <- array(NA, dim = c(niter, ncov_psi, nspecies))
for(iter in 1:niter){
  beta_psi_output[iter,,] <- matrix(beta_psi_output0[iter,], ncov_psi, nspecies, byrow = F)
}

beta_ord_output <- array(NA, dim = c(niter, ncov_ord, d))
for(iter in 1:niter){
  beta_ord_output[iter,,] <- matrix(beta_ord_output0[iter,], ncov_ord, d, byrow = F)
}

beta_theta_output <- array(NA, dim = c(niter, ncov_theta, nspecies))
for(iter in 1:niter){
  beta_theta_output[iter,,] <- 
    matrix(beta_theta_output0[iter,], ncov_theta, nspecies, byrow = F)
}

# REPARAMETRISE PROBABILITIES ----------

niter <- nrow(L_output)

L_output_reparam <- L_output
U_output_reparam <- U_output
# E_output_reparam <- E_output
beta_ord_output_reparam <- beta_ord_output

d <- dim(L_output)[2]

for (iter in 1:niter) {
  print(iter)
  
  if(d == 1){
    
    L1 <- L_output[iter,1,1]
    L_output_reparam[iter,1,] <- L_output[iter,1,] / L1
    U_output_reparam[iter,,1] <- U_output[iter,,1] * L1
    beta_ord_output_reparam[iter,,1] <- beta_ord_output[iter,,1] * L1
    
  } else {
    
    L_current <- L_output[iter,,]
    # E_current <- E_output[iter,,]
    U_current <- U_output[iter,,]
    beta_ord_current <- beta_ord_output[iter,,]
    
    qr_decomp <- qr(L_current)
    Q_current <- qr.Q(qr_decomp)
    R_current <- qr.R(qr_decomp)
    
    Q2 <- Q_current %*% diag(diag(R_current), nrow = d)
    invQ2 <- diag(1 / diag(R_current), nrow = d) %*% t(Q_current)
    
    betapsiord_new <- beta_ord_current %*% Q2
    # E_new <- E_current %*% Q2
    L_new <- invQ2 %*% L_current
    U_new <- U_current %*% Q2
    
    L_output_reparam[iter,,] <- L_new
    # E_output_reparam[iter,] <- E_new
    U_output_reparam[iter,,] <- U_new
    beta_ord_output_reparam[iter,,] <- betapsiord_new
  }
  
}

# OCCUPANCY PROBABILITY --------

# create classes names
{
  speciesNames <- colnames(data1$OTU)
  speciesNames <- make.unique(speciesNames)
  
  classesNames <- sapply(speciesNames, function(x){
    strsplit(x, split = "_")[[1]][1]
  })
}
##

siteNames <- data1$info$Site[!duplicated(data1$info$Site)]

#
occProbs_output <- array(NA, dim = c(niter, nsites, nspecies))
for (i in 1:nsites) {
  for (iter in 1:niter) {
    occProbs_output[iter,i,] <-  logistic(
      X_psi[i,] %*% beta_psi_output[iter,,] + 
        U_output[iter,i,] %*% L_output[iter,,]
    )
  }
}


# extract species names from OTUs
species <- colnames(data1$OTU)

aa <- str_split(species, "_")
spp <- rep(NA, length(aa))

for(i in 1:length(aa)) {
  tmp <- aa[[i]][4:5]
  spp[i] <- paste(tmp[1], tmp[2])
}
spp

occ <- cbind(
  (apply(occProbs_output, c(3), function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>% t),
  apply(occProbs_output, c(3), mean)) |>
  as.data.frame() |>
  mutate(species = spp, meanocc = V3) |>
  select(-V3)

occ 

occ$phylum <- classesNames
orderSpecies <- order(occ$phylum, occ$meanocc)

occ$species <- factor(occ$species,
                      levels = occ$species[orderSpecies])

occ_plot <- ggplot(occ, aes(x = species, y = meanocc, col=phylum)) +  
  geom_point() + 
  geom_errorbar(data = occ, 
                aes(x = species,
                    ymin = `2.5%`,
                    ymax = `97.5%`,
                    color = phylum)) +
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0,
                             size = 8)) +
  coord_flip() + ylab("Occupancy Probability") + xlab("Species")

occ_plot

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "occupancy_occPlus_FISH_20241018_5factors.jpeg"), width = 8,
       height = 8,
       units = "in",
       dpi = 300)

write.csv(occ, here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "occupancy_occPlus_FISH_20241018_5factors.csv"))

# occupancy by site for spatial predictions

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


occProbs_mode_by_site <- matrix(NA, nrow = nsites, ncol = nspecies)
for (i in 1:nsites) {
  for (j in 1:nspecies) {
    occProbs_mode_by_site[i,j] <- getmode(occProbs_output[,i,j])
  } 
}
str(occProbs_mode_by_site)


occProbs_mean_by_site <- matrix(NA, nrow = nsites, ncol = nspecies)
for (i in 1:nsites) {
  for (j in 1:nspecies) {
    occProbs_mean_by_site[i,j] <- mean(occProbs_output[,i,j])
  } 
}
str(occProbs_mean_by_site)

occProbs_CI_by_site <- array(NA, dim = c(2, nsites, nspecies))
for (i in 1:nsites) {
  for (j in 1:nspecies) {
    occProbs_CI_by_site[,i,j] <- quantile(occProbs_output[,i,j], probs = c(0.025, 0.975))
  } 
}

str(occProbs_CI_by_site)
head(occProbs_CI_by_site)

occProbs_CI_lower_by_site <- matrix(NA, nrow = nsites, ncol = nspecies)
occProbs_CI_upper_by_site <- matrix(NA, nrow = nsites, ncol = nspecies)

for (i in 1:nsites) {
  for (j in 1:nspecies) {
    occProbs_CI_lower_by_site[i, j] <- quantile(occProbs_output[, i, j], probs = 0.025)
    occProbs_CI_upper_by_site[i, j] <- quantile(occProbs_output[, i, j], probs = 0.975)
  }
}

occProbs_SD_by_site <- matrix(NA, nrow = nsites, ncol = nspecies)
for (i in 1:nsites) {
  for (j in 1:nspecies) {
    occProbs_SD_by_site[i,j] <- sd(occProbs_output[,i,j])
  } 
}

str(occProbs_SD_by_site)

occ <- cbind(
  occProbs_mean_by_site,
  occProbs_CI_lower_by_site,
  occProbs_CI_upper_by_site,
  occProbs_SD_by_site,
  occProbs_mode_by_site) |>
  as.data.frame() |>
  mutate(site = siteNames)

oldnames <- as.character(names(occ))
newnames <- as.character(c(paste0("Mean ", spp), paste0("lower95CI ", spp), paste0("upper95CI ", spp), paste0("SD ", spp), paste0("Mode ", spp)))

names(occ) <- c(newnames, "site")

occ

write.csv(occ, here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "occupancybysite_occPlus_FISH_20241018_5factors.csv"))

# DETECTION PROBABILITIES -------

# create classes names
{
  speciesNames <- colnames(data1$OTU)
  speciesNames <- make.unique(speciesNames)
  
  classesNames <- sapply(speciesNames, function(x){
    strsplit(x, split = "_")[[1]][1]
  })
}

# 1 stage probabilities output
{
  
  beta_theta_output1 <- aperm(beta_theta_output, c(1,3,2))
  str(beta_theta_output1)
  
  params_CI_theta <- apply(beta_theta_output1, 2, function(x) {
    logistic(quantile(x, probs = c(0.025, 0.975)))
  }) %>% t 
  
  params_CI_theta0 <- apply(theta0_output, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>% t 
  
  params_mean_theta <- apply(beta_theta_output1, 2, function(x) {
    logistic(mean(x))
  })
  
  params_CI <- data.frame(theta = params_CI_theta,
                          theta0 = params_CI_theta0,
                          mean1st = params_mean_theta)
  
  params_CI$species <- spp #1:ncol(data$OTU)
  
  data_plot <- params_CI
  
  data_plot$phylum <- classesNames
  
  orderSpecies <- order(data_plot$phylum, data_plot$theta.2.5.)
  
  data_plot$species <- factor(data_plot$species,
                              levels = data_plot$species[orderSpecies])
  
  theta_plot <- ggplot() +  
    geom_errorbar(data = data_plot, 
                  aes(x = species,
                      ymin = theta.2.5.,
                      ymax = theta.97.5.,
                      color = phylum)) + 
    geom_errorbar(data = data_plot, 
                  aes(x = species,
                      ymin = theta0.2.5.,
                      ymax = theta0.97.5.,
                      color = phylum), linetype = "dashed") + 
    # geom_hline(aes(yintercept = 0),
    # color = "red") +
    theme_bw() + 
    theme(
      axis.text = element_text(angle = 0,
                               size = 8)
    ) + coord_flip() + ylab("Collection Probability") + xlab("Species")
  
  theta_plot
  
  ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "1stagerate_FISH_20241018_5factors.jpeg"), width = 8,
         height = 8,
         units = "in",
         dpi = 300)
  
}

row.names(data_plot) <- NULL

write.csv(data_plot, here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "1stagerate_occPlus_FISH_20241018_5factors.csv"))

# 2 stage probabilities output
{
  
  params_CI_p <- apply(p_output, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>% t 
  
  params_CI_q <- apply(q_output, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>% t 
  
  params_mean_p <- apply(p_output, 2, function(x) {
    logistic(mean(x))
  })
  
  params_CI <- data.frame(p = params_CI_p,
                          q = params_CI_q,
                          mean2nd = params_mean_p)
  
  params_CI$species <- spp #1:ncol(data$OTU)
  
  data_plot <- params_CI
  
  data_plot$phylum <- classesNames
  
  orderSpecies <- order(data_plot$phylum, data_plot$p.2.5.)
  
  data_plot$species <- factor(data_plot$species,
                              levels = data_plot$species[orderSpecies])
  
  pq_plot <- ggplot() +  
    geom_errorbar(data = data_plot, 
                  aes(x = species,
                      ymin = p.2.5.,
                      ymax = p.97.5.,
                      color = phylum)) + 
    geom_errorbar(data = data_plot, 
                  aes(x = species,
                      ymin = q.2.5.,
                      ymax = q.97.5.,
                      color = phylum), linetype = "dashed") + 
    theme(
      axis.text = element_text(angle = 0,
                               size = 8)
    ) + coord_flip() + ylab("Amplification Probability") + xlab("Species")
  
  pq_plot
  
  ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "2stagerate_FISH_20241018_5factors.jpeg"), width = 8,
         height = 8,
         units = "in",
         dpi = 300)
  
}

row.names(data_plot) <- NULL

write.csv(data_plot, here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "2stagerate_occPlus_FISH_20241018_5factors.csv"))

# ORDINATION COVARIATE OUTPUT ------

param <- "beta_ord"

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws)),drop=F]

data_plot <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% 
  as.data.frame

data_plot$x <- colnames(fitmodel$X_ord)

# Factor 1
data_plot1 <- data_plot[1:(nrow(data_plot)/5),]

ggplot(data_plot1, aes(x = x,
                       ymin = `2.5%`,
                       ymax = `97.5%`)) + geom_errorbar() + 
  geom_hline(aes(yintercept = 0), color = "red") +
  theme(
    axis.text = element_text(angle = 45)
  ) + theme_bw() + 
  coord_flip() +
  xlab("Covariate") #+
#scale_y_continuous(breaks = seq(-5, 1, 1), limits = c(-5.5, 1))


ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor1_FISH_20241018_5factors.jpeg"), width = 5,
       height = 2,
       units = "in",
       dpi = 300)

# Factor 2
data_plot2 <- data_plot[((nrow(data_plot)/5)+1):(nrow(data_plot)/5+13),]

ggplot(data_plot2, aes(x = x,
                       ymin = `2.5%`,
                       ymax = `97.5%`)) + geom_errorbar() + 
  geom_hline(aes(yintercept = 0), color = "red") +
  theme(
    axis.text = element_text(angle = 45)
  ) + theme_bw() + 
  coord_flip() +
  xlab("Covariate") #+
#scale_y_continuous(breaks = seq(-5, 1, 1), limits = c(-5, 1))


ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor2_FISH_20241018_5factors.jpeg"), width = 5,
       height = 2,
       units = "in",
       dpi = 300)


# Factor 3
data_plot3 <- data_plot[((nrow(data_plot)/5)+14):(nrow(data_plot)/5+26),]

ggplot(data_plot3, aes(x = x,
                       ymin = `2.5%`,
                       ymax = `97.5%`)) + geom_errorbar() + 
  geom_hline(aes(yintercept = 0), color = "red") +
  theme(
    axis.text = element_text(angle = 45)
  ) + theme_bw() + 
  coord_flip() +
  xlab("Covariate") #+
#scale_y_continuous(breaks = seq(-5, 1, 1), limits = c(-5, 1))


ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor3_FISH_20241018_5factors.jpeg"), width = 5,
       height = 2,
       units = "in",
       dpi = 300)

# Factor 4
data_plot4 <- data_plot[((nrow(data_plot)/5)+27):(nrow(data_plot)/5+39),]

ggplot(data_plot4, aes(x = x,
                       ymin = `2.5%`,
                       ymax = `97.5%`)) + geom_errorbar() + 
  geom_hline(aes(yintercept = 0), color = "red") +
  theme(
    axis.text = element_text(angle = 45)
  ) + theme_bw() + 
  coord_flip() +
  xlab("Covariate") #+
#scale_y_continuous(breaks = seq(-5, 1, 1), limits = c(-5, 1))


ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor4_FISH_20241018_5factors.jpeg"), width = 5,
       height = 2,
       units = "in",
       dpi = 300)

# Factor 5
data_plot5 <- data_plot[((nrow(data_plot)/5)+40):((nrow(data_plot)/5)+52),]

ggplot(data_plot5, aes(x = x,
                       ymin = `2.5%`,
                       ymax = `97.5%`)) + geom_errorbar() + 
  geom_hline(aes(yintercept = 0), color = "red") +
  theme(
    axis.text = element_text(angle = 45)
  ) + theme_bw() + 
  coord_flip() +
  xlab("Covariate") #+
#scale_y_continuous(breaks = seq(-5, 1, 1), limits = c(-5, 1))


ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor5_FISH_20241018_5factors.jpeg"), width = 5,
       height = 2,
       units = "in",
       dpi = 300)




# FACTOR SCORES -------

param <- "U"

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
siteNames <- data$info$Site[!duplicated(data$info$Site)]


factor1 <- samples_subset[,1:101]

data_plot1 <- apply(factor1, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% 
  t %>% 
  as.data.frame %>% 
  mutate(site = siteNames) |>
  rename(F1.lower = `2.5%`, F1.upper = `97.5%`)

data_plot1$F1.mean <- apply(factor1, 2, mean)


coves <- read_tsv(here("data", "OTUtable_12S_toSP_GLG23_20240528.txt"))
coves <- coves %>% select(site, river) %>% distinct(site, river) %>% filter(!is.na(river))
data_plot1 <- left_join(data_plot1, coves, by = 'site')

ggplot(data_plot1, aes(x = site,
                       ymin = F1.lower,
                       ymax = F1.upper 
)) + geom_errorbar() + 
  #col = Regional_Division)) + geom_errorbar() + 
  coord_flip() + 
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0)
  ) + xlab("")

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor1scores_FISH_20241018_5factors.jpeg"), width = 8,
       height = 10,
       units = "in",
       dpi = 300)


factor2 <- samples_subset[,102:202]

data_plot2 <- apply(factor2, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% 
  t %>% 
  as.data.frame %>% 
  mutate(site = siteNames) |>
  rename(F2.lower = `2.5%`, F2.upper = `97.5%`)

data_plot2$F2.mean <- apply(factor2, 2, mean)

data_plot2 <- left_join(data_plot2, coves, by = 'site')

ggplot(data_plot2, aes(x = site,
                       ymin = F2.lower,
                       ymax = F2.upper
)) + geom_errorbar() + 
  # col = Regional_Division)) + geom_errorbar() +  
  coord_flip() + 
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0)
  ) + xlab("")

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor2scores_FISH_20241018_5factors.jpeg"), width = 8,
       height = 10,
       units = "in",
       dpi = 300)

factor3 <- samples_subset[,203:303]

data_plot3 <- apply(factor3, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% 
  t %>% 
  as.data.frame %>% 
  mutate(site = siteNames) |>
  rename(F3.lower = `2.5%`, F3.upper = `97.5%`)

data_plot3$F3.mean <- apply(factor3, 2, mean)

data_plot3 <- left_join(data_plot3, coves, by = 'site')

ggplot(data_plot3, aes(x = site,
                       ymin = F3.lower,
                       ymax = F3.upper
)) + geom_errorbar() + 
  # col = Regional_Division)) + geom_errorbar() +  
  coord_flip() + 
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0)
  ) + xlab("")

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor3scores_FISH_20241018_5factors.jpeg"), width = 8,
       height = 10,
       units = "in",
       dpi = 300)


factor4 <- samples_subset[,304:404]

data_plot4 <- apply(factor4, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% 
  t %>% 
  as.data.frame %>% 
  mutate(site = siteNames) |>
  rename(F4.lower = `2.5%`, F4.upper = `97.5%`)

data_plot4$F4.mean <- apply(factor4, 2, mean)

data_plot4 <- left_join(data_plot4, coves, by = 'site')

ggplot(data_plot4, aes(x = site,
                       ymin = F4.lower,
                       ymax = F4.upper
)) + geom_errorbar() + 
  # col = Regional_Division)) + geom_errorbar() +  
  coord_flip() + 
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0)
  ) + xlab("")

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor4scores_FISH_20241018_5factors.jpeg"), width = 8,
       height = 10,
       units = "in",
       dpi = 300)


factor5 <- samples_subset[,405:505]

data_plot5 <- apply(factor5, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% 
  t %>% 
  as.data.frame %>% 
  mutate(site = siteNames) |>
  rename(F5.lower = `2.5%`, F5.upper = `97.5%`)

data_plot5$F5.mean <- apply(factor5, 2, mean)

data_plot5 <- left_join(data_plot5, coves, by = 'site')

ggplot(data_plot5, aes(x = site,
                       ymin = F5.lower,
                       ymax = F5.upper
)) + geom_errorbar() + 
  # col = Regional_Division)) + geom_errorbar() +  
  coord_flip() + 
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0)
  ) + xlab("")

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor5scores_FISH_20241018_5factors.jpeg"), width = 8,
       height = 10,
       units = "in",
       dpi = 300)



# plot scores on a single plot
factorscores <- as.data.frame(cbind(data_plot1, data_plot2, data_plot3, data_plot4, data_plot5))[,-c(3,5,8,10,13,15,18,20)]
factorscores <- left_join(factorscores, data$info %>% select(Site, streamOrigin), by = c("site" = "Site")) %>% distinct()

# order Regional Divisions from high to low latitude
factorscores$river <- factor(factorscores$river, 
                             levels=c('nu', 'dl', 
                                      'lc', 'dy'))

# two-dimensional

# F1 ~ F2
ggplot(factorscores, aes(F1.mean, F2.mean, col=river)) +
  #ggplot(factorscores, aes(F1.mean, F2.mean, col=streamOrigin)) +
  geom_point() +
  geom_label_repel(aes(label = site),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  geom_errorbar(aes(xmin = F1.lower, xmax = F1.upper)) +
  geom_errorbar(aes(ymin = F2.lower, ymax = F2.upper)) + 
  stat_ellipse(level = 0.9, type = "norm") +
  xlab ("Factor 1 Scores") +
  ylab ("Factor 2 Scores") +
  theme_bw() +
  theme(
    legend.text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18)
  )
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores12_FISH_20241018_byRegion.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores12_FISH_20241018_bystreamOrigin.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)

# F1 ~ F3
#ggplot(factorscores, aes(F1.mean, F3.mean, col=river)) +
ggplot(factorscores, aes(F1.mean, F3.mean, col=streamOrigin)) +
  geom_point() +
  geom_label_repel(aes(label = site),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  geom_errorbar(aes(xmin = F1.lower, xmax = F1.upper)) +
  geom_errorbar(aes(ymin = F3.lower, ymax = F3.upper)) + 
  stat_ellipse(level = 0.9, type = "norm") +
  xlab ("Factor 1 Scores") +
  ylab ("Factor 3 Scores") +
  theme_bw() +
  theme(
    legend.text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18)
  )
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores13_FISH_20241018_byRegion.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores13_FISH_20241018_bystreamOrigin.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)

# F2 ~ F3
ggplot(factorscores, aes(F2.mean, F3.mean, col=river)) +
  #ggplot(factorscores, aes(F2.mean, F3.mean, col=streamOrigin)) +
  geom_point() +
  geom_label_repel(aes(label = site),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  geom_errorbar(aes(xmin = F2.lower, xmax = F2.upper)) +
  geom_errorbar(aes(ymin = F3.lower, ymax = F3.upper)) + 
  stat_ellipse(level = 0.9, type = "norm") +
  xlab ("Factor 2 Scores") +
  ylab ("Factor 3 Scores") +
  theme_bw() +
  theme(
    legend.text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18)
  )
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores23_FISH_20241018_byRegion.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores23_FISH_20241018_bystreamOrigin.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)


# three-dimensional
library(scatterplot3d)
library(car)

color_map <- c("nu" = "springgreen4", "dl" = "royalblue", 
               "lc" = "darkorange", "dy" = "violetred")
factorscores$color <- color_map[as.character(factorscores$river)]
color_map <- c("insideGLGS" = "tomato", "outsideGLGS" = "cyan3")
factorscores$color <- color_map[as.character(factorscores$streamOrigin)]

# F123

# Open a PNG file to save the plot
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores123_FISH_20241018_5factors_byRiver.jpeg"), width = 800, height = 800)
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores123_FISH_20241018_5factors_bystreamOrigin.jpeg"), width = 800, height = 800)

# Create the 3D scatter plot
s3d <- scatterplot3d(factorscores$F1.mean, factorscores$F2.mean, factorscores$F3.mean,
                     color = factorscores$color, pch = 16, cex.symbols = 2,
                     xlab = "Factor 1", ylab = "Factor 2", zlab = "Factor 3", 
                     angle = 45)

# Add error bars
with(factorscores, {
  # X-axis error bars
  arrows(s3d$xyz.convert(F1.lower, F2.mean, F3.mean)$x, s3d$xyz.convert(F1.lower, F2.mean, F3.mean)$y,
         s3d$xyz.convert(F1.upper, F2.mean, F3.mean)$x, s3d$xyz.convert(F1.upper, F2.mean, F3.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Y-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F2.lower, F3.mean)$x, s3d$xyz.convert(F1.mean, F2.lower, F3.mean)$y,
         s3d$xyz.convert(F1.mean, F2.upper, F3.mean)$x, s3d$xyz.convert(F1.mean, F2.upper, F3.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Z-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F2.mean, F3.lower)$x, s3d$xyz.convert(F1.mean, F2.mean, F3.lower)$y,
         s3d$xyz.convert(F1.mean, F2.mean, F3.upper)$x, s3d$xyz.convert(F1.mean, F2.mean, F3.upper)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
})

with(factorscores, {
  text(s3d$xyz.convert(F1.mean, F2.mean, F3.mean),
       labels = site,                             
       col = color,                               
       cex = 0.9,                                 
       pos = 4)                                   
})

# Add a legend based on region or streamOrigin
legend("topright", legend = names(color_map), col = color_map, pch = 16, cex = 1.5)

# Close the PNG device to save the file
dev.off()


# F124

# Open a PNG file to save the plot
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores124_FISH_20241018_5factors_byRegion.jpeg"), width = 800, height = 800)
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores124_FISH_20241018_5factors_bystreamOrigin.jpeg"), width = 800, height = 800)

# Create the 3D scatter plot
s3d <- scatterplot3d(factorscores$F1.mean, factorscores$F2.mean, factorscores$F4.mean,
                     color = factorscores$color, pch = 16, cex.symbols = 2,
                     xlab = "Factor 1", ylab = "Factor 2", zlab = "Factor 4", 
                     angle = 45)

# Add error bars
with(factorscores, {
  # X-axis error bars
  arrows(s3d$xyz.convert(F1.lower, F2.mean, F4.mean)$x, s3d$xyz.convert(F1.lower, F2.mean, F4.mean)$y,
         s3d$xyz.convert(F1.upper, F2.mean, F4.mean)$x, s3d$xyz.convert(F1.upper, F2.mean, F4.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Y-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F2.lower, F4.mean)$x, s3d$xyz.convert(F1.mean, F2.lower, F4.mean)$y,
         s3d$xyz.convert(F1.mean, F2.upper, F4.mean)$x, s3d$xyz.convert(F1.mean, F2.upper, F4.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Z-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F2.mean, F4.lower)$x, s3d$xyz.convert(F1.mean, F2.mean, F4.lower)$y,
         s3d$xyz.convert(F1.mean, F2.mean, F4.upper)$x, s3d$xyz.convert(F1.mean, F2.mean, F4.upper)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
})

with(factorscores, {
  text(s3d$xyz.convert(F1.mean, F2.mean, F4.mean),
       labels = site,                             
       col = color,                               
       cex = 0.9,                                 
       pos = 4)                                   
})

# Add a legend based on region or streamOrigin
legend("topright", legend = names(color_map), col = color_map, pch = 16, cex = 1.5)

# Close the PNG device to save the file
dev.off()


# F125

# Open a PNG file to save the plot
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores125_FISH_20241018_5factors_byRegion.jpeg"), width = 800, height = 800)
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores125_FISH_20241018_5factors_bystreamOrigin.jpeg"), width = 800, height = 800)

# Create the 3D scatter plot
s3d <- scatterplot3d(factorscores$F1.mean, factorscores$F2.mean, factorscores$F5.mean,
                     color = factorscores$color, pch = 16, cex.symbols = 2,
                     xlab = "Factor 1", ylab = "Factor 2", zlab = "Factor 5", 
                     angle = 45)

# Add error bars
with(factorscores, {
  # X-axis error bars
  arrows(s3d$xyz.convert(F1.lower, F2.mean, F5.mean)$x, s3d$xyz.convert(F1.lower, F2.mean, F5.mean)$y,
         s3d$xyz.convert(F1.upper, F2.mean, F5.mean)$x, s3d$xyz.convert(F1.upper, F2.mean, F5.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Y-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F2.lower, F5.mean)$x, s3d$xyz.convert(F1.mean, F2.lower, F5.mean)$y,
         s3d$xyz.convert(F1.mean, F2.upper, F5.mean)$x, s3d$xyz.convert(F1.mean, F2.upper, F5.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Z-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F2.mean, F5.lower)$x, s3d$xyz.convert(F1.mean, F2.mean, F5.lower)$y,
         s3d$xyz.convert(F1.mean, F2.mean, F5.upper)$x, s3d$xyz.convert(F1.mean, F2.mean, F5.upper)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
})

with(factorscores, {
  text(s3d$xyz.convert(F1.mean, F2.mean, F5.mean),
       labels = site,                             
       col = color,                               
       cex = 0.9,                                 
       pos = 4)                                   
})

# Add a legend based on region or streamOrigin
legend("topright", legend = names(color_map), col = color_map, pch = 16, cex = 1.5)

# Close the PNG device to save the file
dev.off()


# F234

# Open a PNG file to save the plot
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores234_FISH_20241018_5factors_byRegion.jpeg"), width = 800, height = 800)
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores234_FISH_20241018_5factors_bystreamOrigin.jpeg"), width = 800, height = 800)

# Create the 3D scatter plot
s3d <- scatterplot3d(factorscores$F2.mean, factorscores$F3.mean, factorscores$F4.mean,
                     color = factorscores$color, pch = 16, cex.symbols = 2,
                     xlab = "Factor 2", ylab = "Factor 3", zlab = "Factor 4", 
                     angle = 45)

# Add error bars
with(factorscores, {
  # X-axis error bars
  arrows(s3d$xyz.convert(F2.lower, F3.mean, F4.mean)$x, s3d$xyz.convert(F2.lower, F3.mean, F4.mean)$y,
         s3d$xyz.convert(F2.upper, F3.mean, F4.mean)$x, s3d$xyz.convert(F2.upper, F3.mean, F4.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Y-axis error bars
  arrows(s3d$xyz.convert(F2.mean, F3.lower, F4.mean)$x, s3d$xyz.convert(F2.mean, F3.lower, F4.mean)$y,
         s3d$xyz.convert(F2.mean, F3.upper, F4.mean)$x, s3d$xyz.convert(F2.mean, F3.upper, F4.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Z-axis error bars
  arrows(s3d$xyz.convert(F2.mean, F3.mean, F4.lower)$x, s3d$xyz.convert(F2.mean, F3.mean, F4.lower)$y,
         s3d$xyz.convert(F2.mean, F3.mean, F4.upper)$x, s3d$xyz.convert(F2.mean, F3.mean, F4.upper)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
})

with(factorscores, {
  text(s3d$xyz.convert(F2.mean, F3.mean, F4.mean),
       labels = site,                             
       col = color,                               
       cex = 0.9,                                 
       pos = 4)                                   
})

# Add a legend based on region or streamOrigin
legend("topright", legend = names(color_map), col = color_map, pch = 16, cex = 1.5)

# Close the PNG device to save the file
dev.off()


# F235

# Open a PNG file to save the plot
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores235_FISH_20241018_5factors_byRegion.jpeg"), width = 800, height = 800)
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores235_FISH_20241018_5factors_bystreamOrigin.jpeg"), width = 800, height = 800)

# Create the 3D scatter plot
s3d <- scatterplot3d(factorscores$F2.mean, factorscores$F3.mean, factorscores$F5.mean,
                     color = factorscores$color, pch = 16, cex.symbols = 2,
                     xlab = "Factor 2", ylab = "Factor 3", zlab = "Factor 5", 
                     angle = 45)

# Add error bars
with(factorscores, {
  # X-axis error bars
  arrows(s3d$xyz.convert(F2.lower, F3.mean, F5.mean)$x, s3d$xyz.convert(F2.lower, F3.mean, F5.mean)$y,
         s3d$xyz.convert(F2.upper, F3.mean, F5.mean)$x, s3d$xyz.convert(F2.upper, F3.mean, F5.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Y-axis error bars
  arrows(s3d$xyz.convert(F2.mean, F3.lower, F5.mean)$x, s3d$xyz.convert(F2.mean, F3.lower, F5.mean)$y,
         s3d$xyz.convert(F2.mean, F3.upper, F5.mean)$x, s3d$xyz.convert(F2.mean, F3.upper, F5.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Z-axis error bars
  arrows(s3d$xyz.convert(F2.mean, F3.mean, F5.lower)$x, s3d$xyz.convert(F2.mean, F3.mean, F5.lower)$y,
         s3d$xyz.convert(F2.mean, F3.mean, F5.upper)$x, s3d$xyz.convert(F2.mean, F3.mean, F5.upper)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
})

with(factorscores, {
  text(s3d$xyz.convert(F2.mean, F3.mean, F5.mean),
       labels = site,                             
       col = color,                               
       cex = 0.9,                                 
       pos = 4)                                   
})

# Add a legend based on region or streamOrigin
legend("topright", legend = names(color_map), col = color_map, pch = 16, cex = 1.5)

# Close the PNG device to save the file
dev.off()


# F134

# Open a PNG file to save the plot
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores134_FISH_20241018_5factors_byRegion.jpeg"), width = 800, height = 800)
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores134_FISH_20241018_5factors_bystreamOrigin.jpeg"), width = 800, height = 800)

# Create the 3D scatter plot
s3d <- scatterplot3d(factorscores$F1.mean, factorscores$F3.mean, factorscores$F4.mean,
                     color = factorscores$color, pch = 16, cex.symbols = 2,
                     xlab = "Factor 1", ylab = "Factor 3", zlab = "Factor 4", 
                     angle = 45)

# Add error bars
with(factorscores, {
  # X-axis error bars
  arrows(s3d$xyz.convert(F1.lower, F3.mean, F4.mean)$x, s3d$xyz.convert(F1.lower, F3.mean, F4.mean)$y,
         s3d$xyz.convert(F1.upper, F3.mean, F4.mean)$x, s3d$xyz.convert(F1.upper, F3.mean, F4.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Y-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F3.lower, F4.mean)$x, s3d$xyz.convert(F1.mean, F3.lower, F4.mean)$y,
         s3d$xyz.convert(F1.mean, F3.upper, F4.mean)$x, s3d$xyz.convert(F1.mean, F3.upper, F4.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Z-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F3.mean, F4.lower)$x, s3d$xyz.convert(F1.mean, F3.mean, F4.lower)$y,
         s3d$xyz.convert(F1.mean, F3.mean, F4.upper)$x, s3d$xyz.convert(F1.mean, F3.mean, F4.upper)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
})

with(factorscores, {
  text(s3d$xyz.convert(F1.mean, F3.mean, F4.mean),
       labels = site,                             
       col = color,                               
       cex = 0.9,                                 
       pos = 4)                                   
})

# Add a legend based on region or streamOrigin
legend("topright", legend = names(color_map), col = color_map, pch = 16, cex = 1.5)

# Close the PNG device to save the file
dev.off()


# F135

# Open a PNG file to save the plot
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores135_FISH_20241018_5factors_byRegion.jpeg"), width = 800, height = 800)
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores135_FISH_20241018_5factors_bystreamOrigin.jpeg"), width = 800, height = 800)

# Create the 3D scatter plot
s3d <- scatterplot3d(factorscores$F1.mean, factorscores$F3.mean, factorscores$F5.mean,
                     color = factorscores$color, pch = 16, cex.symbols = 2,
                     xlab = "Factor 1", ylab = "Factor 3", zlab = "Factor 5", 
                     angle = 45)

# Add error bars
with(factorscores, {
  # X-axis error bars
  arrows(s3d$xyz.convert(F1.lower, F3.mean, F5.mean)$x, s3d$xyz.convert(F1.lower, F3.mean, F5.mean)$y,
         s3d$xyz.convert(F1.upper, F3.mean, F5.mean)$x, s3d$xyz.convert(F1.upper, F3.mean, F5.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Y-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F3.lower, F5.mean)$x, s3d$xyz.convert(F1.mean, F3.lower, F5.mean)$y,
         s3d$xyz.convert(F1.mean, F3.upper, F5.mean)$x, s3d$xyz.convert(F1.mean, F3.upper, F5.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Z-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F3.mean, F5.lower)$x, s3d$xyz.convert(F1.mean, F3.mean, F5.lower)$y,
         s3d$xyz.convert(F1.mean, F3.mean, F5.upper)$x, s3d$xyz.convert(F1.mean, F3.mean, F5.upper)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
})

with(factorscores, {
  text(s3d$xyz.convert(F1.mean, F3.mean, F5.mean),
       labels = site,                             
       col = color,                               
       cex = 0.9,                                 
       pos = 4)                                   
})

# Add a legend based on region or streamOrigin
legend("topright", legend = names(color_map), col = color_map, pch = 16, cex = 1.5)

# Close the PNG device to save the file
dev.off()


# F345

# Open a PNG file to save the plot
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores345_FISH_20241018_5factors_byRegion.jpeg"), width = 800, height = 800)
png(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorscores345_FISH_20241018_5factors_bystreamOrigin.jpeg"), width = 800, height = 800)

# Create the 3D scatter plot
s3d <- scatterplot3d(factorscores$F3.mean, factorscores$F4.mean, factorscores$F5.mean,
                     color = factorscores$color, pch = 16, cex.symbols = 2,
                     xlab = "Factor 3", ylab = "Factor 4", zlab = "Factor 5", 
                     angle = 45)

# Add error bars
with(factorscores, {
  # X-axis error bars
  arrows(s3d$xyz.convert(F3.lower, F4.mean, F5.mean)$x, s3d$xyz.convert(F3.lower, F4.mean, F5.mean)$y,
         s3d$xyz.convert(F3.upper, F4.mean, F5.mean)$x, s3d$xyz.convert(F3.upper, F4.mean, F5.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Y-axis error bars
  arrows(s3d$xyz.convert(F3.mean, F4.lower, F5.mean)$x, s3d$xyz.convert(F3.mean, F4.lower, F5.mean)$y,
         s3d$xyz.convert(F3.mean, F4.upper, F5.mean)$x, s3d$xyz.convert(F3.mean, F4.upper, F5.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Z-axis error bars
  arrows(s3d$xyz.convert(F3.mean, F4.mean, F5.lower)$x, s3d$xyz.convert(F3.mean, F4.mean, F5.lower)$y,
         s3d$xyz.convert(F3.mean, F4.mean, F5.upper)$x, s3d$xyz.convert(F3.mean, F4.mean, F5.upper)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
})

with(factorscores, {
  text(s3d$xyz.convert(F3.mean, F4.mean, F5.mean),
       labels = site,                             
       col = color,                               
       cex = 0.9,                                 
       pos = 4)                                   
})

# Add a legend based on region or streamOrigin
legend("topright", legend = names(color_map), col = color_map, pch = 16, cex = 1.5)

# Close the PNG device to save the file
dev.off()



# FACTOR LOADINGS -----

param <- "L"

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
str(samples_subset)

factor1 <- samples_subset[,seq(1,ncol(samples_subset),5)]

data_plot1 <- apply(factor1, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% 
  t %>% 
  as.data.frame %>%
  mutate(species = spp) %>%
  rename(F1.lower = `2.5%`, F1.upper = `97.5%`)

data_plot1$F1.mean <- apply(factor1, 2, mean)

data_plot1 %>% 
  mutate(x = species) %>%
  #filter(`2.5%` > - 20) %>% 
  ggplot(aes(x = x,
             # y = trueParams,
             ymin = F1.lower,
             ymax = F1.upper)) + 
  geom_errorbar() + 
  coord_flip() + 
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0),
    axis.text.y = element_text(size = 7)
  )

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor1loadings_FISH_20241018_5factors.jpeg"), width = 6,
       height = 8,
       units = "in",
       dpi = 300)



factor2 <- samples_subset[,seq(2,ncol(samples_subset),5)]

data_plot2 <- apply(factor2, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% 
  t %>% 
  as.data.frame  %>%
  mutate(species = spp) |>
  rename(F2.lower = `2.5%`, F2.upper = `97.5%`)

data_plot2$F2.mean <- apply(factor2, 2, mean)

data_plot2 %>% 
  mutate(x = species) %>%
  #filter(`2.5%` > - 20) %>% 
  ggplot(aes(x = x,
             # y = trueParams,
             ymin = F2.lower,
             ymax = F2.upper)) + 
  geom_errorbar() + 
  coord_flip() + 
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0),
    axis.text.y = element_text(size = 7)
  )

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor2loadings_FISH_20241018_5factors.jpeg"), width = 6,
       height = 8,
       units = "in",
       dpi = 300)


factor3 <- samples_subset[,seq(3,ncol(samples_subset),5)]

data_plot3 <- apply(factor3, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% 
  t %>% 
  as.data.frame  %>%
  mutate(species = spp) |>
  rename(F3.lower = `2.5%`, F3.upper = `97.5%`)

data_plot3$F3.mean <- apply(factor3, 2, mean)

data_plot3 %>% 
  mutate(x = species) %>%
  #filter(`2.5%` > - 20) %>% 
  ggplot(aes(x = x,
             # y = trueParams,
             ymin = F3.lower,
             ymax = F3.upper)) + 
  geom_errorbar() + 
  coord_flip() + 
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0),
    axis.text.y = element_text(size = 7)
  )

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor3loadings_FISH_20241018_5factors.jpeg"), width = 6,
       height = 8,
       units = "in",
       dpi = 300)


factor4 <- samples_subset[,seq(4,ncol(samples_subset),5)]

data_plot4 <- apply(factor4, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% 
  t %>% 
  as.data.frame  %>%
  mutate(species = spp) |>
  rename(F4.lower = `2.5%`, F4.upper = `97.5%`)

data_plot4$F4.mean <- apply(factor4, 2, mean)

data_plot4 %>% 
  mutate(x = species) %>%
  #filter(`2.5%` > - 20) %>% 
  ggplot(aes(x = x,
             # y = trueParams,
             ymin = F4.lower,
             ymax = F4.upper)) + 
  geom_errorbar() + 
  coord_flip() + 
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0),
    axis.text.y = element_text(size = 7)
  )

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor4loadings_FISH_20241018_5factors.jpeg"), width = 6,
       height = 8,
       units = "in",
       dpi = 300)


factor5 <- samples_subset[,seq(5,ncol(samples_subset),5)]

data_plot5 <- apply(factor5, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% 
  t %>% 
  as.data.frame  %>%
  mutate(species = spp) |>
  rename(F5.lower = `2.5%`, F5.upper = `97.5%`)

data_plot5$F5.mean <- apply(factor5, 2, mean)

data_plot5 %>% 
  mutate(x = species) %>%
  #filter(`2.5%` > - 20) %>% 
  ggplot(aes(x = x,
             # y = trueParams,
             ymin = F5.lower,
             ymax = F5.upper)) + 
  geom_errorbar() + 
  coord_flip() + 
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0),
    axis.text.y = element_text(size = 7)
  )

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factor5loadings_FISH_20241018_5factors.jpeg"), width = 6,
       height = 8,
       units = "in",
       dpi = 300)



# plot loadings on a single plot

factorloadings <- as.data.frame(cbind(data_plot1, data_plot2))[,-3]

spp_list <- read_tsv(here("data","SP_INFO_20240528.txt"))
sp_status <- spp_list %>% select(species, status)
factorloadings <- left_join(factorloadings, sp_status, by = "species")

ggplot(factorloadings, aes(F1.mean, F2.mean, col=status)) +
  geom_point() +
  geom_errorbar(aes(ymin = F2.lower, ymax = F2.upper)) + 
  geom_errorbar(aes(xmin = F1.lower, xmax = F1.upper)) +
  xlab ("Factor 1 Loadings") +
  ylab ("Factor 2 Loadings") +
  theme_bw() +
  theme(
    legend.text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18)
  )

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "factorloadings12_FISH_20241018_5factors.jpeg"), width = 8,
       height = 5,
       units = "in",
       dpi = 300)


# FOR COMPARISON WITH MSOM COEFFICIENTS -----
Cov_names <- colnames(fitmodel$X_ord)
BL_output <- array(NA, dim = c(niter, ncov_ord, nspecies))

for(iter in 1:niter){
  B_iter <- beta_ord_output[iter,,]  # ncov_ord x d
  L_iter <- L_output[iter,,]         # d x nspecies
  BL_iter <- B_iter %*% L_iter
  BL_output[iter,,] <- BL_iter
}

niter <- dim(BL_output)[1]
ncov_ord <- dim(BL_output)[2]
nspecies <- dim(BL_output)[3]

effects_summary_BL <- data.frame()

for (covariate in 1:ncov_ord) {
  cov_effect <- BL_output[, covariate, ]
  cov_mean <- apply(cov_effect, 1, mean)
  cov_lci <- apply(cov_effect, 1, quantile, probs = 0.025)
  cov_uci <- apply(cov_effect, 1, quantile, probs = 0.975)
  effects_summary_BL <- rbind(effects_summary_BL, data.frame(
    variable = Cov_names[covariate], 
    mean = mean(cov_mean),
    lci = quantile(cov_lci, 0.025),
    uci = quantile(cov_uci, 0.975)
  ))
}

effects_long_BL <- effects_summary_BL %>%
  pivot_longer(cols = c(mean, lci, uci), names_to = ".value", values_to = "value")

ggplot(effects_long_BL, aes(x = variable, ymin = lci, ymax = uci)) +
  geom_errorbar() +
  #geom_point(aes(y = mean), size = 2, color = "blue") +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  theme_bw() +
  coord_flip() +
  labs(title = "occPlus - Posterior Distribution of Environmental Variable Effects",
       y = "Effect Size",
       x = "Covariate") +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  ) +
  scale_y_continuous(breaks = seq(-15, 2, 5), limits = c(-15, 4))

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "occPlus_Covariates_FISH_20241018.jpeg"), width = 5,
       height = 2,
       units = "in",
       dpi = 300)


effects_perSP_BL <- data.frame(species = character(),
                               covariate = character(),
                               mean_effect = numeric(),
                               ci_lower = numeric(),
                               ci_upper = numeric(),
                               stringsAsFactors = FALSE)
for (j in 1:dim(BL_output)[3]) {
  for (cov in 1:dim(BL_output)[2]) {
    effects <- BL_output[, cov, j]
    mean_effect <- mean(effects)
    ci_lower <- quantile(effects, 0.025)
    ci_upper <- quantile(effects, 0.975)
    effects_perSP_BL <- rbind(effects_perSP_BL, data.frame(
      species = spp[j],
      covariate = Cov_names[cov],
      mean_effect = mean_effect,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    ))
  }
}
write.csv(effects_perSP_BL, here("outputs","occPlus_outputs_Viorel covs","FISH_outputs_20241018_5factors","envCov_perSP_table_occPlus_FISH.csv"))


save(list = ls(), file = here("data","Rdata files","occPlus_FISH_20241018_5factors_Vioreldata.Rdata"))


