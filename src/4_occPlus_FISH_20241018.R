
library(here)
library(tidyverse)
library(stringr)
library(conflicted)
library("rstan") 
library(RColorBrewer)
library(ggrepel)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

conflicts_prefer(cowplot::align_plots)
conflicts_prefer(tibble::add_case)
conflicts_prefer(terra::allNA)
conflicts_prefer(terra::area)
conflicts_prefer(beepr::beep)
conflicts_prefer(stats::chisq.test)
conflicts_prefer(pracma::clear)
conflicts_prefer(dplyr::collapse)
conflicts_prefer(terra::compare)
conflicts_prefer(purrr::cross)
conflicts_prefer(janitor::crosstab)
conflicts_prefer(terra::delaunay)
conflicts_prefer(janitor::fisher.test)
conflicts_prefer(ggeffects::get_title)
conflicts_prefer(raster::getData)
conflicts_prefer(gstat::idw)
conflicts_prefer(pracma::integral)
conflicts_prefer(purrr::is_empty)
conflicts_prefer(lubridate::is.Date)
conflicts_prefer(terra::is.empty)
conflicts_prefer(dplyr::lag)
conflicts_prefer(scales::ordinal)
conflicts_prefer(pracma::quad)
conflicts_prefer(pracma::rat)
conflicts_prefer(arsenal::relrisk)
conflicts_prefer(janitor::remove_empty_cols)
conflicts_prefer(janitor::remove_empty_rows)
conflicts_prefer(tidyr::replace_na)
conflicts_prefer(terra::rotate)
conflicts_prefer(terra::shift)
conflicts_prefer(terra::size)
conflicts_prefer(terra::spin)
conflicts_prefer(lubridate::stamp)
conflicts_prefer(pracma::std)
conflicts_prefer(tictoc::tic)
conflicts_prefer(tictoc::toc)
conflicts_prefer(rstan::traceplot)
conflicts_prefer(terra::union)
conflicts_prefer(terra::varnames)
conflicts_prefer(terra::`varnames<-`)
conflicts_prefer(viridis::viridis_pal)
conflicts_prefer(terra::where.max)
conflicts_prefer(terra::where.min)
conflicts_prefer(tidyr::extract)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(base::intersect)
conflicts_prefer(foreign::read.dbf)
conflicts_prefer(foreign::write.dbf)


# ROOT FOLDER IS "Gaoligongsghan"

source(here("src", "5_occPlus_function_20241008.R"))

# data ---

load(here("data", "gaoligongshan_FISH_data_20240912.rda"))

data$info$streamOrigin <- as.factor(data$info$streamOrigin)
data$info$river <- as.factor(data$info$river)
data$info$Reservoir_downstream <- as.factor(data$info$Reservoir_downstream)
data$info$water_velocity[is.na(data$info$water_velocity)] <- mean(data$info$water_velocity, na.rm = T)

# run model ---
results <- runOccPlus(
  #data, d = 2,
  #data, d = 3,
  #data, d = 4,
  data, d = 5,
  #data, d = 6,
  #occCovariates = c("latitude", "longitude"),
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
  detCovariates = c("Filter_vol","PrevDayRain")
)

d2_results <- results
d3_results <- results
d4_results <- results
d5_results <- results
d6_results <- results

# comparison
# The loo comparison is used to select the optimal value of d
library(loo)

log_lik2 <- extract_log_lik(d2_results$vb_fit, parameter_name = "log_lik", merge_chains = FALSE)
log_lik3 <- extract_log_lik(d3_results$vb_fit, parameter_name = "log_lik", merge_chains = FALSE)
log_lik4 <- extract_log_lik(d4_results$vb_fit, parameter_name = "log_lik", merge_chains = FALSE)
log_lik5 <- extract_log_lik(d5_results$vb_fit, parameter_name = "log_lik", merge_chains = FALSE)
log_lik6 <- extract_log_lik(d6_results$vb_fit, parameter_name = "log_lik", merge_chains = FALSE)

loo2 <- loo::loo(log_lik2)
loo3 <- loo::loo(log_lik3)
loo4 <- loo::loo(log_lik4)
loo5 <- loo::loo(log_lik5)
loo6 <- loo::loo(log_lik6)

comparison <- loo_compare(loo2, loo3, loo4, loo5, loo6)
print(comparison)
#   elpd_diff se_diff
# model4    0.0       0.0 
# model3   -0.6      17.6 
# model5  -99.7      17.3 
# model2 -108.8      21.2 
# model1 -191.6      24.8 
detach("package:loo", unload = TRUE)

# analyze output ----
results <- d5_results

matrix_of_draws <- results$matrix_of_draws

nsites <- length(unique(data$info$Site))
nspecies <- ncol(data$OTU)
niter <- nrow(matrix_of_draws)
ncov_psi <- ncol(results$X_psi)
ncov_ord <- ncol(results$X_ord)
ncov_theta <- ncol(results$X_theta)

X_psi <- results$X_psi
X_ord <- results$X_ord
X_theta <- results$X_theta
d <- results$d

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
  speciesNames <- colnames(data$OTU)
  speciesNames <- make.unique(speciesNames)
  
  classesNames <- sapply(speciesNames, function(x){
    strsplit(x, split = "_")[[1]][1]
  })
}
##

siteNames <- data$info$Site[!duplicated(data$info$Site)]

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
species <- colnames(data$OTU)

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
  speciesNames <- colnames(data$OTU)
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
  # Note: In the theta_plot, if the dashed line (false positive probability) 
  # appears to the right of the solid line (true detection probability), 
  # this indicates a random bug occurred during model fitting. 
  # In such cases, re-run the model until the dashed and solid lines are in the correct order.
  
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



# FOR COMPARISON WITH MSOM COEFFICIENTS -----
Cov_names <- colnames(results$X_ord)
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


