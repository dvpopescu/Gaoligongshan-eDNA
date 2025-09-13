
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

source(here("src", "5_occPlus_function_singlespecies.R"))

# data ---

load(here("data", "gaoligongshan_bullfrog_data_20250127.rda"))

data$info$streamOrigin <- as.factor(data$info$streamOrigin)
data$info$Regional_Division <- as.factor(data$info$Regional_Division)
data$info$river <- as.factor(data$info$river)

# run model ---
results <- runOccPlus(
  data, d = 0,
  occCovariates = c("annual_precipitation", 
                    "mean_temperature", 
                    "Percent_of_forest_cover_1km",
                    "Human_modification_1km",
                    "streamOrigin",
                    "Patch_density_1km",
                    "Edge_density_1km",
                    "river",
                    "Regional_Division"
                    ),
  ordCovariates = c(),
  detCovariates = c("Filter_vol","PrevDayRain")
)

# analyze output ----
matrix_of_draws <- results$matrix_of_draws

nsites <- length(unique(data$info$Site))
nspecies <- 1
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
  matrix_of_draws[,"theta0"]
p_output <- 
  matrix_of_draws[,"p"]
q_output <- 
  matrix_of_draws[,"q"]

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


# OCCUPANCY PROBABILITY --------

# create classes names
{
  speciesNames <- names(data$OTU)
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
species <- names(data$OTU)

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

write.csv(occ, here("outputs", "output_bullfrog", "occupancy_occPlus_bullfrog_20250217.csv"))

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

write.csv(occ, here("outputs", "output_bullfrog", "occupancybysite_occPlus_bullfrog_20250217.csv"))

# DETECTION PROBABILITIES -------

# create classes names
{
  speciesNames <- names(data$OTU)
  speciesNames <- make.unique(speciesNames)
  
  classesNames <- sapply(speciesNames, function(x){
    strsplit(x, split = "_")[[1]][1]
  })
}

# 1 stage probabilities output

  
  beta_theta_output1 <- aperm(beta_theta_output, c(1,3,2))
  str(beta_theta_output1)
  
  params_CI_theta <- apply(beta_theta_output1, 2, function(x) {
    logistic(quantile(x, probs = c(0.025, 0.975)))
  }) %>% t 
  
  theta0_output <- matrix(theta0_output, ncol = 1)
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

row.names(data_plot) <- NULL

write.csv(data_plot, here("outputs", "output_bullfrog", "1stagerate_occPlus_bullfrog_20250217.csv"))

# 2 stage probabilities output

  p_output <- matrix(p_output, ncol = 1)
  params_CI_p <- apply(p_output, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>% t 
  
  q_output <- matrix(q_output, ncol = 1)
  params_CI_q <- apply(q_output, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>% t 
  
  params_mean_p <- apply(p_output, 2, mean)
  
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
  
row.names(data_plot) <- NULL

write.csv(data_plot, here("outputs", "output_bullfrog", "2stagerate_occPlus_bullfrog_20250217.csv"))


param <- "beta_psi"

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws)),drop=F]

data_plot <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% 
  as.data.frame

data_plot$x <- colnames(results$X_psi)

ggplot(data_plot, aes(x = x,
                       ymin = `2.5%`,
                       ymax = `97.5%`)) + geom_errorbar() + 
  geom_hline(aes(yintercept = 0), color = "red") +
  theme(
    axis.text = element_text(angle = 45)
  ) + theme_bw() + 
  coord_flip() +
  xlab("Covariate") #+
#scale_y_continuous(breaks = seq(-5, 1, 1), limits = c(-5.5, 1))
