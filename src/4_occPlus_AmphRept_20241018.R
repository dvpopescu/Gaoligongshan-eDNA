
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

load(here("data", "gaoligongshan_AmphRept_data_20241014.rda"))

data$info$streamOrigin <- as.factor(data$info$streamOrigin)
data$info$Regional_Division <- as.factor(data$info$Regional_Division)
data$info$river <- as.factor(data$info$river)

# run model ---
detach("package:loo", unload = TRUE)

results <- runOccPlus(
  #data, d = 2,
  data, d = 3,
  #data, d = 4,
  #data, d = 5,
  #data, d = 6,
  #occCovariates = c(),
  ordCovariates = c("annual_precipitation", 
                    "mean_temperature", 
                    "Percent_of_forest_cover_1km",
                    "Human_modification_1km",
                    "Patch_density_1km",
                    "Edge_density_1km",
                    "streamOrigin", 
                    "river",
                    "Regional_Division"
                    ),
  detCovariates = c("Filter_vol","PrevDayRain")
)

d2_results <- results
d3_results <- results
d4_results <- results
d5_results <- results
d6_results <- results

# comparison
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
# model2    0.0       0.0 
# model3  -19.8       9.3 
# model4  -83.0      10.8 
# model1  -94.8      16.6 
# model5 -155.3      14.1 

# analyze output ----
results <- d3_results

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

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "occupancy_occPlus_AmphRept_20241018_3factors.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)

write.csv(occ, here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "occupancy_occPlus_AmphRept_20241018_3factors.csv"))

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

write.csv(occ, here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "occupancybysite_occPlus_AmphRept_20241018_3factors.csv"))

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
  
  ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "1stagerate_AmphRept_20241018_3factors.jpeg"), width = 8,
         height = 6,
         units = "in",
         dpi = 300)
  
}

row.names(data_plot) <- NULL

write.csv(data_plot, here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "1stagerate_occPlus_AmphRept_20241018_3factors.csv"))

# 2 stage probabilities output
{
  
  params_CI_p <- apply(p_output, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>% t 
  
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
  
  ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "2stagerate_AmphRept_20241018_3factors.jpeg"), width = 8,
         height = 6,
         units = "in",
         dpi = 300)
  
}

row.names(data_plot) <- NULL

write.csv(data_plot, here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "2stagerate_occPlus_AmphRept_20241018_3factors.csv"))

# ORDINATION COVARIATE OUTPUT ------

param <- "beta_ord"

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws)),drop=F]

data_plot <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% 
  as.data.frame

data_plot$x <- colnames(results$X_ord)

# Factor 1
data_plot1 <- data_plot[1:(nrow(data_plot)/3),]

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


ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor1_AmphRept_20241018_3factors.jpeg"), width = 5,
       height = 2,
       units = "in",
       dpi = 300)

# Factor 2
data_plot2 <- data_plot[((nrow(data_plot)/3)+1):(nrow(data_plot)/3+12),]

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


ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor2_AmphRept_20241018_3factors.jpeg"), width = 5,
       height = 2,
       units = "in",
       dpi = 300)

# Factor 3
data_plot3 <- data_plot[((nrow(data_plot)/3)+13):(nrow(data_plot)/3+24),]

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


ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor3_AmphRept_20241018_3factors.jpeg"), width = 5,
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


zones <- read_tsv(here("data", "23GLG_zones_info_20240903.txt"))
data_plot1 <- left_join(data_plot1, zones, by = 'site')

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

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor1scores_AmphRept_20241018_3factors.jpeg"), width = 8,
       height = 10,
       units = "in",
       dpi = 300)

site_info <- data$info %>% select(Site, latitude, streamOrigin, river, altitude) %>% mutate(site = Site) %>% select(-Site)
data_plot1 <- left_join(data_plot1, site_info, by = 'site')
data_plot1 <- data_plot1 %>%
  mutate(site = factor(site, levels = unique(site[order(-altitude)])))
ggplot(data_plot1, aes(
  x = site,
  ymin = F1.lower,
  ymax = F1.upper,
  color = streamOrigin
)) +
  geom_errorbar(width = 0.2) +  
  coord_flip() +  
  theme_bw() +  
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 0, size = 12),  
    strip.text = element_text(size = 12, face = "bold")  
  ) +
  xlab("Site") + 
  ylab("F1") + 
  scale_color_manual(
    values = c("insideGLGS" = "tomato", "outsideGLGS" = "cyan3"),
    name = "Stream Origin"
  ) +
  facet_grid(rows = vars(Regional_Division), scales = "free_y", space = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor1scores_AmphRept_20241018_3factors_sorted.jpeg"), width = 5,
       height = 6,
       units = "in",
       dpi = 300)
ggplot(data_plot1, aes(
  x = site,
  ymin = F1.lower,
  ymax = F1.upper,
  color = river
)) +
  geom_errorbar(width = 0.2) +  
  coord_flip() +  
  theme_bw() +  
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 0, size = 12),  
    strip.text = element_text(size = 12, face = "bold")  
  ) +
  xlab("Site") + 
  ylab("F1 Estimate") + 
  facet_grid(rows = vars(Regional_Division), scales = "free_y")
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor1scores_AmphRept_20241018_3factors_sortedbyriver.jpeg"), width = 5,
       height = 6,
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

data_plot2 <- left_join(data_plot2, zones, by = 'site')

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

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor2scores_AmphRept_20241018_3factors.jpeg"), width = 8,
       height = 10,
       units = "in",
       dpi = 300)

site_info <- data$info %>% select(Site, latitude, streamOrigin, river, altitude) %>% mutate(site = Site) %>% select(-Site)
data_plot2 <- left_join(data_plot2, site_info, by = 'site')
data_plot2 <- data_plot2 %>%
  mutate(site = factor(site, levels = unique(site[order(-altitude)])))
ggplot(data_plot2, aes(
  x = site,
  ymin = F2.lower,
  ymax = F2.upper,
  color = streamOrigin
)) +
  geom_errorbar(width = 0.2) +  
  coord_flip() +  
  theme_bw() +  
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 0, size = 12),  
    strip.text = element_text(size = 12, face = "bold")  
  ) +
  xlab("Site") + 
  ylab("F2 Estimate") + 
  scale_color_manual(
    values = c("insideGLGS" = "tomato", "outsideGLGS" = "cyan3"),
    name = "Stream Origin"
  ) +
  facet_grid(rows = vars(Regional_Division), scales = "free_y", space = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor2scores_AmphRept_20241018_3factors_sorted.jpeg"), width = 5,
       height = 6,
       units = "in",
       dpi = 300)
ggplot(data_plot2, aes(
  x = site,
  ymin = F2.lower,
  ymax = F2.upper,
  color = river
)) +
  geom_errorbar(width = 0.2) +  
  coord_flip() +  
  theme_bw() +  
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 0, size = 12),  
    strip.text = element_text(size = 12, face = "bold")  
  ) +
  xlab("Site") + 
  ylab("F2 Estimate") + 
  facet_grid(rows = vars(Regional_Division), scales = "free_y")
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor2scores_AmphRept_20241018_3factors_sortedbyriver.jpeg"), width = 5,
       height = 6,
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

data_plot3 <- left_join(data_plot3, zones, by = 'site')

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

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor3scores_AmphRept_20241018_3factors.jpeg"), width = 8,
       height = 10,
       units = "in",
       dpi = 300)
site_info <- data$info %>% select(Site, latitude, streamOrigin, river) %>% mutate(site = Site) %>% select(-Site)
data_plot3 <- left_join(data_plot3, site_info, by = 'site')
data_plot3 <- data_plot3 %>%
  mutate(site = factor(site, levels = unique(site[order(-latitude)])))
ggplot(data_plot3, aes(
  x = site,
  ymin = F3.lower,
  ymax = F3.upper,
  color = streamOrigin
)) +
  geom_errorbar(width = 0.2) +  
  coord_flip() +  
  theme_bw() +  
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 0, size = 12),  
    strip.text = element_text(size = 12, face = "bold")  
  ) +
  xlab("Site") + 
  ylab("F3 Estimate") + 
  scale_color_manual(
    values = c("insideGLGS" = "tomato", "outsideGLGS" = "cyan3"),
    name = "Stream Origin"
  ) +
  facet_grid(rows = vars(Regional_Division), scales = "free_y")
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor3scores_AmphRept_20241018_3factors_sorted.jpeg"), width = 5,
       height = 6,
       units = "in",
       dpi = 300)
ggplot(data_plot3, aes(
  x = site,
  ymin = F3.lower,
  ymax = F3.upper,
  color = river
)) +
  geom_errorbar(width = 0.2) +  
  coord_flip() +  
  theme_bw() +  
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 0, size = 12),  
    strip.text = element_text(size = 12, face = "bold")  
  ) +
  xlab("Site") + 
  ylab("F3 Estimate") + 
  facet_grid(rows = vars(Regional_Division), scales = "free_y")
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor3scores_AmphRept_20241018_3factors_sortedbyriver.jpeg"), width = 5,
       height = 6,
       units = "in",
       dpi = 300)



# plot scores on a single plot
factorscores <- as.data.frame(cbind(data_plot1, data_plot2, data_plot3))[,-c(9,11,12,15,17,18)]
factorscores <- left_join(factorscores, data$info %>% select(Site, streamOrigin), by = c("site" = "Site")) %>% distinct()

# order Regional Divisions from high to low latitude
factorscores$Regional_Division <- factor(factorscores$Regional_Division, 
                                         levels=c('East', 
                                                  'Northwest',  
                                                  'Southwest'))


# two-dimensional

# F1 ~ F2
#ggplot(factorscores, aes(F1.mean, F2.mean, col=Regional_Division)) +
ggplot(factorscores, aes(F1.mean, F2.mean, col=streamOrigin)) +
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
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factorscores12_AmphRept_20241018_byRegion.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factorscores12_AmphRept_20241018_bystreamOrigin.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)

# F1 ~ F3
ggplot(factorscores, aes(F1.mean, F3.mean, col=Regional_Division)) +
#ggplot(factorscores, aes(F1.mean, F3.mean, col=streamOrigin)) +
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
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factorscores13_AmphRept_20241018_byRegion.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factorscores13_AmphRept_20241018_bystreamOrigin.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)

# F2 ~ F3
#ggplot(factorscores, aes(F2.mean, F3.mean, col=Regional_Division)) +
ggplot(factorscores, aes(F2.mean, F3.mean, col=streamOrigin)) +
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
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factorscores23_AmphRept_20241018_byRegion.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factorscores23_AmphRept_20241018_bystreamOrigin.jpeg"), width = 8,
       height = 6,
       units = "in",
       dpi = 300)


# three-dimensional
library(scatterplot3d)
library(car)

color_map <- c("East" = "darkorange", "Northwest" = "limegreen", "Southwest" = "dodgerblue1")
factorscores$color <- color_map[as.character(factorscores$Regional_Division)]
color_map <- c("insideGLGS" = "tomato", "outsideGLGS" = "cyan3")
factorscores$color <- color_map[as.character(factorscores$streamOrigin)]

# F123

# Open a PNG file to save the plot
png(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factorscores123_AmphRept_20241018_4factors_byRegion.jpeg"), width = 800, height = 800)
png(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factorscores123_AmphRept_20241018_4factors_bystreamOrigin.jpeg"), width = 800, height = 800)

# Create the 3D scatter plot
s3d <- scatterplot3d(factorscores$F1.mean, factorscores$F3.mean, factorscores$F2.mean,
                     color = factorscores$color, pch = 16, cex.symbols = 2,
                     xlab = "Factor 1", ylab = "Factor 3", zlab = "Factor 2", 
                     angle = 45)

# Add error bars
with(factorscores, {
  # X-axis error bars
  arrows(s3d$xyz.convert(F1.lower, F3.mean, F2.mean)$x, s3d$xyz.convert(F1.lower, F3.mean, F2.mean)$y,
         s3d$xyz.convert(F1.upper, F3.mean, F2.mean)$x, s3d$xyz.convert(F1.upper, F3.mean, F2.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Y-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F3.lower, F2.mean)$x, s3d$xyz.convert(F1.mean, F3.lower, F2.mean)$y,
         s3d$xyz.convert(F1.mean, F3.upper, F2.mean)$x, s3d$xyz.convert(F1.mean, F3.upper, F2.mean)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
  
  # Z-axis error bars
  arrows(s3d$xyz.convert(F1.mean, F3.mean, F2.lower)$x, s3d$xyz.convert(F1.mean, F3.mean, F2.lower)$y,
         s3d$xyz.convert(F1.mean, F3.mean, F2.upper)$x, s3d$xyz.convert(F1.mean, F3.mean, F2.upper)$y,
         col = color, angle = 90, length = 0.05, lwd = 1)
})

with(factorscores, {
  text(s3d$xyz.convert(F1.mean, F3.mean, F2.mean),
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

factor1 <- samples_subset[,seq(1,ncol(samples_subset),3)]

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

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor1loadings_AmphRept_20241018_3factors.jpeg"), width = 6,
       height = 5,
       units = "in",
       dpi = 300)


factor2 <- samples_subset[,seq(2,ncol(samples_subset),3)]

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

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor2loadings_AmphRept_20241018_3factors.jpeg"), width = 6,
       height = 5,
       units = "in",
       dpi = 300)

spp_list <- read_tsv(here("data","SP_withBODYSIZE.txt"))
spp_list <- spp_list %>% 
  mutate(IUCN = ifelse(is.na(IUCN), "NA", IUCN)) %>% 
  mutate(IUCN = ifelse(status == "dom" | status == "wild_or_dom", "domestic", IUCN))
spp_list <- spp_list %>% 
  mutate(China_redlist = ifelse(is.na(China_redlist), "NA", China_redlist)) %>% 
  mutate(China_redlist = ifelse(status == "dom" | status == "wild_or_dom", "domestic", China_redlist))
spp_list <- spp_list %>% mutate(ORDER = str_extract(otuID, "(?<=_)[^_]+(?=_)")) %>% 
  filter(class == "Amphibia" | class == "Reptilia") %>% select(species, ORDER, IUCN)

data_plot1 <- left_join(data_plot1, spp_list, by = "species")

IUCN_order <- c("CR", "EN", "VU", "NT", "LC", "DD", "NA")
IUCN_colors <- c("CR" = "brown", "EN" = "red", "VU" = "orange", "NT" = "dodgerblue", 
                 "LC" = "mediumseagreen", "DD" = "grey60", "NA" = "grey40")
order_labels <- c("Anura","Caudata","Squamata","Testudines")
dataplot_12 <- left_join(data_plot1, data_plot2, by = "species")
dataplot_long <- dataplot_12 %>%
  pivot_longer(
    cols = c(F1.lower, F1.upper, F1.mean, F2.lower, F2.upper, F2.mean),
    names_to = c("Type", ".value"),
    names_sep = "\\."
  )

ggplot(dataplot_long, aes(x = mean, y = species, color = IUCN)) +
  geom_point(size = 1) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +
  facet_grid(ORDER~Type, scales = "free_y", space = "free", labeller = labeller(ORDER = order_labels)) +
  labs(
    x = "Factor scores",
    y = NULL, 
    title = ""
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8, hjust = 1),
    strip.background = element_blank(),
    strip.placement = "outside") +
  scale_color_manual(values = IUCN_colors, breaks = IUCN_order) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")
ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor1&2loadings_AmphRept_20241018_3factors_sorted.jpeg"), width = 6,
       height = 8,
       units = "in",
       dpi = 300)






factor3 <- samples_subset[,seq(3,ncol(samples_subset),3)]

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

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factor3loadings_AmphRept_20241018_3factors.jpeg"), width = 6,
       height = 5,
       units = "in",
       dpi = 300)


# plot loadings on a single plot

factorloadings <- as.data.frame(cbind(data_plot1, data_plot2, data_plot3))[,-c(3, 7)]

factorloadings <- cbind(factorloadings, classesNames)

ggplot(factorloadings, aes(F1.mean, F2.mean, col=classesNames)) +
  geom_point() +
  geom_errorbar(aes(ymin = F2.lower, ymax = F2.upper)) + 
  geom_errorbar(aes(xmin = F1.lower, xmax = F1.upper)) +
  xlab ("Factor 1 Loadings") +
  ylab ("Factor 2 Loadings") +
  theme_bw()

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "factorloadings12_AmphRept_20241018_3factors.jpeg"), width = 8,
       height = 5,
       units = "in",
       dpi = 300)


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

ggsave(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "occPlus_Covariates_AmphRept_20241018.jpeg"), width = 5,
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
write.csv(effects_perSP_BL, here("outputs","occPlus_outputs_Viorel covs","AmphRept_outputs_20241018_3factors","envCov_perSP_table_occPlus_AmphRept.csv"))


save(list = ls(), file = here("data","Rdata files","occPlus_AmphRept_20241018_3factors_Vioreldata.Rdata"))


