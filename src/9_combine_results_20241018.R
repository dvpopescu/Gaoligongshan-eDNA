library(tidyverse)
library(here)
library(ggplot2)
library(ggrepel)
library(reshape2)

#######################################################################################################
# Plot detection and occupancy probabilities of MSOM and occPlus

# MSOM -------------------------------------------------------

# prepare data
# detection rates
msom_det_mammal <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Mammal_outputs_20241018", "detection_preds_msom_Mammal_20241018.csv"))
msom_det_ave <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Ave_outputs_20241018", "detection_preds_msom_Ave_20241018.csv"))
msom_det_herp <- read.csv(here("outputs", "msom_outputs_Viorel covs", "AmphRept_outputs_20241018", "detection_preds_msom_AmphRept_20241018.csv"))
msom_det_fish <- read.csv(here("outputs", "msom_outputs_Viorel covs", "FISH_outputs_20241018", "detection_preds_msom_FISH_20241018.csv"))

msom_det_combined <- bind_rows(msom_det_mammal, msom_det_ave, msom_det_herp, msom_det_fish) %>%
  mutate(SORT = factor(phylum, levels = c("Reptilia", "Aves", "Mammalia", "Amphibia", "Teleostei", "Actinopteri"))) %>%
  arrange(SORT, mean) %>% 
  mutate(species = factor(species, levels = unique(species))) %>% 
  mutate(CIrange = upper95CI - lower95CI)

# occupancy rates
msom_occ_mammal <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Mammal_outputs_20241018", "occupancy_preds_msom_Mammal_20241018.csv"))
msom_occ_ave <- read.csv(here("outputs", "msom_outputs_Viorel covs", "Ave_outputs_20241018", "occupancy_preds_msom_Ave_20241018.csv"))
msom_occ_herp <- read.csv(here("outputs", "msom_outputs_Viorel covs", "AmphRept_outputs_20241018", "occupancy_preds_msom_AmphRept_20241018.csv"))
msom_occ_fish <- read.csv(here("outputs", "msom_outputs_Viorel covs", "FISH_outputs_20241018", "occupancy_preds_msom_FISH_20241018.csv"))

msom_occ_combined <- bind_rows(msom_occ_mammal, msom_occ_ave, msom_occ_herp, msom_occ_fish) %>%
  mutate(SORT = factor(phylum, levels = c("Reptilia", "Aves", "Mammalia", "Amphibia", "Teleostei", "Actinopteri"))) %>%
  arrange(SORT, psi.mean.ms) %>% 
  mutate(species = factor(species, levels = unique(species))) %>% 
  mutate(CIrange = upperCI - lowerCI)


# occPlus -------------------------------------------------------

# prepare data
# 1st stage detection rates
OP_det1_mammal <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Mammal_outputs_20241018_4factors", "1stagerate_occPlus_Mammal_20241018_4factors.csv"))
OP_det1_ave <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Ave_outputs_20241018_2factors", "1stagerate_occPlus_Ave_20241018_2factors.csv"))
OP_det1_herp <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "1stagerate_occPlus_AmphRept_20241018_3factors.csv"))
OP_det1_fish <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "1stagerate_occPlus_FISH_20241018_5factors.csv"))

OP_det1_combined <- bind_rows(OP_det1_mammal, OP_det1_ave, OP_det1_herp, OP_det1_fish) %>%
  mutate(SORT = factor(phylum, levels = c("Reptilia", "Aves", "Mammalia", "Amphibia", "Teleostei", "Actinopteri"))) %>%
  arrange(SORT, mean1st) %>% 
  mutate(species = factor(species, levels = unique(species))) %>% 
  mutate(CIrange = theta.97.5. - theta.2.5.)

# 2 stage detection rates
OP_det2_mammal <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Mammal_outputs_20241018_4factors", "2stagerate_occPlus_Mammal_20241018_4factors.csv"))
OP_det2_ave <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Ave_outputs_20241018_2factors", "2stagerate_occPlus_Ave_20241018_2factors.csv"))
OP_det2_herp <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "2stagerate_occPlus_AmphRept_20241018_3factors.csv"))
OP_det2_fish <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "2stagerate_occPlus_FISH_20241018_5factors.csv"))

OP_det2_combined <- bind_rows(OP_det2_mammal, OP_det2_ave, OP_det2_herp, OP_det2_fish) %>%
  mutate(SORT = factor(phylum, levels = c("Reptilia", "Aves", "Mammalia", "Amphibia", "Teleostei", "Actinopteri"))) %>%
  arrange(SORT, mean2nd) %>% 
  mutate(species = factor(species, levels = unique(species))) %>% 
  mutate(CIrange = p.97.5. - p.2.5.)

# occPlus occupancy
OP_occ_mammal <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Mammal_outputs_20241018_4factors", "occupancy_occPlus_Mammal_20241018_4factors.csv"))
OP_occ_ave <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "Ave_outputs_20241018_2factors", "occupancy_occPlus_Ave_20241018_2factors.csv"))
OP_occ_herp <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "AmphRept_outputs_20241018_3factors", "occupancy_occPlus_AmphRept_20241018_3factors.csv"))
OP_occ_fish <- read.csv(here("outputs", "occPlus_outputs_Viorel covs", "FISH_outputs_20241018_5factors", "occupancy_occPlus_FISH_20241018_5factors.csv"))

OP_occ_combined <- bind_rows(OP_occ_mammal, OP_occ_ave, OP_occ_herp, OP_occ_fish) %>%
  mutate(SORT = factor(phylum, levels = c("Reptilia", "Aves", "Mammalia", "Amphibia", "Teleostei", "Actinopteri"))) %>%
  arrange(SORT, meanocc) %>% 
  mutate(species = factor(species, levels = unique(species))) %>% 
  mutate(CIrange = X97.5. - X2.5.)


# plotting for each species -----------------------------------------------------------

# msom detection
ggplot(msom_det_combined) +  
  geom_errorbar(aes(x = species,
                    ymin = lower95CI,
                    ymax = upper95CI,
                    color = phylum)) +
  geom_point(aes(x = species,
                 y = mean,
                 color = phylum)) +
  scale_color_manual(values = c(
    "Reptilia" = "black",
    "Aves" = 'deepskyblue4',
    "Mammalia" = 'palevioletred4',
    "Amphibia" = 'magenta3',
    "Teleostei" = 'blue',
    "Actinopteri" = 'grey'
  ),
  breaks = c("Actinopteri", "Teleostei", "Amphibia", "Mammalia", "Aves", "Reptilia")) +
  coord_flip() + 
  ylab("MSOM Detection Probability") + 
  xlab("Species") +
  theme(
    legend.title = element_text(size = 28), 
    legend.text = element_text(size = 24)   
  )

ggsave(here("outputs","combine_results","detection_preds_msom_20241018.jpeg"), width = 8,
       height = 30,
       units = "in",
       dpi = 300)

# msom occupancy
ggplot(msom_occ_combined) +  
  geom_errorbar(aes(x = species,
                    ymin = lowerCI,
                    ymax = upperCI,
                    color = phylum)) +
  geom_point(aes(x = species,
                 y = psi.mean.ms,
                 color = phylum)) +
  scale_color_manual(values = c(
    "Reptilia" = "black",
    "Aves" = 'deepskyblue4',
    "Mammalia" = 'palevioletred4',
    "Amphibia" = 'magenta3',
    "Teleostei" = 'blue',
    "Actinopteri" = 'grey'
  ),
  breaks = c("Actinopteri", "Teleostei", "Amphibia", "Mammalia", "Aves", "Reptilia")) +
  coord_flip() + 
  ylab("MSOM Occupancy Probability") + 
  xlab("Species") +
  theme(
    legend.title = element_text(size = 28), 
    legend.text = element_text(size = 24)   
  )

ggsave(here("outputs","combine_results","occupancy_msom_20241018.jpeg"), width = 8,
       height = 30,
       units = "in",
       dpi = 300)

# occPlus 1 stage detection
ggplot(OP_det1_combined) +  
  geom_errorbar(aes(x = species,
                    ymin = theta.2.5.,
                    ymax = theta.97.5.,
                    color = phylum,
                    linetype = "True positive rate")) + 
  geom_errorbar(aes(x = species,
                    ymin = theta0.2.5.,
                    ymax = theta0.97.5.,
                    color = phylum,
                    linetype = "False positive rate")) + 
  geom_point(aes(x = species,
                 y = mean1st,
                 color = phylum)) +
  scale_color_manual(values = c(
    "Reptilia" = "black",
    "Aves" = 'deepskyblue4',
    "Mammalia" = 'palevioletred4',
    "Amphibia" = 'magenta3',
    "Teleostei" = 'blue',
    "Actinopteri" = 'grey'
  ),
  breaks = c("Actinopteri", "Teleostei", "Amphibia", "Mammalia", "Aves", "Reptilia")) +
  scale_linetype_manual(
    name = "Legend", 
    values = c("True positive rate" = "solid", "False positive rate" = "dashed"), 
    labels = c("False positive rate", "True positive rate") 
  ) +
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0,
                             size = 8)
  ) + coord_flip() + ylab("Stage 1 Detection Probability") + xlab("Species") +
  theme(
    legend.title = element_text(size = 28), 
    legend.text = element_text(size = 24)   
  )


ggsave(here("outputs","combine_results","stage1_detection_occ_20241018.jpeg"), width = 8,
       height = 30,
       units = "in",
       dpi = 300)

# occPlus 2 stage detection
ggplot(OP_det2_combined) +  
  geom_errorbar(aes(x = species,
                    ymin = p.2.5.,
                    ymax = p.97.5.,
                    color = phylum,
                    linetype = "True positive rate")) + 
  geom_errorbar(aes(x = species,
                    ymin = q.2.5.,
                    ymax = q.97.5.,
                    color = phylum,
                    linetype = "False positive rate")) + 
  geom_point(aes(x = species,
                 y = mean2nd,
                 color = phylum)) +
  scale_color_manual(values = c(
    "Reptilia" = "black",
    "Aves" = 'deepskyblue4',
    "Mammalia" = 'palevioletred4',
    "Amphibia" = 'magenta3',
    "Teleostei" = 'blue',
    "Actinopteri" = 'grey'
  ),
  breaks = c("Actinopteri", "Teleostei", "Amphibia", "Mammalia", "Aves", "Reptilia")) +
  scale_linetype_manual(
    name = "Legend", 
    values = c("True positive rate" = "solid", "False positive rate" = "dashed"), 
    labels = c("False positive rate", "True positive rate") 
  ) +
  theme_bw() + 
  theme(
    axis.text = element_text(angle = 0,
                             size = 8)
  ) + coord_flip() + ylab("Stage 2 Detection Probability") + xlab("Species") +
  theme(
    legend.title = element_text(size = 28), 
    legend.text = element_text(size = 24)   
  )


ggsave(here("outputs","combine_results","stage2_detection_occ_20241018.jpeg"), width = 8,
       height = 30,
       units = "in",
       dpi = 300)

# occPlus occupancy
ggplot(OP_occ_combined, aes(x = species, y = meanocc, col=phylum)) +  
  geom_point() + 
  geom_errorbar(aes(x = species,
                    ymin = `X2.5.`,
                    ymax = `X97.5.`,
                    color = phylum)) +
  scale_color_manual(values = c(
    "Reptilia" = "black",
    "Aves" = 'deepskyblue4',
    "Mammalia" = 'palevioletred4',
    "Amphibia" = 'magenta3',
    "Teleostei" = 'blue',
    "Actinopteri" = 'grey'
  ),
  breaks = c("Actinopteri", "Teleostei", "Amphibia", "Mammalia", "Aves", "Reptilia")) +
  coord_flip() + 
  ylab("OccPlus Occupancy Probability") + 
  xlab("Species") +
  theme(
    legend.title = element_text(size = 28), 
    legend.text = element_text(size = 24)   
  )


ggsave(here("outputs","combine_results","occupancy_occPlus_20241018.jpeg"), width = 8,
       height = 30,
       units = "in",
       dpi = 300)


#######################################################################################################
# Drawing a scatterplot to visualize the correlation between the predicted probabilities of MSOM and OccPlus

# Load the main species list once
spp_list <- read_tsv(here("data","SP_withBODYSIZE.txt"))

#' Create a Scatterplot of MSOM vs. OccPlus Probabilities
#'
#' This function reads MSOM and OccPlus data for one or more taxonomic groups,
#' merges them, and generates a scatterplot comparing their predicted probabilities.
#'
#' @param msom_files A named list where names are groups (e.g., "AmphRept") and values are file paths to MSOM CSVs.
#' @param occplus_files A named list with the same structure for OccPlus CSVs.
#' @param color_palette A named vector for ggplot's scale_color_manual.
#' @param output_filename The name for the saved JPEG file.
#' @return A ggplot object.
create_correlation_scatterplot <- function(msom_files, occplus_files, color_palette, output_filename) {
  
  # Read and combine all MSOM files
  msom_data_list <- lapply(msom_files, function(file) read.csv(here(file)))
  # For multiple files, we need to join them by site
  msom_combined <- reduce(msom_data_list, left_join, by = "site")
  
  # Read and combine all OccPlus files
  occplus_data_list <- lapply(occplus_files, function(file) read.csv(here(file)))
  occplus_combined <- reduce(occplus_data_list, left_join, by = "site")
  
  # Convert to long format
  msom_long <- msom_combined %>% 
    dplyr::select(site, starts_with("Mean.")) %>% 
    pivot_longer(cols = starts_with("Mean."), names_to = "species", values_to = "msom_prob")
  
  occplus_long <- occplus_combined %>% 
    dplyr::select(site, starts_with("Mean")) %>% 
    pivot_longer(cols = starts_with("Mean"), names_to = "species", values_to = "occ_prob")
  
  # Merge the two long-format datasets and join with species info
  merged_data <- inner_join(msom_long, occplus_long, by = c("site", "species")) %>% 
    mutate(species = gsub("^Mean.", "", species), 
           species = gsub("\\.", " ", species)) %>% 
    left_join(spp_list, by = "species")
  
  # Create the plot
  p <- ggplot(merged_data, aes(msom_prob, occ_prob, color = class)) +
    geom_point(alpha = 0.6) + # Added alpha for better visualization of dense points
    scale_color_manual(values = color_palette) +
    xlab("MSOM predicted occupancy probability") +
    ylab("OccPlus predicted occupancy probability") +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") + 
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) + # Use coord_fixed for a square plot
    theme_bw() + # A clean theme
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
  
  # Save the plot
  ggsave(here("outputs", "combine_results", output_filename), 
         plot = p, width = 10, height = 8, units = "in", dpi = 300)
  
  return(p)
}

# 1. Generate the plot for FISH
create_correlation_scatterplot(
  msom_files = c(FISH = "outputs/msom_outputs_Viorel covs/FISH_outputs_20241018/occupancy_bysite_msom_FISH_20241018.csv"),
  occplus_files = c(FISH = "outputs/occPlus_outputs_Viorel covs/FISH_outputs_20241018_5factors/occupancybysite_occPlus_FISH_20241018_5factors.csv"),
  color_palette = c("Teleostei" = 'blue', "Actinopteri" = 'grey'),
  output_filename = "predicted_prob_MSOM_vs_OCC_bySPbySITE_FISH_20241018.jpeg"
)

# 2. Generate the plot for TERRESTRIAL vertebrates
create_correlation_scatterplot(
  msom_files = c(
    AmphRept = "outputs/msom_outputs_Viorel covs/AmphRept_outputs_20241018/occupancy_bysite_msom_AmphRept_20241018.csv",
    Mammal = "outputs/msom_outputs_Viorel covs/Mammal_outputs_20241018/occupancy_bysite_msom_Mammal_20241018.csv",
    Aves = "outputs/msom_outputs_Viorel covs/Ave_outputs_20241018/occupancy_bysite_msom_Ave_20241018.csv"
  ),
  occplus_files = c(
    AmphRept = "outputs/occPlus_outputs_Viorel covs/AmphRept_outputs_20241018_3factors/occupancybysite_occPlus_AmphRept_20241018_3factors.csv",
    Mammal = "outputs/occPlus_outputs_Viorel covs/Mammal_outputs_20241018_4factors/occupancybysite_occPlus_Mammal_20241018_4factors.csv",
    Aves = "outputs/occPlus_outputs_Viorel covs/Ave_outputs_20241018_2factors/occupancybysite_occPlus_Ave_20241018_2factors.csv"
  ),
  color_palette = c("Reptilia" = "black", "Aves" = 'deepskyblue4', "Mammalia" = 'palevioletred4', "Amphibia" = 'magenta3'),
  output_filename = "predicted_prob_MSOM_vs_OCC_bySPbySITE_TERR_20241018.jpeg"
)



#######################################################################################################
# FOCUSED ANALYSIS AND PLOTTING

# --- 1. DATA PREPARATION (Run Once) ---

# Load and prepare the main species information list
spp_list <- read_tsv(here("data","SP_withBODYSIZE.txt")) %>%
  mutate(
    # Create a more robust 'status' column specifically for plotting
    plot_status = case_when(
      status %in% c("dom", "wild_or_dom") ~ "domestic",
      class == "Teleostei" & status == "invasive" ~ "invasive",
      class == "Teleostei" & status != "invasive" ~ "native", #<-- CRITICAL FIX!
      TRUE ~ NA_character_ # All others are NA for this specific context
    ),
    
    # Clean IUCN and China_redlist
    IUCN = ifelse(is.na(IUCN) | status %in% c("dom", "wild_or_dom"), "domestic", IUCN),
    China_redlist = ifelse(is.na(China_redlist) | status %in% c("dom", "wild_or_dom"), "domestic", China_redlist)
  )

# --- 2. FUNCTION DEFINITIONS ---

#' Prepare Occupancy Data for Inside vs. Outside Comparison
#'
#' This function takes model output, processes it to long format, aggregates by river origin,
#' and prepares it for scatter plotting.
#'
#' @param model_occ_df A dataframe with model occupancy predictions.
#' @return A processed dataframe ready for plotting.
prepare_origin_comparison_data <- function(model_occ_df) {
  site_info <- model_occ_df %>% dplyr::select(site, riverOrigin)
  
  # Process and pivot data to long format
  model_occ_long <- model_occ_df %>%
    dplyr::select(site, starts_with("Mean.")) %>% 
    pivot_longer(cols = starts_with("Mean."), names_to = "species", values_to = "Mean") %>%
    mutate(species = gsub("^Mean\\.", "", species),
           species = gsub("\\.", " ", species)) %>%   
    left_join(site_info, by = "site")
  
  # This entire block is now one continuous pipe.
  model_occ_wide <- model_occ_long %>% 
    group_by(species, riverOrigin) %>% 
    summarise(Mean_prob = mean(Mean), .groups = "drop") %>%
    pivot_wider(names_from = riverOrigin, values_from = Mean_prob) %>%
    rename(Mean.in = insideGLGS, Mean.out = outsideGLGS) %>%
    left_join(spp_list, by = "species")
  
  return(model_occ_wide)
}

#' Create an Inside vs. Outside Occupancy Scatterplot
#'
#' Generates a scatterplot comparing occupancy inside vs. outside PAs.
#'
#' @param data The processed dataframe from `prepare_origin_comparison_data`.
#' @param color_var The column to use for coloring points (e.g., "IUCN", "status").
#' @param color_palette A named vector of colors.
#' @param color_limits A vector specifying the order for the color legend.
#' @param color_labels A named vector for the color legend labels.
#' @param legend_title The title for the color legend (e.g., "IUCN Status", "Status").
#' @param is_mammal_bodysize A boolean to indicate if it's the special mammal bodysize plot.
#' @param output_filename The name for the saved JPEG file.
create_origin_scatterplot <- function(data, color_var, color_palette, color_limits, color_labels, legend_title, is_mammal_bodysize = FALSE, output_filename) {
  
  # Base plot
  p <- ggplot(data, aes(y = Mean.in, x = Mean.out, color = !!sym(color_var))) +
    geom_abline(intercept = 0, slope = 1) +
    xlab("Occupancy probability — river origin outside") +
    ylab("Occupancy probability — river origin inside") +
    scale_color_manual(name = legend_title, values = color_palette, limits = color_limits, labels = color_labels) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18)
    ) + 
    guides(color = guide_legend(override.aes = list(size = 5)))
  # Add specific layers for different plot types
  if (is_mammal_bodysize) {
    data <- data %>%
      mutate(sized_class = cut(AdultBodyMass, 
                               breaks = c(-Inf, 500, 5000, Inf), 
                               labels = c("< 0.5 kg", "0.5 kg to 5 kg", "> 5 kg")))
    
    p <- p +
      geom_point(data = data, aes(size = sized_class)) +
      scale_size_manual(values = c("< 0.5 kg" = 3, "0.5 kg to 5 kg" = 5, "> 5 kg" = 7), name = "Adult Body Mass") +
      geom_text_repel(aes(label = species), size = 3, box.padding = 0.5, max.overlaps = 10)
    width <- 12
    height <- 10
  } else {
    p <- p + geom_point(size = 3)
    width <- 8
    height <- 6
  }
  
  ggsave(here("outputs", "combine_results", "riverOrigin", output_filename), 
         plot = p, width = width, height = height, units = "in", dpi = 300)
  
  return(p)
}


# --- 3. ANALYSIS AND PLOTTING ---

# Define color palettes and orders once
IUCN_colors_mammal <- c("EN" = "red", "VU" = "orange", "NT" = "dodgerblue", "LC" = "mediumseagreen", "DD" = "grey80", "NA" = "grey60", "domestic" = "black")
IUCN_order_mammal <- c("EN", "VU", "NT", "LC", "DD", "NA", "domestic")
IUCN_labels_mammal <- c("Endangered", "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient", "Not Applicable", "Domestic")

IUCN_colors_ave <- c("VU" = "orange", "NT" = "dodgerblue", "LC" = "mediumseagreen", "domestic" = "black")
IUCN_order_ave <- c("VU", "NT", "LC", "domestic")
IUCN_labels_ave <- c("Vulnerable", "Near Threatened", "Least Concern", "Domestic")

IUCN_colors_herp <- c("CR" = "brown", "EN" = "red", "VU" = "orange", "NT" = "dodgerblue", "LC" = "mediumseagreen", "DD" = "grey60", "NA" = "grey40")
IUCN_order_herp <- c("CR", "EN", "VU", "NT", "LC", "DD", "NA")
IUCN_labels_herp <- c("Critically Endangered", "Endangered", "Vulnerable", "Near Threatened", "Least Concern", "Data Deficient", "Not Applicable")

colors_fish <- c("native" = "royalblue1", "invasive" = "palevioletred")
order_fish <- c("native", "invasive")
labels_fish <- c("Native", "Non-native")


## --- Generate scatterplots for each group ---
groups_to_plot <- list(
  Mammal_bodysize = list(
    msom_file = "outputs/msom_outputs_Viorel covs/Mammal_outputs_20241018/occupancy_bysite_msom_Mammal_20241018.csv",
    occplus_file = "outputs/occPlus_outputs_Viorel covs/Mammal_outputs_20241018_4factors/occupancybysite_occPlus_Mammal_20241018_4factors.csv",
    color_var = "IUCN", 
    palette = IUCN_colors_mammal, 
    limits = IUCN_order_mammal, 
    labels = IUCN_labels_mammal, 
    legend_title = "IUCN Status", 
    is_bodysize = TRUE
  ),
  Ave = list(
    msom_file = "outputs/msom_outputs_Viorel covs/Ave_outputs_20241018/occupancy_bysite_msom_Ave_20241018.csv",
    occplus_file = "outputs/occPlus_outputs_Viorel covs/Ave_outputs_20241018_2factors/occupancybysite_occPlus_Ave_20241018_2factors.csv",
    color_var = "IUCN", 
    palette = IUCN_colors_ave, 
    limits = IUCN_order_ave, 
    labels = IUCN_labels_ave, 
    legend_title = "IUCN Status", 
    is_bodysize = FALSE
  ),
  AmphRept = list(
    msom_file = "outputs/msom_outputs_Viorel covs/AmphRept_outputs_20241018/occupancy_bysite_msom_AmphRept_20241018.csv",
    occplus_file = "outputs/occPlus_outputs_Viorel covs/AmphRept_outputs_20241018_3factors/occupancybysite_occPlus_AmphRept_20241018_3factors.csv",
    color_var = "IUCN", 
    palette = IUCN_colors_herp, 
    limits = IUCN_order_herp, 
    labels = IUCN_labels_herp, 
    legend_title = "IUCN Status", 
    is_bodysize = FALSE
  ),
  FISH = list(
    msom_file = "outputs/msom_outputs_Viorel covs/FISH_outputs_20241018/occupancy_bysite_msom_FISH_20241018.csv",
    occplus_file = "outputs/occPlus_outputs_Viorel covs/FISH_outputs_20241018_5factors/occupancybysite_occPlus_FISH_20241018_5factors.csv",
    color_var = "plot_status", 
    palette = colors_fish, 
    limits = order_fish, 
    labels = labels_fish, 
    legend_title = "Status", 
    is_bodysize = FALSE
  )
)

for (group_name in names(groups_to_plot)) {
  params <- groups_to_plot[[group_name]]
  
  # Process MSOM data
  msom_df <- read.csv(here(params$msom_file))
  msom_processed_data <- prepare_origin_comparison_data(msom_df)
  
  # Process OccPlus data
  occplus_df <- read.csv(here(params$occplus_file)) %>% inner_join(msom_df[,c("site", "riverOrigin")], by = "site")
  occplus_processed_data <- prepare_origin_comparison_data(occplus_df)
  
  # Create MSOM plot
  group_label <- ifelse(group_name == "Mammal_bodysize", "Mammal", group_name)
  create_origin_scatterplot(
    data = msom_processed_data,
    color_var = params$color_var,
    color_palette = params$palette,
    color_limits = params$limits,
    color_labels = params$labels,
    legend_title = params$legend_title, # <-- Passing the new parameter
    is_mammal_bodysize = params$is_bodysize,
    output_filename = paste0("scatterplot_occupancy_", group_label, "_MSOM_20241018_groupbyOrigin", ifelse(params$is_bodysize, "_bodysize", ""), ".jpeg")
  )
  
  # Create OccPlus plot
  create_origin_scatterplot(
    data = occplus_processed_data,
    color_var = params$color_var,
    color_palette = params$palette,
    color_limits = params$limits,
    color_labels = params$labels,
    legend_title = params$legend_title, # <-- Passing the new parameter
    is_mammal_bodysize = params$is_bodysize,
    output_filename = paste0("scatterplot_occupancy_", group_label, "_OccPlus_20241018_groupbyOrigin", ifelse(params$is_bodysize, "_bodysize", ""), ".jpeg")
  )
}



#######################################################################################################
# regressions of altitude and latitude by occupancy corrected species richness

glgs <- read_tsv(here("data", "OTUtable_12S_toSP_GLG23_20240528.txt"))

## --- Generate Richness vs. Covariate Plots ---
glgs_covs <- glgs %>% 
  dplyr::select(site, longitude:DistFromPAedge) %>% 
  filter(status == "sample") %>% 
  summarise(across(longitude:DistFromPAedge, first), .by = c(site)) %>% 
  mutate(riverOrigin = if_else(is.na(DistFromPAedge), "outsideGLGS", "insideGLGS"))

occplus_files_all <- c(
  "outputs/occPlus_outputs_Viorel covs/AmphRept_outputs_20241018_3factors/occupancybysite_occPlus_AmphRept_20241018_3factors.csv",
  "outputs/occPlus_outputs_Viorel covs/Mammal_outputs_20241018_4factors/occupancybysite_occPlus_Mammal_20241018_4factors.csv",
  "outputs/occPlus_outputs_Viorel covs/Ave_outputs_20241018_2factors/occupancybysite_occPlus_Ave_20241018_2factors.csv",
  "outputs/occPlus_outputs_Viorel covs/FISH_outputs_20241018_5factors/occupancybysite_occPlus_FISH_20241018_5factors.csv"
)

occplus_all_data <- occplus_files_all %>%
  lapply(function(file) read.csv(here(file))) %>%
  reduce(left_join, by = "site") %>%
  left_join(glgs_covs, by = "site") %>%
  mutate(summed_probs = rowSums(across(starts_with("Mean."))))

# Plot by latitude
lat_plot <- ggplot(occplus_all_data, aes(latitude, summed_probs)) +
  geom_point(aes(color = riverOrigin)) +
  geom_smooth(method = "gam", alpha = .2, aes(fill = riverOrigin)) +
  xlab(paste0("S", strrep(" ", 20), "Latitude (ºN)", strrep(" ", 20), "N")) +
  ylab("Summed species occupancies") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )
ggsave(here("outputs", "combine_results", "OccPlus_spp_rich_bylatitude.jpeg"), 
       plot = lat_plot, width = 8, height = 6, units = "in", dpi = 300)

# Plot by altitude
alt_plot <- ggplot(occplus_all_data, aes(altitude, summed_probs)) +
  geom_point(aes(color = riverOrigin)) +
  geom_smooth(method = "gam", alpha = .2, aes(fill = riverOrigin))   +
  xlab(paste0("low", strrep(" ", 18), "Altitude (m)", strrep(" ", 18), "high")) +
  ylab("Summed species occupancies") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )
ggsave(here("outputs", "combine_results", "OccPlus_spp_rich_byaltitude.jpeg"), 
       plot = alt_plot, width = 8, height = 6, units = "in", dpi = 300)


