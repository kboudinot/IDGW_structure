# Load packages
library(adegenet)
library(memgene)

library(geosphere)
library(sf)

library(tidyverse)
library(dplyr)

library(doParallel)
library(foreach)

sample_sex <- read.csv("sample_metadata_regions.csv") %>% select(Sample, Sex)
sample_metadata_DAU <- read.csv("sample_metadata_DAU.csv") 

metadata <- left_join(sample_metadata_DAU, sample_sex)

metadata$Sex[metadata$Sex == "M"] <- "Male"
metadata$Sex[metadata$Sex == "F"] <- "Female"
metadata$Sex[metadata$Sex == "U"] <- "Unknown"
metadata$Sex[is.na(metadata$Sex)] <- "Unknown"

metadata <- metadata %>%
  mutate(Sample = str_remove(Sample, "^CluIDFG")) %>%        # remove prefix
  mutate(Sample = str_replace_all(Sample, "_", ".")) %>%     # replace underscores
  mutate(Sample = str_remove(Sample, "\\.(initial|f1).*$")) %>%  # remove .initial or .f1 + anything after
  filter(!is.na(Region)) %>%
  distinct()

structure_data_rel_filt <- read.table("structure_data_3.txt", header = T) # from filtering script

length(structure_data_rel_filt$Indiv)

genos_meta_rel_filt <- right_join(metadata, structure_data_rel_filt, by = c("Sample" = "Indiv"))

genind_rel_filt@pop <- as.factor(genos_meta_rel_filt$Sex)

genind_rel_filt_sep <- seppop(genind_rel_filt)

genind_rel_filt_F <- genind_rel_filt_sep$Female
genind_rel_filt_M <- genind_rel_filt_sep$Male

# Filtered for relatedness
genD_rf_all <- dist(genind_rel_filt, method = "euclidean")
## Male only
genD_rf_f <- dist(genind_rel_filt_F, method = "euclidean")
## Female only
genD_rf_m <- dist(genind_rel_filt_M, method = "euclidean")

# Filtered for related
pts <- st_as_sf(genos_meta_rel_filt, coords = c("SCell_long", "SCell_lat"), crs = 4326)
coords_rf_all <- st_coordinates(pts)

meta_rel_filt_F <- genos_meta_rel_filt %>% filter(Sex == "Female")

pts <- st_as_sf(meta_rel_filt_F, coords = c("SCell_long", "SCell_lat"), crs = 4326)
coords_rf_f <- st_coordinates(pts)

meta_rel_filt_M <- genos_meta_rel_filt %>% filter(Sex == "Male")

pts <- st_as_sf(meta_rel_filt_M, coords = c("SCell_long", "SCell_lat"), crs = 4326)
coords_rf_m <- st_coordinates(pts)

# Set up datasets as a named list
# genD matrix and coords
dataset_list <- list(
  rel_filt_all = list(genD = genD_rf_all, coords = coords_rf_all),
  rel_filt_m = list(genD = genD_rf_m, coords = coords_rf_m),
  rel_filt_f = list(genD = genD_rf_f, coords = coords_rf_f)
)

n_cores <- 3
cl <- makeCluster(n_cores)
clusterCall(cl, function(libs) .libPaths(libs), .libPaths())  # <-- add this
registerDoParallel(cl)

start_time <- Sys.time()

results <- foreach(
  name = names(dataset_list),
  .packages = "memgene",
  .errorhandling = "pass"
) %dopar% {
  d <- dataset_list[[name]]
  mgQuick(
    d$genD, d$coords,
    forwardPerm = 999,
    finalPerm = 999,
    longlat = TRUE,
    verbose = FALSE
  )
}

end_time <- Sys.time()
elapsed_time <- end_time - start_time

names(results) <- names(dataset_list)
stopCluster(cl)
sapply(results, function(x) inherits(x, "error"))

saveRDS(results, "related_kept_memgene.RData")

# Results

results <- readRDS("related_filt_memgene.RData")

results$rel_filt_all$P
results$rel_filt_m$P
results$rel_filt_f$P

results$rel_filt_all$RsqAdj
results$rel_filt_m$RsqAdj
results$rel_filt_f$RsqAdj

# Visualize results

library(ggplot2)
library(patchwork)
library(sf)
library(ggnewscale)
library(ggspatial)

region_palette <- c(
  "#DDCC77", "#117733", "#CC6677", "#332288", "#44AA99", "#88CCEE", "gray"
)

# --- Reusable function for one spatial panel ---
make_spatial_plot <- function(name, panel_label, mg_axis = 1, size_limits) {
  result <- results[[name]]
  coords <- dataset_list[[name]]$coords
  
  plot_data <- data.frame(
    X = coords[, 1],
    Y = coords[, 2],
    mg = result$memgene[, mg_axis]
  )
  plot_data$mg_sign <- ifelse(plot_data$mg >= 0, "positive", "negative")
  plot_data$mg_size <- abs(plot_data$mg)
  
  points_sf <- st_as_sf(plot_data, coords = c("X", "Y"), crs = 4326)
  points_sf <- st_transform(points_sf, st_crs(regions_sf))
  
  ggplot() +
    geom_sf(data = regions_sf, aes(fill = NAME), color = "black", alpha = 0.3) +
    scale_fill_manual(values = region_palette, guide = "none") +
    new_scale_fill() +
    geom_sf(data = points_sf, aes(size = mg_size, fill = mg_sign),
            shape = 21, color = "black", stroke = 0.5) +   # alpha removed -> fully opaque
    scale_fill_manual(values = c("positive" = "black", "negative" = "white")) +
    scale_size_continuous(range = c(0.5, 3), limits = size_limits) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    labs(title = panel_label, size = "Magnitude", fill = "Sign") +
    theme_minimal(base_size = 9) +
    theme(axis.text = element_text(size = 10))+
    annotation_scale(location = "tr", width_hint = 0.3,
                     unit_category = "metric", style = "bar",
                     text_cex = 0.6, height = unit(0.15, "cm"))
}

# --- Names of the three datasets to compare, in this order ---
group_names  <- c("rel_filt_all", "rel_filt_m", "rel_filt_f")
panel_labels <- c("A) All", "B) Male", "C) Female")

# --- mg1 figure ---
mg1_vals   <- unlist(lapply(group_names, function(n) abs(results[[n]]$memgene[, 1])))
mg1_limits <- c(0, max(mg1_vals))

mg1_plots <- Map(function(n, l) make_spatial_plot(n, l, mg_axis = 1, size_limits = mg1_limits),
                 group_names, panel_labels)

mg1_tiled <- wrap_plots(mg1_plots, ncol = 3, guides = "collect") +
  plot_annotation(title = "Spatial MEMGENE Axis 1 (mg1)")

# --- mg2 figure ---
mg2_vals   <- unlist(lapply(group_names, function(n) abs(results[[n]]$memgene[, 2])))
mg2_limits <- c(0, max(mg2_vals))

mg2_plots <- Map(function(n, l) make_spatial_plot(n, l, mg_axis = 2, size_limits = mg2_limits),
                 group_names, panel_labels)

mg2_tiled <- wrap_plots(mg2_plots, ncol = 3, guides = "collect") +
  plot_annotation(title = "Spatial MEMGENE Axis 2 (mg2)")

mg1_tiled
mg2_tiled

ggsave("mg1_tiled.png", mg1_tiled, width = 12, height = 4.5, dpi = 600, bg = "white")
ggsave("mg2_tiled.png", mg2_tiled, width = 12, height = 4.5, dpi = 600, bg = "white")
