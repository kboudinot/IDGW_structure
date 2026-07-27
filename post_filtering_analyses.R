# IDGW structure 

# install packages 
library(pophelper)
library(tidyverse)
library(gridExtra)
library(grid)
library(sf)
library(scatterpie)
library(rnaturalearth)

# load structure outputs
sfiles <- list.files("/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/STRUCTURE/Structure_out_3", full.names = T)
slist <- readQ(files=sfiles, filetype = "structure")
summariseQ(tabulateQ(slist)) # k = 1-10, 5 runs each

# load metadata
metadata_structure <- read.csv("/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/sample_metadata_DAU.csv")
structure_data_df # from filtering script

# clean data for plotting
metadata_matched3 <- metadata_structure %>% 
  filter(Sample %in% structure_data_df$Indiv)

duplicated_samples3 <- metadata_matched3 %>%
  count(Sample) %>%
  filter(n > 1) %>%
  pull(Sample)

# pull all rows from metadata_matched that are duplicates
metadata_duplicates3 <- metadata_matched3 %>%
  filter(Sample %in% duplicated_samples3)

# remove duplicates in metadata, keeping the first occurrence
metadata_unique3<- metadata_matched3 %>%
  group_by(Sample) %>%
  slice(1) %>%    
  ungroup()

# match the order of structure_data$Indiv
metadata_ordered3 <- structure_data_df %>% 
  select(Indiv) %>%
  left_join(metadata_unique3, by = c("Indiv" = "Sample"))

# colors 
custom_structure_10 <- c(
  "#F0E442", "#009E73", "#E69F00", "#56B4E9", "#D55E00", "#CC79A7", "#0072B2","#FFFFFF", "#BBBBBB", "#000000")

                      
# plot k 2,3
k <- alignK(slist[c(11,16)])
k_update <- alignK(slist[c(11,16,31)])

Region <- metadata_ordered3$Region
grplab_df <- data.frame(
  Region = dplyr::recode(
    metadata_ordered3$Region,
    "magic_valley" = "Magic Valley",
    "southwest"    = "Southwest",
    "upper_snake"  = "Upper Snake",
    "clearwater"   = "Clearwater",
    "salmon"       = "Salmon",
    "panhandle"    = "Panhandle"
  )
)

p_sub <- plotQ(
  k,
  imgoutput = "join",
  clustercol = custom_structure_10,
  exportplot = FALSE,
  returnplot = TRUE,
  grplab = grplab_df,
  selgrp = "Region",
  subsetgrp = c(
    "Panhandle",
    "Clearwater",
    "Salmon",
    "Upper Snake",
    "Southwest",
    "Magic Valley"
  ),
  ordergrp = TRUE,
  basesize = 16,
  divsize = 0.25,
  grplabangle = 0.1,
  grplabsize = 4,
  grplabheight = 0.1,
  sharedindlab = FALSE,
  splab = paste0("K=", sapply(k, ncol))
)

grid.arrange(p_sub$plot[[1]])


library(gridExtra)


runs_to_plot <- alignK(slist[c(6,11,16,21,26,31,36,41,46)])
p_all <- plotQ(
  runs_to_plot,
  imgoutput = "join",
  clustercol = custom_structure_10,
  exportplot = FALSE,
  returnplot = TRUE,
  grplab = grplab_df,
  selgrp = "Region",
  subsetgrp = c(
    "Panhandle",
    "Clearwater",
    "Salmon",
    "Upper Snake",
    "Southwest",
    "Magic Valley"
  ),
  ordergrp = TRUE,
  basesize = 12,
  divsize = 0.25,
  grplabangle = 0.1,
  grplabsize = 4,
  grplabheight = 0.1,
  sharedindlab = FALSE,
  splab = paste0("K=", sapply(runs_to_plot, ncol))
)

# Plot all K values (K = 2–10)
grid.arrange(p_all$plot[[1]])





# PCA ----
library(adegenet)
library(tidyverse)

# data
structure_data_df3 
metadata_ordered3

# Drop the individual ID column from the genotype matrix
geno_df_3 <- structure_data_df3[, -which(colnames(structure_data_df3) == "Indiv")]

# Get base locus names (remove _allele1 / _allele2)
loci <- unique(sub("_allele[12]$", "", colnames(geno_df_3)))

# Combine allele1 and allele2 per locus
geno_df_3 <- sapply(loci, function(locus) {
  a1 <- geno_df_3[[paste0(locus, "_allele1")]]
  a2 <- geno_df_3[[paste0(locus, "_allele2")]]
  paste0(a1, a2)
})

# Convert to data.frame and keep rownames
geno_df_3 <- as.data.frame(geno_df_3)
rownames(geno_df_3) <- rownames(structure_data_df3)

# Build genind
genind_obj_3 <- df2genind(
  geno_df_3,
  sep = "",                 # split "12" -> allele 1 and allele 2
  ploidy = 2,
  type = "codom",
  ind.names = structure_data_df3$Indiv,
  NA.char = 0
)

# Population vector
pop_vec <- metadata_ordered3$Region
pop(genind_obj_3) <- factor(pop_vec)

##pca
geno_scaled <- scaleGen(genind_obj_3, center = T, scale = F, NA.method = "mean")
pca_res <- dudi.pca(geno_scaled, scannf = FALSE, nf = 3)
eig_values <- pca_res$eig
var_explained <- eig_values / sum(eig_values) * 100
var_explained[1:5]  # first 5 PCs

pca_df <- as.data.frame(pca_res$li)

pca_df$Indiv <- indNames(genind_obj_3)

pca_df$Population <- metadata_ordered3$Region

# rename regions
pca_df$Population <- recode(
  pca_df$Population,
  "panhandle"    = "Panhandle",
  "clearwater"   = "Clearwater",
  "salmon"       = "Salmon",
  "upper_snake"  = "Upper Snake",
  "southwest"    = "Southwest",
  "magic_valley" = "Magic Valley"
)

# plotting order
region_order <- c(
  "Panhandle",
  "Clearwater",
  "Salmon",
  "Upper Snake",
  "Southwest",
  "Magic Valley"
)

# sample sizes
pop_counts <- pca_df %>%
  count(Population, name = "N")

# Add counts to dataframe
pca_df <- pca_df %>%
  left_join(pop_counts, by = "Population") %>%
  mutate(
    Pop_Label = paste0(Population, " (n=", N, ")"),
    Pop_Label = factor(
      Pop_Label,
      levels = paste0(
        region_order,
        " (n=",
        pop_counts$N[match(region_order, pop_counts$Population)],
        ")"
      )
    )
  )

# colors
custom_palette <- c(
  "#DDCC77",# panhandle
  "#117733", 
  "#CC6677", # salmon
  "#332288",   # clearwater
  "#44AA99", # upper snake
  "#88CCEE" # southwest # magic valley 
)

cols <- setNames(custom_palette, levels(pca_df$Pop_Label))


# plot
ggplot(
  pca_df,
  aes(
    x = Axis1,
    y = Axis2,
    fill = Pop_Label
  )
) +
  geom_point(
    size = 5,
    shape = 21,
    color = "black",
    stroke = 0.6,
    alpha = 1
  ) +
  scale_fill_manual(values = cols) +
  coord_cartesian(ylim = c(-10, 10)) +
  coord_cartesian(
    ylim = range(pca_df$Axis2, na.rm = TRUE) + c(-0.5,0.5)
  ) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 2), "%)"),
    fill = "Region"
  ) +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 20),
    axis.ticks = element_line(linewidth = 0.6),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 22)
  )



# FST ----
library(hierfstat)
# genind
genind_obj_3@pop 

head(genind_obj_3@pop) # confirm same order as metadata
region_labels <- metadata_ordered3 %>%
  transmute(Region = str_to_title(Region)) %>%  
  pull(Region)

genind_obj_3@pop <- as.factor(region_labels)
levels(genind_obj_3@pop)
head(genind_obj_3@pop)

# convert genind to hierfstat format
hf3 <- genind2hierfstat(genind_obj_3)
hf3$pop <- droplevels(factor(hf3$pop))

table(genind_obj_3@pop)


# pairwise fst
pair_fst3 <- pairwise.WCfst(hf3)

new_order <- c("Panhandle", "Clearwater", "Salmon", "Upper_snake", "Southwest", "Magic_valley")

# Reorder rows and columns
pair_fst3_ordered <- pair_fst3[new_order, new_order]

#write.csv(pair_fst3_ordered, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/output/pairwise_fst_regions.csv")

# add boot.ppbetas
# a matrix with upper limit of the bootstrap CI above the diagonal and lower limit below the diagonal
fst_boot <- boot.ppbetas(
  dat     = hf3,
  nboot   = 999,
  quant   = c(0.025, 0.975),
  diploid = TRUE,
  digits  = 4
)

# check significance
x <- melt(fst_boot) %>% drop_na()

# significant if lower CI > 0
snp_signif <- x %>% filter(value > 0)



#### data not filtered for relatedness
geno_split_2



# He ----
genpop <- genind2genpop(genind_obj_3)
He <- adegenet::Hs(genpop) %>% as.data.frame()
colnames(He) <- c("He")
mean(He$He)
sd(He$He)
#write.csv(He, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/output/He_regions.csv")

genind_obj_3@pop
panhandle_genind <- genind_obj_3[pop="Panhandle"]
clearwater_genind <- genind_obj_3[pop="Clearwater"]
salmon_genind <- genind_obj_3[pop="Salmon"]
uppersnake_genind <- genind_obj_3[pop="Upper_snake"]
southwest_genind <- genind_obj_3[pop="Southwest"]
magicvalley_genind <- genind_obj_3[pop="Magic_valley"]

panhandle_Hstest  <- Hs.test(panhandle_genind, genind_obj_3, n.sim = 499, alter = "two-sided")
clearwater_Hstest  <- Hs.test(clearwater_genind, genind_obj_3, n.sim = 499, alter = "two-sided")
salmon_Hstest  <- Hs.test(salmon_genind, genind_obj_3, n.sim = 499, alter = "two-sided")
uppersnake_Hstest  <- Hs.test(uppersnake_genind, genind_obj_3, n.sim = 499, alter = "two-sided")
southwest_Hstest  <- Hs.test(southwest_genind, genind_obj_3, n.sim = 499, alter = "two-sided")
magicvalley_Hstest  <- Hs.test(magicvalley_genind, genind_obj_3, n.sim = 499, alter = "two-sided")

Hstests <- list(panhandle_Hstest,
                clearwater_Hstest,
                salmon_Hstest,
                uppersnake_Hstest,
                southwest_Hstest,
                magicvalley_Hstest)

hs_test_df <- bind_rows(
  lapply(Hstests, function(result) {
    data.frame(
      Region      = as.character(result$call$x),
      Observation = result$obs,
      Variance    = result$expvar[3],
      P_value     = result$pvalue
    )
  })
)
hs_test_df$SD <- sqrt(hs_test_df$Variance)
hs_test_df[, c("Region", "Observation", "SD", "P_value")]
# sig test adjusted
hs_test_df %>% filter(P_value <= 0.0083)  # bonferroni alpha should be 0.05/6 == 0.0083

# IBD ----
library(tidyverse)
library(vcfR)
library(adegenet)
library(vegan)
library(geosphere)

# load data
genind_obj_3
metadata_ordered3

# calculate genetic distance 
gen_dist <- dist(genind_obj_3, method = "euclidean")

# calculate geo distance
coords <- metadata_ordered3 %>%
  select(SCell_long, SCell_lat) %>%
  as.matrix()

# geodesic distance in meters
geo_dist <- distm(coords, fun = distGeo)

# convert to distance object
geo_dist <- as.dist(geo_dist)

# run mantel test
mantel_result <- mantel(
  gen_dist,
  geo_dist,
  method = "pearson",
  permutations = 9999
)

correlog <- vegan::mantel.correlog(gen_dist, D.geo = geo_dist, nperm=9999, cutoff = TRUE)
mantel_corr_df <- data.frame()

mantel_corr_df <- data.frame(
  class = correlog$mantel.res[, "class.index"],
  n_dist = correlog$mantel.res[, "n.dist"],
  mantel_r = correlog$mantel.res[, "Mantel.cor"],
  pvalue = correlog$mantel.res[, "Pr(Mantel)"]
)

mantel_corr_df


library(tidyverse)
df <- as.data.frame(correlog$mantel.res)

# Check column names
colnames(df)
# [1] "class.index" "n.dist" "Mantel.cor" "Pr(Mantel)" "Pr.corrected"

# Prepare dataframe
df <- df %>%
  filter(!is.na(Mantel.cor)) %>%                  # remove NAs
  mutate(
    `Distance Class (Kilometers)` = class.index / 1000, # convert m to km
  )

# Plot
p_mantel <- ggplot(df, aes(x = `Distance Class (Kilometers)`, y = Mantel.cor)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
  geom_path(color = "black") +
  theme_minimal() +
  labs(x = "Distance class (km)", y = "Mantel correlation") +
  scale_x_continuous(breaks = seq(0, max(df$`Distance Class (Kilometers)`), by = 50)) +
  scale_y_continuous(breaks = seq(floor(min(df$Mantel.cor)), ceiling(max(df$Mantel.cor)), by = 0.05)) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(color = "black", size = 30),  
    axis.text.y = element_text(color = "black", size = 30),  
    axis.title.x = element_text(size = 30),                  
    axis.title.y = element_text(size = 30),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  )

p_mantel



## study area map 

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)

idaho <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(name == "Idaho")

shapefile_path <- "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/region_shapefiles"
region_files <- list.files(shapefile_path, pattern = "\\.shp$", full.names = TRUE)

# Read and combine
regions_sf <- do.call(rbind, lapply(region_files, st_read))
regions_sf <- regions_sf %>%
  mutate(NAME = gsub(" Region", "", NAME))

# Convert to sf points
samples_sf <- st_as_sf(metadata, coords = c("SCell_long", "SCell_lat"), crs = 4326)

# colors
custom_palette <- c(
  "#DDCC77",# panhandle
  "#117733", # clearwater
  "#CC6677", # salmon
  "#332288",   # upper snake
  "#44AA99", # southwest 
  "#88CCEE",  # magic valley 
  "gray"
)



ggplot() +
  geom_sf(data = regions_sf, aes(fill = NAME), color = "black", alpha = 0.3) +
  geom_sf(data = samples_sf, aes(color = Region), size = 3, alpha = 0.8) +
  scale_fill_manual(values = custom_palette) +
  scale_color_manual(values = custom_palette) +
  theme_minimal() +
  labs(
    fill = "Region",
    color = "Region"
  )


# ---- Count samples per region ----
region_counts <- metadata_ordered3 %>%
  mutate(Region = recode(Region,
                         "panhandle"    = "Panhandle",
                         "clearwater"   = "Clearwater",
                         "salmon"       = "Salmon",
                         "upper_snake"  = "Upper Snake",
                         "southwest"    = "Southwest",
                         "magic_valley" = "Magic Valley")) %>%
  group_by(Region) %>%
  summarise(n_samples = n(), .groups = "drop")

regions_sf <- st_transform(regions_sf, crs = 4326)


region_centroids <- st_centroid(regions_sf) %>%
  select(Region, geometry) %>%
  left_join(region_counts, by = "Region")  # now n_samples exists

region_centroids_coords <- region_centroids %>%
  st_coordinates() %>%
  as.data.frame() %>%
  bind_cols(region_centroids %>% st_drop_geometry() %>% select(Region, n_samples)) %>%
  rename(lon = X, lat = Y)

# ---- Plot ----

#regions_sf <- regions_sf %>%
mutate(Region = recode(Region,
                       "Panhandle Region"   = "Panhandle",
                       "Clearwater Region"  = "Clearwater",
                       "Salmon Region"      = "Salmon",
                       "Upper Snake Region" = "Upper Snake",
                       "Southwest Region"   = "Southwest",
                       "Magic Valley Region"= "Magic Valley",
                       "Southeast Region"   = "Southeast"))  # optional

region_centroids_coords <- region_centroids_coords %>%
  mutate(Region = recode(Region,
                         "Panhandle Region"   = "Panhandle",
                         "Clearwater Region"  = "Clearwater",
                         "Salmon Region"      = "Salmon",
                         "Upper Snake Region" = "Upper Snake",
                         "Southwest Region"   = "Southwest",
                         "Magic Valley Region"= "Magic Valley",
                         "Southeast Region"   = "Southeast"))

# Define correct order
region_order <- c("Panhandle", "Clearwater", "Salmon", "Upper Snake", 
                  "Southwest", "Magic Valley", "Southeast")

# Re-order factor in both data frames
regions_sf <- regions_sf %>%
  mutate(Region = factor(NAME, levels = region_order))

region_centroids_coords <- region_centroids_coords %>%
  mutate(Region = factor(Region, levels = region_order))

# Now plot
idaho_map <- ggplot() +
  geom_sf(data = regions_sf, aes(fill = NAME), color = "black", alpha = 0.6) +
  geom_point(data = region_centroids_coords,
             aes(x = lon, y = lat, size = n_samples, fill = Region),
             shape = 21, color = "black", alpha = 1, na.rm = T) +
  scale_fill_manual(values = custom_palette) +
  scale_size_continuous(range = c(4, 20)) +
  theme_minimal() +
  labs(fill = "Region", size = "Number of samples")+
  theme(
    axis.title = element_blank()  # remove the "x" and "y" labels
  )

# --- USA map ---
usa <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(!name %in% c("Alaska", "Hawaii"))  # remove Alaska and Hawaii

# Highlight Idaho in the inset
idaho_sf <- usa %>% filter(name == "Idaho")

# Inset map of USA
usa_inset <- ggplot() +
  geom_sf(data = usa, fill = "gray90", color = "white") +
  geom_sf(data = idaho_sf, fill = "gray90", color = "black") +
  theme_void()

# Combine main map + inset
final_map <- ggdraw() +
  draw_plot(idaho_map) +
  draw_plot(usa_inset, x = 0.32, y = 0.67, width = 0.35, height = 0.35)

final_map








# --- K2 plot ---### Plot for K=2 across state
# identify K for each run
K_values <- sapply(slist, ncol)
# subset K = 2 runs
k2_runs <- slist[K_values == 2]
# align clusters across K = 2 runs
k2_aligned <- alignK(k2_runs)
# stack runs into array: individuals × clusters × runs
Q_array_k2 <- simplify2array(lapply(k2_aligned, as.matrix))

# average across runs (3rd dimension)
Q_mean_k2 <- apply(Q_array_k2, c(1, 2), mean)

# back to data.frame
Qmat_k2 <- as.data.frame(Q_mean_k2)
Qmat_k2$Indiv <- metadata_ordered3$Indiv
Q_indiv_k2 <- Qmat_k2 %>%
  left_join(
    metadata_ordered3 %>%
      select(Indiv, Region, SCell_long, SCell_lat),
    by = "Indiv"
  )
# average by Scell
Q_cell_k2 <- Q_indiv_k2 %>%
  group_by(SCell_long, SCell_lat) %>%
  summarise(
    across(starts_with("Cluster"), mean),
    n = n(),
    .groups = "drop"
  )

Q_cell_k2 <- Q_cell_k2 %>%
  mutate(
    pie_size = scales::rescale(n, to = c(0.1, 0.2))
  )


cluster_cols_k2 <- grep("^Cluster", colnames(Q_cell_k2), value = TRUE)
custom_structure_2 <- c(
  "#F0E442", "#009E73")
p_k2 <- ggplot() +
  geom_sf(
    data = regions_sf,
    fill = "grey95",
    color = "grey40",
    linewidth = 0.4,
    alpha = 0.6
  ) +
  geom_scatterpie(
    data = Q_cell_k2,
    aes(
      x = SCell_long,
      y = SCell_lat,
      r = pie_size
    ),
    cols = cluster_cols_k2,
    pie_scale = 1
  ) +
  scale_fill_manual(values = custom_structure_2) +
  coord_sf(
    xlim = c(-117, -110),
    ylim = c(42, 49)
  ) +
  theme_minimal() +
  labs(
    fill = "Cluster"
  ) +
  theme(
    axis.title = element_blank(),
    legend.position = "bottom"
  )



# --- K3 plot ---
# identify K for each run
K_values <- sapply(slist, ncol)
# subset K = 3 runs
k3_runs <- slist[K_values == 3]
# align clusters across K = 3 runs
k3_aligned <- alignK(k3_runs)
# stack runs into array: individuals × clusters × runs
Q_array_k3 <- simplify2array(lapply(k3_aligned, as.matrix))

# average across runs (3rd dimension)
Q_mean_k3 <- apply(Q_array_k3, c(1, 2), mean)

# back to data.frame
Qmat_k3 <- as.data.frame(Q_mean_k3)
Qmat_k3$Indiv <- metadata_ordered3$Indiv
Q_indiv_k3 <- Qmat_k3 %>%
  left_join(
    metadata_ordered3 %>%
      select(Indiv, Region, SCell_long, SCell_lat),
    by = "Indiv"
  )
# average by Scell
Q_cell_k3 <- Q_indiv_k3 %>%
  group_by(SCell_long, SCell_lat) %>%
  summarise(
    across(starts_with("Cluster"), mean),
    n = n(),
    .groups = "drop"
  )

Q_cell_k3 <- Q_cell_k3 %>%
  mutate(
    pie_size = scales::rescale(n, to = c(0.1, 0.2))
  )

cluster_cols_k3 <- grep("^Cluster", colnames(Q_cell_k3), value = TRUE)
custom_structure_3 <- c(
  "#F0E442", "#009E73","#E69F00")
p_k3 <- ggplot() +
  geom_sf(
    data = regions_sf,
    fill = "grey95",
    color = "grey40",
    linewidth = 0.4,
    alpha = 0.6
  ) +
  geom_scatterpie(
    data = Q_cell_k3,
    aes(
      x = SCell_long,
      y = SCell_lat,
      r = pie_size
    ),
    cols = cluster_cols_k3,
    pie_scale = 1
  )+
  scale_fill_manual(values = custom_structure_3) +
  coord_sf(
    xlim = c(-117, -110),
    ylim = c(42, 49)
  ) +
  theme_minimal() +
  labs(
    title = "K = 3",
    fill = "Cluster"
  ) +
  theme(
    axis.title = element_blank(),
    legend.position = "bottom"
  )

# --- Combine K2 + K3 maps with shared legend ---
library(patchwork)
library(cowplot)
main_map <- wrap_plots(
  p_k2,
  p_k3,
  guides = "collect"
) &
  theme(legend.position = "bottom")

# --- Final combined figure with USA inset and region key ---
final_map <- ggdraw() +
  draw_plot(main_map, x = 0, y = 0, width = 0.8, height = 1) +          # main K2/K3 maps
  draw_plot(usa_inset, x = 0.55, y = 0.67, width = 0.25, height = 0.25)  # USA inset

# --- Display the figure ---
final_map




# study area map plot

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)

# data 
metadata_ordered3

metadata <- metadata_ordered3

# Recode regions for nicer names
metadata <- metadata %>%
  mutate(
    Region = recode(
      Region,
      "magic_valley" = "Magic Valley",
      "southwest"    = "Southwest",
      "upper_snake"  = "Upper Snake",
      "clearwater"   = "Clearwater",
      "panhandle"    = "Panhandle",
      "salmon"       = "Salmon"
    )
  )

idaho <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(name == "Idaho")

shapefile_path <- "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/region_shapefiles"

# Read all shapefiles in the folder
region_files <- list.files(shapefile_path, pattern = "\\.shp$", full.names = TRUE)

# Read and combine
regions_sf <- do.call(rbind, lapply(region_files, st_read))
regions_sf <- regions_sf %>%
  mutate(NAME = gsub(" Region", "", NAME))

# Convert to sf points
samples_sf <- st_as_sf(metadata, coords = c("SCell_long", "SCell_lat"), crs = 4326)

# ---- Count samples per region ----
region_counts <- metadata_ordered3 %>%
  mutate(
    Region = recode(
      Region,
      "panhandle"    = "Panhandle",
      "clearwater"   = "Clearwater",
      "salmon"       = "Salmon",
      "upper_snake"  = "Upper Snake",
      "southwest"    = "Southwest",
      "magic_valley" = "Magic Valley"
    )
  ) %>%
  group_by(Region) %>%
  summarise(
    n_samples = n(),
    .groups = "drop"
  )

# Transform shapefile CRS
regions_sf <- regions_sf %>%
  rename(Region = NAME) %>%
  st_transform(crs = 4326)

# Add sample counts to regions
regions_sf <- regions_sf %>%
  left_join(region_counts, by = "Region")

# Create region label locations
region_centroids <- regions_sf %>%
  st_point_on_surface()

# Check result
region_centroids %>%
  select(Region, n_samples)

region_centroids_coords <- region_centroids %>%
  st_coordinates() %>%
  as.data.frame() %>%
  bind_cols(region_centroids %>% st_drop_geometry() %>% select(Region, n_samples)) %>%
  rename(lon = X, lat = Y)

# ---- Plot ----
custom_palette <- c(
  "#DDCC77",# panhandle
  "#117733", # clearwater
  "#CC6677", # salmon
  "#332288",   # upper snake
  "#44AA99", # southwest 
  "#88CCEE",  # magic valley 
  "gray"
)
regions_sf <- regions_sf %>%
mutate(Region = recode(Region,
                       "Panhandle Region"   = "Panhandle",
                       "Clearwater Region"  = "Clearwater",
                       "Salmon Region"      = "Salmon",
                       "Upper Snake Region" = "Upper Snake",
                       "Southwest Region"   = "Southwest",
                       "Magic Valley Region"= "Magic Valley",
                       "Southeast Region"   = "Southeast"))  # optional

region_centroids_coords <- region_centroids_coords %>%
  mutate(Region = recode(Region,
                         "Panhandle Region"   = "Panhandle",
                         "Clearwater Region"  = "Clearwater",
                         "Salmon Region"      = "Salmon",
                         "Upper Snake Region" = "Upper Snake",
                         "Southwest Region"   = "Southwest",
                         "Magic Valley Region"= "Magic Valley",
                         "Southeast Region"   = "Southeast"))

# Define correct order
region_order <- c("Panhandle", "Clearwater", "Salmon", "Upper Snake", 
                  "Southwest", "Magic Valley", "Southeast")

# Re-order factor in both data frames
regions_sf <- regions_sf %>%
  mutate(
    Region = factor(Region, levels = region_order)
  )

region_centroids_coords <- region_centroids_coords %>%
  mutate(Region = factor(Region, levels = region_order))

# Custom legend labels with sample sizes
region_labels <- c(
  "Panhandle"    = "Panhandle (n = 512)",
  "Clearwater"   = "Clearwater (n = 497)",
  "Salmon"       = "Salmon (n = 224)",
  "Upper Snake"  = "Upper Snake (n = 102)",
  "Southwest"    = "Southwest (n = 221)",
  "Magic Valley" = "Magic Valley (n = 83)",
  "Southeast"    = "Southeast (n = 1)"
)

# Now plot
idaho_map <- ggplot() +
  geom_sf(
    data = regions_sf,
    aes(fill = Region),
    color = "black",
    alpha = 0.3
  ) +
  scale_fill_manual(
    values = custom_palette,
    labels = region_labels
  ) +
  theme_minimal() +
  labs(
    fill = "Region"
  ) +
 theme(
    axis.title = element_blank(),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18)
  )

# --- USA map ---
usa <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(!name %in% c("Alaska", "Hawaii"))

# Highlight Idaho in the inset
idaho_sf <- usa %>% filter(name == "Idaho")

# Inset map
usa_inset <- ggplot() +
  geom_sf(data = usa, fill = "gray90", color = "white") +
  geom_sf(data = idaho_sf, fill = "gray90", color = "black") +
  theme_void()

# Combine main map + inset
final_map <- ggdraw() +
  draw_plot(idaho_map) +
  draw_plot(usa_inset, x = 0.32, y = 0.67, width = 0.35, height = 0.35)

final_map






# unfiltered data FST, He
library(adegenet)
library(hierfstat)
library(tidyverse)

########################################################
# Load metadata and filter samples
########################################################

metadata_structure <- read.csv(
  "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/sample_metadata_DAU.csv"
)

# Keep only samples with metadata
filt2_df <- filt2_df %>%
  filter(Indiv_short %in% metadata_structure$Sample)

# Keep only samples with complete metadata
metadata_complete <- metadata_structure %>%
  filter(Sample %in% filt2_df$Indiv_short) %>%
  filter(if_all(everything(), ~ !is.na(.)))

filt2_df <- filt2_df %>%
  filter(Indiv_short %in% metadata_complete$Sample)


# Split each genotype into allele1 and allele2
geno_split_2 <- filt2_df %>%
  select(Indiv, Indiv_long, Indiv_short) %>%
  bind_cols(
    lapply(
      filt2_df %>% select(-Indiv, -Indiv_long, -Indiv_short),
      function(x) {
        alleles <- stringr::str_split_fixed(x, "", 2)
        
        tibble(
          allele1 = alleles[, 1],
          allele2 = alleles[, 2]
        )
      }
    ) |>
      setNames(names(filt2_df)[-(1:3)]) |>
      purrr::imap_dfc(~tibble(
        !!paste0(.y, "_allele1") := .x$allele1,
        !!paste0(.y, "_allele2") := .x$allele2
      ))
  )

# Make allele1 >= allele2
loci <- unique(sub("_allele[12]$", "", names(geno_split_2)[-(1:3)]))

for (locus in loci) {
  
  a1 <- paste0(locus, "_allele1")
  a2 <- paste0(locus, "_allele2")
  
  allele1 <- geno_split_2[[a1]]
  allele2 <- geno_split_2[[a2]]
  
  geno_split_2[[a1]] <- pmax(allele1, allele2)
  geno_split_2[[a2]] <- pmin(allele1, allele2)
}

# Remove ID columns
geno_df_unfilt <- geno_split_2 %>%
  select(-Indiv, -Indiv_long, -Indiv_short)

# Combine allele1 and allele2 back into one genotype string
loci <- unique(sub("_allele[12]$", "", colnames(geno_df_unfilt)))

geno_df_unfilt <- sapply(loci, function(locus) {
  
  a1 <- geno_df_unfilt[[paste0(locus, "_allele1")]]
  a2 <- geno_df_unfilt[[paste0(locus, "_allele2")]]
  
  paste0(a1, a2)
})

# Replace NANA with NA
geno_df_unfilt[geno_df_unfilt == "NANA"] <- NA

# Convert back to data.frame
geno_df_unfilt <- as.data.frame(geno_df_unfilt)

# Remove prop_missing if present
geno_df_unfilt <- geno_df_unfilt %>%
  select(-any_of("prop_missing"))

# Check columns
dim(geno_df_unfilt)
colnames(geno_df_unfilt)

# Remove prop_missing column
geno_df_unfilt <- geno_df_unfilt %>%
  select(-any_of("prop_missing"))

# Check
table(is.na(geno_df_unfilt))


# Confirm genotype lengths (should only have length 2 or NA)
table(nchar(unlist(geno_df_unfilt)), useNA = "ifany")



# Create genind
unfilt_genind <- df2genind(
  geno_df_unfilt,
  sep = "",
  ploidy = 2,
  type = "codom",
  ind.names = filt2_df$Indiv_short,
  NA.char = "NA"
)

# Check
nInd(unfilt_genind)
nLoc(unfilt_genind)


# Combine with metadata 
# metadata
metadata_structure <- read.csv("/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/sample_metadata_DAU.csv")
metadata_structure$Sample
indNames(unfilt_genind) <- make.unique(filt2_df$Indiv_short)

region_lookup <- setNames(
  metadata_structure$Region,
  metadata_structure$Sample
)

# make region
matched_region <- region_lookup[indNames(unfilt_genind)]
# check
sum(is.na(matched_region)) # should be 0
length(matched_region) # 1640
# should match 
nInd(unfilt_genind) #1640

# assign to genind
pop(unfilt_genind) <- as.factor(matched_region)
unfilt_pop <-as.factor(matched_region)
# check
table(pop(unfilt_genind))

# remove samples with NA region 
na_inds <- is.na(pop(unfilt_genind)) 
unfilt_genind <- unfilt_genind[!na_inds, ]
pop_vector <- pop(unfilt_genind)

# convert to integer factor (hierfstat prefers numeric/pop as first column)
pop_numeric <- as.numeric(as.factor(pop_vector))

unfilt_genind
unfilt_hf <- genind2hierfstat(unfilt_genind)


########################################################
# FST
########################################################

unfilt_hf <- genind2hierfstat(unfilt_genind)

# Pairwise Weir & Cockerham FST
unfilt_fst <- pairwise.WCfst(unfilt_hf)

# Region order
new_order <- c(
  "panhandle",
  "clearwater",
  "salmon",
  "upper_snake",
  "southwest",
  "magic_valley",
  "southeast"
)

new_order <- new_order[new_order %in% rownames(unfilt_fst)]

unfilt_fst_ordered <- unfilt_fst[new_order, new_order]

write.csv(
  unfilt_fst_ordered,
  "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/output/fst_full.csv",
  row.names = TRUE
)


########################################################
# FST bootstrap
########################################################

  unfilt_fst_boot <- boot.ppbetas(
  dat = unfilt_hf,
  nboot = 999,
  quant = c(0.025, 0.975),
  diploid = TRUE,
  digits = 4
)

fst_significant <- reshape2::melt(unfilt_fst_boot) %>%
  drop_na() %>%
  filter(value > 0)


########################################################
# Expected heterozygosity (He)
########################################################

unfilt_genpop <- genind2genpop(unfilt_genind)

popNames(unfilt_genpop)
unfilt_genpop <- unfilt_genpop[popNames(unfilt_genpop) != "southeast"]
popNames(unfilt_genpop)

unfilt_He <- adegenet::Hs(unfilt_genpop) %>%
  as.data.frame()

colnames(unfilt_He) <- "He"

# remove regions with no samples
unfilt_He <- unfilt_He %>%
  filter(!is.na(He))

mean(unfilt_He$He)

unfilt_He

unfilt_genind@pop
panhandle_genind <- unfilt_genind[pop="panhandle"]
clearwater_genind <- unfilt_genind[pop="clearwater"]
salmon_genind <- unfilt_genind[pop="salmon"]
uppersnake_genind <- unfilt_genind[pop="upper_snake"]
southwest_genind <- unfilt_genind[pop="southwest"]
magicvalley_genind <- unfilt_genind[pop="magic_valley"]

panhandle_Hstest  <- Hs.test(panhandle_genind, unfilt_genind, n.sim = 499, alter = "two-sided")
clearwater_Hstest  <- Hs.test(clearwater_genind, unfilt_genind, n.sim = 499, alter = "two-sided")
salmon_Hstest  <- Hs.test(salmon_genind, unfilt_genind, n.sim = 499, alter = "two-sided")
uppersnake_Hstest  <- Hs.test(uppersnake_genind, unfilt_genind, n.sim = 499, alter = "two-sided")
southwest_Hstest  <- Hs.test(southwest_genind, unfilt_genind, n.sim = 499, alter = "two-sided")
magicvalley_Hstest  <- Hs.test(magicvalley_genind, unfilt_genind, n.sim = 499, alter = "two-sided")

Hstests <- list(panhandle_Hstest,
                clearwater_Hstest,
                salmon_Hstest,
                uppersnake_Hstest,
                southwest_Hstest,
                magicvalley_Hstest)

hs_test_df <- bind_rows(
  lapply(Hstests, function(result) {
    data.frame(
      Region      = as.character(result$call$x),
      Observation = result$obs,
      Variance    = result$expvar[3],
      P_value     = result$pvalue
    )
  })
)
hs_test_df$SD <- sqrt(hs_test_df$Variance)
hs_test_df[, c("Region", "Observation", "SD", "P_value")]
# sig test adjusted
hs_test_df %>% filter(P_value <= 0.0083)  # bonferroni alpha should be 0.05/6 == 0.0083

