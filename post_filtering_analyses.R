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
 "#ffffea", "#6c584c", "#adc178", "#a98467", "#e7e8c3", "#9a7961", "#dde5b4", "#c5d396", "#aba370", "#58554e")

custom_structure <-  c(
  "#779c46",  # Panhandle
  "#b2a24a",# Salmon
  "#c5d396", #Clearwater
  "#9b8965", # Upper snake
  "#6c584c",  # southwest
  "#58554e"  # magic valley 
)
                      
# plot k 2,3,4
k <- alignK(slist[c(11,16)])
k_update <- alignK(slist[c(11,16,31)])

labels_region$Region <- dplyr::recode(
  labels_region$Region,
  "magic_valley" = "Magic Valley",
  "southwest"    = "Southwest",
  "upper_snake"  = "Upper Snake",
  "clearwater"   = "Clearwater",
  "salmon"       = "Salmon",
  "panhandle"    = "Panhandle"
)

p_sub <- plotQ(
  k,
  imgoutput = "join",
  clustercol = custom_structure_10,
  exportplot = FALSE,
  returnplot = TRUE,
  grplab = labels_region,   # must contain Region column
  selgrp = "Region",        # required, cannot be NULL
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
  grplabsize =4,
  grplabheight = 0.1,         # hides the group label panel
  sharedindlab = FALSE,
  splab = paste0("K=", sapply(k, ncol))
)


grid.arrange(p_sub$plot[[1]])


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
  "#ffffea",   # Panhandle
  "#adc178",   # Clearwater
  "#6c584c",   # Salmon
  "#a98467",   # Upper Snake
  "#e7e8c3",   # Southwest
  "#58554e"    # Magic Valley
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

panhandle_Hstest  <- Hs.test(panhandle_genind, genind_obj_3_filt, n.sim = 499, alter = "two-sided")
clearwater_Hstest  <- Hs.test(clearwater_genind, genind_obj_3_filt, n.sim = 499, alter = "two-sided")
salmon_Hstest  <- Hs.test(salmon_genind, genind_obj_3_filt, n.sim = 499, alter = "two-sided")
uppersnake_Hstest  <- Hs.test(uppersnake_genind, genind_obj_3_filt, n.sim = 499, alter = "two-sided")
southwest_Hstest  <- Hs.test(southwest_genind, genind_obj_3_filt, n.sim = 499, alter = "two-sided")
magicvalley_Hstest  <- Hs.test(magicvalley_genind, genind_obj_3_filt, n.sim = 499, alter = "two-sided")

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
