# Load packages
library(tidyverse)
library(adegenet) 
library(poppr)
library(ade4)
library(RColorBrewer)
library(plotly)
library(related) 
library(igraph)

# Load raw genotype data
snp_raw <- read.csv("/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/microhaplotype_genotypes_raw.csv", na.strings = "")

# Clean raw data
snp_df <- snp_raw %>%
  mutate(Indiv = Indiv |>              
      basename() |>                      # remove full path
      sub("\\.1\\.bam$", "", x = _)      # remove trailing ".1.bam"
  )

# Code microhaplotypes
# Identify all unique alleles per locus
alleles_by_locus <- snp_df %>%
  select(Locus, Allele1, Allele2) %>%
  pivot_longer(cols = c(Allele1, Allele2), values_to = "Allele") %>%
  filter(!is.na(Allele)) %>%      
  distinct(Locus, Allele) %>%      
  group_by(Locus) %>%
  mutate(Allele_Index = dense_rank(Allele)) %>%
  ungroup()

# Attach numeric indices and make ordered genotype codes
df_encoded <- snp_df %>%
  left_join(alleles_by_locus,
            by = c("Locus", "Allele1" = "Allele")) %>%
  rename(A1_Code = Allele_Index) %>%
  left_join(alleles_by_locus,
            by = c("Locus", "Allele2" = "Allele")) %>%
  rename(A2_Code = Allele_Index) %>%
  mutate(
    Genotype_Code = ifelse(
      is.na(A1_Code) | is.na(A2_Code),
      NA_character_,
      sprintf("%d%d", A1_Code, A2_Code)  # keep A1/A2 order exactly
    )
  )

# Long to wide: one row per individual
genotype_matrix <- df_encoded %>%
  select(Indiv, Locus, Genotype_Code) %>%
  pivot_wider(
    names_from = Locus,
    values_from = Genotype_Code
  ) %>%
  mutate(
    Indiv_long  = Indiv,
    Indiv_short = str_extract(Indiv, "Clu.*|NTC.*"),
    Indiv_short = if_else(is.na(Indiv_short), Indiv, Indiv_short)
  ) %>%
  select(Indiv, Indiv_long, Indiv_short, everything())  

# filter 
genotype_missing <- genotype_matrix %>% 
  filter(if_all(4:343, is.na))   #rows where all genotypes are NA

genotype_missing <- genotype_missing$Indiv 

# Make the dataframe a genind object
snp_genind <- df2genind(
  genotype_matrix %>% select(-Indiv, -Indiv_short, -Indiv_long), 
  ind.names = genotype_matrix$Indiv,   
  ploidy = 2,
  sep = "",
  NA.char = "NA"
)

################################################################################
# Filter for missingness 
genind_loci_filt <- missingno(snp_genind,
                              type = "loci", 
                              cutoff = 0.20, 
                              quiet = FALSE, 
                              freq = FALSE)
genind_ind_filt <- missingno(genind_loci_filt,
                             type = "genotype", 
                             cutoff = 0.20, 
                             quiet = FALSE, 
                             freq = FALSE)

# Vector of samples and loci removed
removed_loci_1 <- setdiff(locNames(snp_genind), locNames(genind_loci_filt))
removed_sample_1 <- setdiff(indNames(snp_genind), indNames(genind_ind_filt)) #330
# add individuals that were removed by df2genind that have no scored loci 
removed_sample_1 <- c(removed_sample_1, genotype_missing)

# Combined genind for initial loci and genotype missingness filter
genind_filt1 <- snp_genind %>%
  missingno(type = "loci", cutoff = 0.20, quiet = FALSE, freq = FALSE) %>%
  missingno(type = "genotype", cutoff = 0.20, quiet = FALSE, freq = FALSE)

          # Genotype missigness PCA
          ind_scale <- scaleGen(genind_ind_filt, center = TRUE, scale = FALSE, NA.method = "mean")
          pca_res <- dudi.pca(ind_scale, scannf = FALSE, nf = 3)  
          
          # missing per individual
          geno <- tab(genind_ind_filt)            # extract SNP matrix
          indMissing <- rowMeans(is.na(geno))   # proportion missing per sample
          
          # Extract PCA scores (li) and percent variance
          pca_scores <- as.data.frame(pca_res$li)
          var_expl <- round(pca_res$eig / sum(pca_res$eig) * 100, 2)
          
          # Add metadata
          pca_scores$pop <- pop(snp_genind)  # or pop(genind) if unfiltered
          pca_scores$missing <- indMissing
          
          
##################################################################################          
# Filter diagnostic and missingness
    # remove samples and loci filtered for initial missingness
    snp_df_filt <- snp_df %>%
      filter(!Indiv %in% removed_sample_1) %>% 
      filter(!Locus %in% removed_loci_1)
    
    # Clean
    snps <- snp_df_filt %>%
      mutate(
        genotype = ifelse(
          is.na(Allele1) & is.na(Allele2),
          NA,
          paste0(Allele1, "/", Allele2)
        )
      ) %>%
      select(Indiv, Locus, genotype) %>% # long name
      pivot_wider(
        names_from = Locus,
        values_from = genotype
      )

# read diagnostic loci key
diagnostic_loci <- read.csv("/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/neutral_diagnostic_loci_new(in).csv", 
                            header = T) %>% 
  select(1:2) %>% 
  mutate(OG_locus_ID = sub("_[^_]*$", "", OG_locus_ID))
diagnostic_dog <- diagnostic_loci %>% 
  filter(Diagnostic == "Dog")
diagnostic_coyote <- diagnostic_loci %>% 
  filter(Diagnostic == "Coyote")


# Subset genotype_matrix: keep Indiv and only matching loci
loci_to_keep <- diagnostic_loci$OG_locus_ID[
  diagnostic_loci$OG_locus_ID != "Clu_chr9_24084144"]
loci_to_keep_dog <- diagnostic_dog$OG_locus_ID
loci_to_keep_dog <- loci_to_keep_dog[loci_to_keep_dog != "Clu_chr9_24084144"]
loci_to_keep_coyote <- diagnostic_coyote$OG_locus_ID
genotype_diagnostic <- snps %>%
  select(Indiv, all_of(loci_to_keep))
genotype_diagnostic_dog <- snps %>%
  select(Indiv, all_of(loci_to_keep_dog))
genotype_diagnostic_coyote <- snps %>%
  select(Indiv, all_of(loci_to_keep_coyote))

# df to genind
genotype_diagnostic_genind <- df2genind(
  genotype_diagnostic %>% select(-Indiv),  # remove column for loci
  ind.names = genotype_diagnostic$Indiv,   # assign Indiv as individual names
  ploidy = 2,
  sep = "/",
  NA.char = "NA")
genotype_diagnostic_dog_genind <- df2genind(
  genotype_diagnostic_dog %>% select(-Indiv),  # remove column for loci
  ind.names = genotype_diagnostic_dog$Indiv,   # assign Indiv as individual names
  ploidy = 2,
  sep = "/",
  NA.char = "NA")
genotype_diagnostic_coyote_genind <- df2genind(
  genotype_diagnostic_coyote %>% select(-Indiv),  # remove column for loci
  ind.names = genotype_diagnostic_coyote$Indiv,   # assign Indiv as individual names
  ploidy = 2,
  sep = "/",
  NA.char = "NA")

# filter for missingness 
genotype_diagnostic_genind_filt_loci <- missingno(genotype_diagnostic_genind, type = "loci", cutoff = 0.20, quiet = FALSE, freq = FALSE)
genotype_diagnostic_genind_filt_genotype <- missingno(genotype_diagnostic_genind_filt_loci, type = "genotype", cutoff = 0.20, quiet = FALSE, freq = FALSE)
genotype_diagnostic_genind <- genotype_diagnostic_genind_filt_genotype
genotype_diagnostic_dog_genind_filt_loci <- missingno(genotype_diagnostic_dog_genind, type = "loci", cutoff = 0.20, quiet = FALSE, freq = FALSE)
genotype_diagnostic_dog_genind_filt_genotype <- missingno(genotype_diagnostic_dog_genind_filt_loci, type = "genotype", cutoff = 0.20, quiet = FALSE, freq = FALSE)
genotype_diagnostic_dog_genind <- genotype_diagnostic_dog_genind_filt_genotype
genotype_diagnostic_coyote_genind_filt_loci <- missingno(genotype_diagnostic_coyote_genind, type = "loci", cutoff = 0.20, quiet = FALSE, freq = FALSE)
genotype_diagnostic_coyote_genind_filt_genotype <- missingno(genotype_diagnostic_coyote_genind_filt_loci, type = "genotype", cutoff = 0.20, quiet = FALSE, freq = FALSE)
genotype_diagnostic_coyote_genind <- genotype_diagnostic_coyote_genind_filt_genotype

# samples removed when filtering based on diagnostic
removed_loci_diagnostic <- locNames(genotype_diagnostic_genind)
#writeLines(removed_loci_diagnostic, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/removed_loci_diagnostic.txt")
removed_sample_2a <- setdiff(indNames(genind_filt1), indNames(genotype_diagnostic_genind_filt_genotype)) # 24 samples removed based on all diagnostic 
#writeLines(removed_sample_2a, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/removed_sample_2a.txt")
removed_sample_2b <- setdiff(indNames(genind_filt1), indNames(genotype_diagnostic_dog_genind_filt_genotype)) # 83 samples removed based on dog diagnostic
#writeLines(removed_sample_2b, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/removed_sample_2b.txt")
removed_sample_2c <- setdiff(indNames(genind_filt1), indNames(genotype_diagnostic_coyote_genind_filt_genotype)) # 6 samples removed based on coyote diagnostic
#writeLines(removed_sample_2c, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/removed_sample_2c.txt")


# Making diagnostic pca
# import metadata
metadata <- read.csv("/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/sample_metadata_regions.csv")
# match to genind
head(indNames(genotype_diagnostic_genind))
head(locNames(genotype_diagnostic_genind))
head(metadata$Sample) # need to add long sample names back 
          
geno_df <- as.data.frame(genotype_diagnostic_genind$tab)
geno_df <- geno_df %>%
  mutate(
  Sample_long = rownames(.),                            # keep the full name
  Sample      = str_extract(Sample_long, "Clu.*|NTC.*") # extract short ID
  )

metadata_df <- geno_df %>%
  left_join(metadata, by = "Sample", relationship = "many-to-many")
          
# create key for harvest and dog 
metadata_key <- metadata_df %>%
  select(Sample, Sample_long, MortYear, GMU, Region) %>%
  mutate(
  key  = if_else(!is.na(MortYear) & !is.na(GMU) & !is.na(Region), "harvest", "not harvest"),
  type = if_else(str_detect(Sample, "CluVARI"), "dog", "NA")  # use "NA" string so unite works
  ) %>%
  unite("key_type", key, type, sep = "_")  # combine into one column
          
key <- metadata_key$key_type[match(indNames(genotype_diagnostic_genind), metadata_df$Sample_long)]
key_dog <- metadata_key$key_type[match(indNames(genotype_diagnostic_dog_genind), metadata_df$Sample_long)]
key_coyote <- metadata_key$key_type[match(indNames(genotype_diagnostic_coyote_genind), metadata_df$Sample_long)]
          
pop(genotype_diagnostic_genind) <- factor(key) 
head(pop(genotype_diagnostic_genind))
pop(genotype_diagnostic_dog_genind) <- factor(key_dog) 
pop(genotype_diagnostic_coyote_genind) <- factor(key_coyote) 
          
# Convert to numeric allele frequency matrix
genotype_diagnostic_genind_scale <- scaleGen(genotype_diagnostic_genind, center = TRUE, scale = FALSE, NA.method = "mean")
genotype_diagnostic_dog_genind_scale <- scaleGen(genotype_diagnostic_dog_genind, center = TRUE, scale = FALSE, NA.method = "mean")
genotype_diagnostic_coyote_genind_scale <- scaleGen(genotype_diagnostic_coyote_genind, center = TRUE, scale = FALSE, NA.method = "mean")
          


# Run PCA using ade4
pca_1 <- dudi.pca(genotype_diagnostic_genind_scale, scannf = FALSE, nf = 3)
pca_2 <- dudi.pca(genotype_diagnostic_dog_genind_scale, scannf = FALSE, nf = 3)
pca_3 <- dudi.pca(genotype_diagnostic_coyote_genind_scale, scannf = FALSE, nf = 3) 
cols <- brewer.pal(n = nPop(genotype_diagnostic_genind), name = "Set1")
cols_2 <- brewer.pal(n = nPop(genotype_diagnostic_dog_genind), name = "Set1")
cols_3 <- brewer.pal(n = nPop(genotype_diagnostic_coyote_genind), name = "Set1")
          
# PCA scatter plot
s.class(pca_1$li,
  pop(genotype_diagnostic_genind),
  col = cols,
  axesell = T,
  cstar = 0,
  cellipse = 0,
  xlim = c(-40, 10),
  ylim = c(-20, 20)
  )

# Add a legend/key on the side
legend(x = -15,  # adjust X coordinate to move horizontally
  y = 15,   # adjust Y coordinate to move vertically
  legend = levels(pop(genotype_diagnostic_genind)),
  col = cols,
  pch = 19,
  bty = "n"   # no box around legend
  )

s.class(pca_2$li, 
  pop(genotype_diagnostic_dog_genind),
  col = cols_2,
  axesell = T,
  cstar = 0,
  cellipse = 0,
  xlim = c(-40, 10),
  ylim = c(-15, 15)
  )

# Add a legend/key on the side
legend(x = -10,  # adjust X coordinate to move horizontally
  y = 10,   # adjust Y coordinate to move vertically
  legend = levels(pop(genotype_diagnostic_dog_genind)), 
  col = cols_2, 
  pch = 19,
  bty = "n"   # no box around legend
  )

s.class(pca_3$li, 
  pop(genotype_diagnostic_coyote_genind),
  col = cols_3,
  axesell = T,
  cstar = 0,
  cellipse = 0,
  xlim = c(-40, 10),
  ylim = c(-15, 15)
  )

# Add a legend/key on the side
legend(x = -10,  # adjust X coordinate to move horizontally
  y = 10,   # adjust Y coordinate to move vertically
  legend = levels(pop(genotype_diagnostic_coyote_genind)), 
  col = cols_3, 
  pch = 19,
  bty = "n"   # no box around legend
  )

#Build a data frame for plotting
pca_df <- pca_1$li %>%
  mutate(Indiv = rownames(pca_1$li),
  Pop   = as.factor(pop(genotype_diagnostic_genind))
  )

# Interactive PCA
p <- plot_ly(data = pca_df,
  x = ~Axis1,
  y = ~Axis2,
  type = "scatter",
  mode = "markers",
  color = ~Pop,
  colors = cols,
  text = ~paste("ID:", Indiv, "<br>Pop:", Pop),
  hoverinfo = "text"
  ) %>%
  layout(title = "Interactive PCA of Genotypes",
  xaxis = list(title = "PC1", range = c(-40, 10)),
  yaxis = list(title = "PC2", range = c(-40, 40)),
  legend = list(title = list(text = "<b>Population</b>"))
  )



# pca in adegenet 
# Build PCA data frame with metadata for plotting

pca_df <- pca_1$li %>%
  mutate(
    Indiv = rownames(pca_1$li),
    Pop   = key[match(rownames(pca_1$li), indNames(genotype_diagnostic_genind))]
  )

# Define colors
n_groups <- length(unique(pca_df$Pop))
cols <- brewer.pal(min(8, n_groups), "Set1")  # adjust palette if more than 8 groups

# Plot PCA without ellipses
ggplot(pca_df, aes(x = Axis1, y = Axis2, color = Pop)) +
  geom_point(size = 3) +
  scale_color_manual(values = cols) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", round(pca_1$eig[1] / sum(pca_1$eig) * 100, 1), "%)"),
    y = paste0("PC2 (", round(pca_1$eig[2] / sum(pca_1$eig) * 100, 1), "%)"),
    color = "Group"
  ) + 
  theme(legend.position = "right")




## filter individuals based on diagnostic pca
# Convert PCA coordinates to data frame
pca1_df <- as.data.frame(pca_1$li) %>%
  mutate(
    Indiv = rownames(pca_1$li),
    Pop   = as.factor(pop(genotype_diagnostic_genind))
  )

# Filter individuals within the specified range
pca1_filtered <- pca1_df %>%
  filter(Axis1 >= -2, Axis1 <= 2,
         Axis2 >= -6, Axis2 <= 6)

# View filtered individuals
keep_indivs <- pca1_filtered$Indiv
removed_sample_3 <- setdiff(indNames(genind_filt1), keep_indivs) # removed 77

# Combine data from all filtering steps
snp_df_filt <- snp_df %>%
  filter(!Locus %in% removed_loci_diagnostic) %>%  # remove diagnostic loci
  filter(!Indiv %in% removed_sample_2a) %>% # remove missing >.2 genotypes all diagnostic
  filter(!Indiv %in% removed_sample_2b) %>% # remove missing >.2 genotypes dog diagnostic
  filter(!Indiv %in% removed_sample_2c) %>% # remove missing >.2 genotypes coyote diagnostic
  filter(!Indiv %in% removed_sample_3) %>% # remove non-wolf samples from diagnostic pca 
  filter(!Indiv %in% removed_sample_1) %>%  # remove missing >.2 genotypes
  filter(!Locus %in% removed_loci_1)  # remove missing >.2 loci

# Filter this df for replicates 
# Identify all unique alleles per locus
alleles_by_locus_filt <- snp_df_filt %>%
  select(Locus, Allele1, Allele2) %>%
  pivot_longer(cols = c(Allele1, Allele2), values_to = "Allele") %>%
  filter(!is.na(Allele)) %>%        
  distinct(Locus, Allele) %>%      
  group_by(Locus) %>%
  mutate(Allele_Index = dense_rank(Allele)) %>%
  ungroup()

# Attach numeric indices and make ordered genotype codes
df_encoded_filt <- snp_df_filt %>%
  left_join(alleles_by_locus,
            by = c("Locus", "Allele1" = "Allele")) %>%
  rename(A1_Code = Allele_Index) %>%
  left_join(alleles_by_locus,
            by = c("Locus", "Allele2" = "Allele")) %>%
  rename(A2_Code = Allele_Index) %>%
  mutate(
    Genotype_Code = ifelse(
      is.na(A1_Code) | is.na(A2_Code),
      NA_character_,
      sprintf("%d%d", A1_Code, A2_Code)  # keep A1/A2 order exactly
    )
  )

# Long to wide: one row per individual
genotype_matrix_filt <- df_encoded_filt %>%
  select(Indiv, Locus, Genotype_Code) %>%
  pivot_wider(
    names_from = Locus,
    values_from = Genotype_Code
  ) %>%
  mutate(
    Indiv_long  = Indiv,
    Indiv_short = str_extract(Indiv, "Clu.*|NTC.*"),
    Indiv_short = if_else(is.na(Indiv_short), Indiv, Indiv_short),
    Indiv_short = str_replace(Indiv_short, "(_initial|_f1).*", "")
  ) %>%
  select(Indiv, Indiv_long, Indiv_short, everything())  # move new columns to front

# Filter for missingness
genotype_matrix_filt$prop_missing <- rowMeans(is.na(genotype_matrix_filt))

# match to short names
filt2_df <- genotype_matrix_filt %>%
  group_by(Indiv_short) %>%
  slice_min(prop_missing, with_ties = F) %>% # Choose replicate with least missing data, NOT allowing ties
  ungroup() #2288 x 284


## look at replicates and calculate error rate 
## Identify replicate samples (triplicates or more)
replicates <- genotype_matrix_filt %>%
  group_by(Indiv_short) %>%
  filter(n() > 2) %>%
  ungroup()

## Create replicate pairs
replicate_map <- replicates %>%
  select(Indiv, Indiv_short) %>%
  group_by(Indiv_short) %>%
  summarise(
    pairs = list(as.data.frame(t(combn(Indiv, 2)))),
    .groups = "drop"
  ) %>%
  unnest(pairs) %>%
  rename(Sample1 = V1, Sample2 = V2)

## Convert genotype matrix to long format
gt_long <- replicates %>%
  select(Indiv, starts_with("Clu_")) %>%
  pivot_longer(
    cols = -Indiv,
    names_to = "Locus",
    values_to = "GT"
  )

## Join genotypes for replicate pairs
paired_gt_df <- replicate_map %>%
  left_join(gt_long, by = c("Sample1" = "Indiv")) %>%
  rename(GT1 = GT) %>%
  left_join(gt_long, by = c("Sample2" = "Indiv", "Locus")) %>%
  rename(GT2 = GT)

## Genotype comparison function
score_gt_match <- function(gt1, gt2){
  
  if (is.na(gt1) & is.na(gt2)) return("BothNA")
  else if (is.na(gt1) | is.na(gt2)) return("OneNA")
  
  gt1 <- as.character(gt1)
  gt2 <- as.character(gt2)
  
  a1 <- sort(strsplit(gt1, "")[[1]])
  a2 <- sort(strsplit(gt2, "")[[1]])
  
  if (identical(a1,a2)) return("Full Match")
  else if (any(a1 %in% a2)) return("Partial Match")
  else return("Mismatch")
}

## Score genotype matches
paired_gt_df <- paired_gt_df %>%
  mutate(MatchType = mapply(score_gt_match, GT1, GT2))

## Calculate pairwise mismatch statistics
pair_errors <- paired_gt_df %>%
  group_by(Sample1, Sample2) %>%
  summarise(
    N = n(),
    Full_Match = sum(MatchType == "Full Match", na.rm = TRUE),
    Partial_Match = sum(MatchType == "Partial Match", na.rm = TRUE),
    Mismatch = sum(MatchType == "Mismatch", na.rm = TRUE),
    Missing_Full = sum(MatchType == "BothNA", na.rm = TRUE),
    Missing_Partial = sum(MatchType == "OneNA", na.rm = TRUE),
    Full_Match_Rate = Full_Match / (Full_Match + Partial_Match + Mismatch),
    Mismatch_Rate = (Mismatch + Partial_Match) / (Full_Match + Partial_Match + Mismatch),
    Inclusive_Missingness = (Missing_Full + Missing_Partial) / N,
    .groups = "drop"
  )

## Extract sample IDs if needed
pair_errors$ID_1 <- sub(".*_", "", pair_errors$Sample1)
pair_errors$ID_2 <- sub(".*_", "", pair_errors$Sample2)

## Summarize mismatch across replicate pairs
pair_err_summ <- pair_errors %>%
  summarise(
    N_Samples = n(),
    Avg_Mismatch = round(mean(Mismatch_Rate, na.rm = TRUE), 4),
    Max_Mismatch = round(max(Mismatch_Rate, na.rm = TRUE), 4),
    Min_Locus_Overlap = min(((1 - Inclusive_Missingness) * N), na.rm = TRUE),
    Avg_Missing = round(mean(Inclusive_Missingness * N, na.rm = TRUE), 1)
  )

#write.csv(paired_gt_df, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/output/paired_gt_df.csv", row.names = F)

overall_error_rate <- sum(pair_errors$Mismatch + pair_errors$Partial_Match) /
  sum(pair_errors$Full_Match +
        pair_errors$Partial_Match +
        pair_errors$Mismatch)


# now filter these for relatedness
# Extract ID column
indiv <- filt2_df$Indiv_short

# Extract only genotype loci (drop ID columns)
loci <- filt2_df %>% 
  select(-Indiv, -Indiv_long, -Indiv_short)

# Split each locus into two allele columns
split_cols <- lapply(names(loci), function(locus_name) {
  x <- as.character(loci[[locus_name]])
  
  allele1 <- ifelse(is.na(x), NA, substr(x, 1, 1))
  allele2 <- ifelse(is.na(x), NA, substr(x, 2, 2))
  
  df_locus <- data.frame(allele1, allele2, stringsAsFactors = FALSE)
  
  names(df_locus) <- paste0(locus_name, c("_allele1", "_allele2"))
  
  df_locus
})

# Combine into a single data frame
geno_split <- do.call(cbind, split_cols)

# Add Indiv column back
geno_split <- cbind(Indiv = indiv, geno_split)

# Ensure allele columns are character
geno_split_2 <- geno_split %>% #2288 x 561
  mutate(across(-Indiv, as.character)) %>%       # make genotype columns character
  select(-prop_missing_allele1, -prop_missing_allele2) %>% 
  # filter(!if_all(2:561, is.na)) %>%              # remove rows where all genotypes are NA
  mutate(across(-Indiv, ~ replace_na(.x, "0"))) # replace remaining NAs with "0" - this is required to run coancestry

#########################################################################################
# Run (must use UI server, shell code to submit r script as job)
related_trioml <- coancestry(geno_split_2, trioml = 1)

# df for results
related_trioml_df <- related_trioml$relatedness

# Save results as csv
write.csv(related_trioml_df, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/output/relatedness_trioml_full.csv", row.names = F)

###########################################################################################
# Load data that ran on the server 
relatedness <- read.csv("/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/output/relatedness_trioml_full.csv") %>%
  dplyr::rename(
    ind1 = ind1.id,
    ind2 = ind2.id,
    relatedness = trioml
  )
head(relatedness)
# Filter for relatedness
# Related indv
related_0.4 <- relatedness %>% 
  filter(relatedness > 0.4) %>% 
  arrange(desc(ind1))

unique(related_0.4$ind1) %>% length()
unique(related_0.4$ind2) %>% length()

df_filt <- related_0.4 %>% 
  select(ind1, ind2, relatedness)

g <- graph_from_data_frame(df_filt[, c("ind1", "ind2")], directed = FALSE)

comp <- components(g)    
related_groups <- split(V(g)$name, comp$membership)

head(related_groups)

# Function to select individuals w most related connections 
break_clusters <- function(df, threshold = 0.4, verbose = TRUE) {
  # Filter by relatedness (optional depending on your threshold)
  df_filt <- df[df$relatedness >= threshold, ]
  removed_ids <- c()
  iteration <- 1
  repeat {
    if (verbose) cat("\nIteration", iteration, "\n")
    # Build the graph
    g <- graph_from_data_frame(df_filt[, c("ind1", "ind2")], directed = FALSE)
    # Identify connected components
    comp <- components(g)
    groups <- split(V(g)$name, comp$membership)
    # Identify groups larger than 1 (i.e., real clusters)
    large_groups <- groups[sapply(groups, length) > 1]
    if (length(large_groups) == 0) {
      if (verbose) cat("No clusters left. Stopping.\n")
      break
    }
    if (verbose) cat("Found", length(large_groups), "clusters.\n")
    # Compute degree of each vertex
    deg <- degree(g)
    # Identify the top-degree individual *in each cluster*
    to_remove <- sapply(large_groups, function(grp_ids) {
      grp_ids[which.max(deg[grp_ids])]
    })
    if (verbose) cat("Removing:", paste(to_remove, collapse=", "), "\n")
    # Store removed IDs
    removed_ids <- c(removed_ids, to_remove)
    # Remove from dataframe
    df_filt <- subset(df_filt, !(ind1 %in% to_remove | ind2 %in% to_remove))
    iteration <- iteration + 1
  }
  return(list(
    cleaned_df = df_filt,
    removed_ids = removed_ids
  ))
}

df <- relatedness %>% 
  select(ind1, ind2, relatedness)

result <- break_clusters(df, threshold = 0.4)

removed_individuals <- result$removed_ids
head(removed_individuals)
length(removed_individuals)


# Combine all individuals from both columns in the relatedness df
all_indivs <- df %>%
  select(ind1, ind2) %>%
  unlist() %>%
  .[!duplicated(.)]

kept_indivs <- setdiff(all_indivs, removed_individuals)

# View the result
kept_indivs

# Cross check that removing those individuals actually leads to no first order pairs
df_clean <- df %>%
  filter(!ind1 %in% removed_individuals,
         !ind2 %in% removed_individuals)

df_clean %>% arrange(desc(relatedness))

length(unique(c(df_clean$ind1, df_clean$ind2)))
length(unique(c(df$ind1, df$ind2)))

df_0.4 <- df %>% filter(relatedness > 0.4)
length(unique(df_0.4$ind1))

length(removed_individuals)

# df of indivs filtered for relatedness
kept_indivs_df <- data.frame(Indiv = kept_indivs, stringsAsFactors = FALSE) 



#######################################################################################################################
# Keep only ones with metadata
metadata_structure <- read.csv("/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/sample_metadata_DAU.csv")
kept_indivs_df$Indiv
metadata_structure

# Pull only metadata rows where Sample matches $Indiv
metadata_matched2 <- metadata_structure %>% 
  filter(Sample %in% kept_indivs_df$Indiv)

duplicated_samples2 <- metadata_matched2 %>%
  count(Sample) %>%
  filter(n > 1) %>%
  pull(Sample)

# Step 2: Pull all rows from metadata_matched that are duplicates
metadata_duplicates2 <- metadata_matched2 %>%
  filter(Sample %in% duplicated_samples2)

# Step 1: Remove duplicates in metadata, keeping the first occurrence
metadata_unique2 <- metadata_matched2 %>%
  group_by(Sample) %>%
  slice(1) %>%    # keep first row for duplicates
  ungroup()

# Match the order of structure_data$Indiv
metadata_ordered2 <- kept_indivs_df %>%
  select(Indiv) %>%
  left_join(metadata_unique2, by = c("Indiv" = "Sample"))

# Filtered df with no missing metadata
structure_clean <- metadata_ordered2 %>%
  filter(if_all(-1, ~ !is.na(.)))

# Format for structure 
structure_data_df <- geno_split_2 %>%
  filter(str_detect(Indiv, paste0(structure_clean$Indiv, collapse = "|"))) 

# Get all loci base names (strip allele1/allele2)
loci <- unique(str_remove(names(structure_data_df)[grepl("allele", names(structure_data_df))],
                          "_allele1|_allele2"))

# For each locus, sort allele1 > allele2
for(locus in loci) {
  a1 <- paste0(locus, "_allele1")
  a2 <- paste0(locus, "_allele2")
  
  structure_data_df[[a1]] <- pmax(structure_data_df[[a1]], structure_data_df[[a2]], na.rm = TRUE)
  structure_data_df[[a2]] <- pmin(structure_data_df[[a1]], structure_data_df[[a2]], na.rm = TRUE)
}

# Final format for STRUCTURE 
structure_data_df3 <- structure_data_df %>%
  mutate(Indiv = str_remove(Indiv, "^CluIDFG")) %>%        # remove prefix
  mutate(Indiv = str_replace_all(Indiv, "_", ".")) %>%     # replace underscores
  mutate(Indiv = str_remove(Indiv, "\\.(initial|f1).*$"))  # remove .initial or .f1 + anything after

# Save csv with final data txt file to use in STRUCTURE analysis
write.table(
  structure_data_df3,
  file = "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/STRUCTURE/structure_data_3.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
