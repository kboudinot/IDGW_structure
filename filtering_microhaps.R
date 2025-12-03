# Load packages
library(tidyverse)
library(adegenet) 
library(stringr)
library(poppr)
library(ade4)
library(RColorBrewer)
library(plotly)
library(related) 
library(igraph)

# Load data
snp_raw <- read.csv("/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/microhaplotype_genotypes_raw.csv", na.strings = "")

# Clean raw data
snp_df <- snp_raw %>%
  mutate(Indiv = Indiv |>              
      basename() |>                      # remove full path
      sub("\\.1\\.bam$", "", x = _)      # remove trailing ".1.bam"
  )

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

# Make the dataframe a genind object
snp_genind <- df2genind(
  genotype_matrix %>% select(-Indiv, -Indiv_short, -Indiv_long), 
  ind.names = genotype_matrix$Indiv,   
  ploidy = 2,
  sep = "",
  NA.char = "NA"
)

# Confirm names 
# indNames(snp_genind) <- indNames(snp_genind) # long names
# #%>% str_extract("Clu.*|NTC.*")

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

# Vector of samples removed
removed_loci_1 <- setdiff(locNames(snp_genind), locNames(genind_loci_filt))
removed_sample_1 <- setdiff(indNames(snp_genind), indNames(genind_ind_filt))

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
          
          # # ggplot
          # ggplot(pca_scores, aes(x = Axis1, y = Axis2, color = missing)) +
          #   geom_point(size = 1, alpha = 0.9) +
          #   scale_color_gradient(
          #     low = "#01579b",  # dark blue
          #     high = "#b3e5fc", # light blue
          #     name = "Missingness"
          #   ) +
          #   labs(
          #     x = paste0("PC1 (", var_expl[1], "%)"),
          #     y = paste0("PC2 (", var_expl[2], "%)"),
          #     title = "PCA of SNP Data",
          #     subtitle = "Colored by genotype missingness"
          #   ) +
          #   theme_bw(base_size = 14) +
          #   theme(
          #     plot.title = element_text(face = "bold"),
          #     panel.grid = element_blank()
          #   )

#genind_filt1 # missingness filtered genind

          
          
# Filter diagnostic and missingness
#snp_df

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

# shorten names
#indNames(genotype_diagnostic_genind) <- str_extract(indNames(genotype_diagnostic_genind), "Clu.*|NTC.*") %>%
 # str_replace("\\.1\\.bam$", "")
#indNames(genotype_diagnostic_dog_genind) <- str_extract(indNames(genotype_diagnostic_dog_genind), "Clu.*|NTC.*") %>%
 # str_replace("\\.1\\.bam$", "")
#indNames(genotype_diagnostic_coyote_genind) <- str_extract(indNames(genotype_diagnostic_coyote_genind), "Clu.*|NTC.*") %>%
 # str_replace("\\.1\\.bam$", "")

#$ Making diagnostic pca
          # import metadata
          metadata <- read.csv("/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/sample_metadata_regions.csv")
          # match to genind
      #    head(indNames(genotype_diagnostic_genind))
    #      head(locNames(genotype_diagnostic_genind))
     #     head(metadata$Sample) # need to add long sample names back 
          
          geno_df <- as.data.frame(genotype_diagnostic_genind$tab)
          geno_df <- geno_df %>%
            mutate(
              Sample_long = rownames(.),                            # keep the full name
              Sample      = str_extract(Sample_long, "Clu.*|NTC.*") # extract short ID
            )
       #   head(geno_df)
          
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
           s.class(
            pca_1$li,
            pop(genotype_diagnostic_genind),
            col = cols,
            axesell = T,
            cstar = 0,
            cellipse = 0,
            xlim = c(-40, 10),
            ylim = c(-20, 20)
          )

          # Add a legend/key on the side
          legend(
            x = -15,  # adjust X coordinate to move horizontally
            y = 15,   # adjust Y coordinate to move vertically
            legend = levels(pop(genotype_diagnostic_genind)),
            col = cols,
            pch = 19,
            bty = "n"   # no box around legend
          )

          # s.class(
          #   pca_2$li, 
          #   pop(genotype_diagnostic_dog_genind),
          #   col = cols_2,
          #   axesell = T,
          #   cstar = 0,
          #   cellipse = 0,
          #   xlim = c(-40, 10),
          #   ylim = c(-15, 15)
          # )
          # # Add a legend/key on the side
          # legend(
          #   x = -10,  # adjust X coordinate to move horizontally
          #   y = 10,   # adjust Y coordinate to move vertically
          #   legend = levels(pop(genotype_diagnostic_dog_genind)), 
          #   col = cols_2, 
          #   pch = 19,
          #   bty = "n"   # no box around legend
          # )
          # 
          # s.class(
          #   pca_3$li, 
          #   pop(genotype_diagnostic_coyote_genind),
          #   col = cols_3,
          #   axesell = T,
          #   cstar = 0,
          #   cellipse = 0,
          #   xlim = c(-40, 10),
          #   ylim = c(-15, 15)
          # )
          # # Add a legend/key on the side
          # legend(
          #   x = -10,  # adjust X coordinate to move horizontally
          #   y = 10,   # adjust Y coordinate to move vertically
          #   legend = levels(pop(genotype_diagnostic_coyote_genind)), 
          #   col = cols_3, 
          #   pch = 19,
          #   bty = "n"   # no box around legend
          # )
          
          # 
          # 
          # Build a data frame for plotting
          pca_df <- pca_1$li %>%
            mutate(
              Indiv = rownames(pca_1$li),
              Pop   = as.factor(pop(genotype_diagnostic_genind))
            )
          # Interactive PCA
          p <- plot_ly(
            data = pca_df,
            x = ~Axis1,
            y = ~Axis2,
            type = "scatter",
            mode = "markers",
            color = ~Pop,
            colors = cols,
            text = ~paste("ID:", Indiv, "<br>Pop:", Pop),
            hoverinfo = "text"
          ) %>%
            layout(
              title = "Interactive PCA of Genotypes",
              xaxis = list(title = "PC1", range = c(-40, 10)),
              yaxis = list(title = "PC2", range = c(-40, 40)),
              legend = list(title = list(text = "<b>Population</b>"))
            )
          p
          # 
          # # Build a data frame for pca_2
          # pca2_df <- pca_2$li %>%
          #   as.data.frame() %>%
          #   mutate(
          #     Indiv = rownames(pca_2$li),
          #     Pop   = as.factor(pop(genotype_diagnostic_dog_genind))
          #   )
          # 
          # # Interactive PCA plot for pca_2
          # p2 <- plot_ly(
          #   data = pca2_df,
          #   x = ~Axis1,
          #   y = ~Axis2,
          #   type = "scatter",
          #   mode = "markers",
          #   color = ~Pop,
          #   colors = cols_2,
          #   text = ~paste("ID:", Indiv, "<br>Pop:", Pop),
          #   hoverinfo = "text"
          # ) %>%
          #   layout(
          #     title = "Interactive PCA - Dog Genotypes",
          #     xaxis = list(title = "PC1", range = c(-40, 10)),
          #     yaxis = list(title = "PC2", range = c(-20, 20)),
          #     legend = list(title = list(text = "<b>Population</b>"))
          #   )
          # p2
          # 
          # # Build a data frame for pca_3
          # pca3_df <- pca_3$li %>%
          #   as.data.frame() %>%
          #   mutate(
          #     Indiv = rownames(pca_3$li),
          #     Pop   = as.factor(pop(genotype_diagnostic_coyote_genind))
          #   )
          # 
          # # Interactive PCA plot for pca_3
          # p3 <- plot_ly(
          #   data = pca3_df,
          #   x = ~Axis1,
          #   y = ~Axis2,
          #   type = "scatter",
          #   mode = "markers",
          #   color = ~Pop,
          #   colors = cols_3,
          #   text = ~paste("ID:", Indiv, "<br>Pop:", Pop),
          #   hoverinfo = "text"
          # ) %>%
          #   layout(
          #     title = "Interactive PCA - Coyote Genotypes",
          #     xaxis = list(title = "PC1", range = c(-40, 10)),
          #     yaxis = list(title = "PC2", range = c(-50, 50)),
          #     legend = list(title = list(text = "<b>Population</b>"))
          #   )
          # p3

## filter individuals based on diagnostic pca----
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
#writeLines(removed_sample_3, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/removed_sample_3.txt")

# Combine data from all filtering steps
snp_df_filt <- snp_df %>%
  filter(!Locus %in% removed_loci_diagnostic) %>%  # remove diagnostic loci
  filter(!Indiv %in% removed_sample_2a) %>% # remove missing >.2 genotypes all diagnostic
  filter(!Indiv %in% removed_sample_2b) %>% # remove missing >.2 genotypes dog diagnostic
  filter(!Indiv %in% removed_sample_2c) %>% # remove missing >.2 genotypes coyote diagnostic
  filter(!Indiv %in% removed_sample_3) %>% # remove non-wolf samples from diagnostic pca 
  filter(!Indiv %in% removed_sample_1) %>%  # remove missing >.2 genotypes
  filter(!Locus %in% removed_loci_1) %>%  # remove missing >.2 loci

length(unique(snp_df_filt$Indiv)) #2614

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
    Indiv_short = if_else(is.na(Indiv_short), Indiv, Indiv_short)
  ) %>%
  select(Indiv, Indiv_long, Indiv_short, everything())  # move new columns to front

#head(genotype_matrix_filt)

# Filter for missingness
genotype_matrix_filt$prop_missing <- rowMeans(is.na(genotype_matrix_filt))

# match to short names
filt2_df <- genotype_matrix_filt %>%
  group_by(Indiv_short) %>%
  slice_min(prop_missing, with_ties = F) %>% # Choose replicate with least missing data, allowing ties
  ungroup()

# unique_loci_filt2_df <- filt2_df %>%
#   select(Locus) %>%
#   distinct() %>% 
#   pull(Locus)

# removed_sample_4 <- setdiff(snp_df$Indiv, filt2_df$Indiv)
# #writeLines(removed_sample_4, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/removed_sample_4.txt")
# 
# filtering_before_relatedness <- snp_df %>%
#   filter(!Locus %in% removed_loci_diagnostic) %>%  # remove diagnostic loci
#   filter(!Indiv %in% removed_sample_2a) %>% # remove missing >.2 genotypes all diagnostic
#   filter(!Indiv %in% removed_sample_2b) %>% # remove missing >.2 genotypes dog diagnostic
#   filter(!Indiv %in% removed_sample_2c) %>% # remove missing >.2 genotypes coyote diagnostic
#   filter(!Indiv %in% removed_sample_3) %>% # remove non-wolf samples from diagnostic pca 
#   filter(!Indiv %in% removed_sample_1) %>%  # remove missing >.2 genotypes
#   filter(!Locus %in% removed_loci_1) %>%  # remove missing >.2 loci
#   filter(!Indiv %in% removed_sample_4) # remove relatedness

# unique_indivs <- filtering_before_relatedness %>%
#   select(Indiv) %>%
#   distinct() %>% 
#   pull(Indiv)
#writeLines(unique_indivs, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/unique_indivs.txt")

# unique_loci <- filtering_before_relatedness %>%
#   select(Locus) %>%
#   distinct() %>% 
#   pull(Locus)
#writeLines(unique_loci, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/unique_loci.txt")

# now filter these for relatedness
#filt2_df

# split alleles in df
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
geno_split <- geno_split %>% 
  mutate(across(-Indiv, as.character)) %>% 
  select(-prop_missing_allele1, -prop_missing_allele2) %>%
  filter(if_all(-Indiv, ~ !is.na(.)))


# Run (must use UI server, shell code to submit r script as job)
related <- coancestry(geno_split, trioml=1, wang = 1)
# related_trioml <- coancestry(geno_split, trioml = 1)

# # df for results
# related_df <- related$relatedness
# related_trioml_df <- related_trioml$relatedness

# Save results as csv
#write.csv(related_df, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/output/relatedness.csv", row.names = F)
write.csv(related_df, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/output/relatedness_trioml.csv", row.names = F)


# Load data that ran on the server 
relatedness <- read.csv("/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/output/relatedness_trioml.csv")

# Make a data frame of relatedness values
trioml_df <- data.frame(trioml = relatedness$trioml)

# Percent above 0.4
trioml_df %>% 
  summarize(percent_above_0.4 = mean(trioml > 0.4) * 100)

# Plot histogram
ggplot(trioml_df, aes(x = trioml)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "white") +
  labs(title = "Histogram of Trioml Relatedness", x = "Relatedness", y = "Count") +
  theme_minimal()



# Filter for relatedness

# Related indv
related_0.4 <- relatedness %>% filter(trioml > 0.4) %>% arrange(desc(ind1.id))

unique(related_0.4$ind1.id) %>% length()
unique(related_0.4$ind2.id) %>% length()


df_filt <- related_0.4 %>% select("ind1.id", "ind2.id", "trioml")

g <- graph_from_data_frame(df_filt[, c("ind1.id", "ind2.id")], directed = FALSE)
comp <- components(g)    # or components(g) in newer igraph

related_groups <- split(V(g)$name, comp$membership)
head(related_groups)


#############################################################
# Function to select individuals w most related connections #
#############################################################
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

df <- relatedness %>% select("ind1.id", "ind2.id", "trioml") #trio
colnames(df) <- c("ind1", "ind2", "relatedness")

result <- break_clusters(df, threshold = 0.4)

clean_df <- result$cleaned_df
removed_individuals <- result$removed_ids
head(removed_individuals)
length(removed_individuals)

# Combine all individuals from both columns in the relatedness df
all_indivs <- unique(c(relatedness$ind1.id, relatedness$ind2.id))

# Keep only those NOT in removed_individuals
kept_indivs <- setdiff(all_indivs, removed_individuals)

# View the result
kept_indivs
writeLines(kept_indivs, "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/kept_indivs.txt")




### Cross check that removing those individuals actually leads to no first order pairs
length(df$ind1)
df_ind1_rm <- df %>% filter(!ind1 %in% removed_individuals)
length(df_ind1_rm$ind1)
df_ind2_rm <- df_ind1_rm %>% filter(!ind2 %in% removed_individuals)
length(df_ind2_rm$ind1)

df_ind2_rm %>% arrange(desc(relatedness))

unique(c(df_ind2_rm$ind1, df_ind2_rm$ind2)) %>% length()

unique(c(df$ind1, df$ind2)) %>% length()

df_0.4 <- df %>% filter(relatedness > 0.4) 

unique(df_0.4$ind1) %>% length()
length(removed_individuals)


## Final samples for STRUCTURE
kept_indivs_df <- data.frame(Indiv = kept_indivs, stringsAsFactors = FALSE) #511
# kept_indivs_clean <- kept_indivs_df %>%
#   mutate(Indiv = gsub("(_init|_f1_2)$", "", Indiv))

structure_data <- geno_split %>%
  filter(str_detect(Indiv, paste0(kept_indivs_df$Indiv, collapse = "|"))) #511

write.csv(structure_data, 
          "/mnt/ceph/kboudinot/wolves/data/neutral/R_files/pop_gen/data/structure_data.csv", 
          row.names = FALSE)
