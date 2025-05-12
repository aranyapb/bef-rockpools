
# Creating 6 raw datasets to use in data analysis - 

# alpha_all_raw, alpha_active_raw, alpha_passive_raw - Alpha scale (individual 
# pools within an inselberg) data for all taxa, only active dispersers and only 
# passive dispersers respectively. 
# Includes the following columns (for each pool) - Inselberg, Pool, Alpha richness,
# Alpha standard error (ASE), Gamma diversity, Biomass, Depth, Gamma, 
# Gamma standard error (GSE), Climate PC1 and Climate PC2. 

# gamma_all_raw, gamma_active_raw, gamma_passive raw - Gamma scale (between 
# inselbergs) data for all taxa, only active dispersers and only 
# passive dispersers respectively. 
# Includes the following columns (for each inselberg) - Inselberg, Pool, 
# Gamma diversity, Median biomass, Median Depth, Mean log biomass, Mean log depth, 
# Gamma, Gamma standard error (GSE),Prominence, Climate PC1 and Climate PC2. 

# Note: Run code in order. 

#-------------------------------------------------------------------------------
# COMMUNITY BIOMASS 
#-------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(readr)

# ALL TAXA 

BM_df <- read_csv("Data-raw/pool-community-biomass-data.csv")

BM_df %>% ncol() # All taxa = 436
BM_df %>% nrow() # 241 rows 

# pools to remove - App14 has a missing depth, Kam11 is too deep and Kor13 is eutrophicated (most likely). 
remove_pools <- c("App14", "Kam11", "Kor13")
BM_df <- BM_df %>% filter(!Pool %in% remove_pools) 

# Calculating total community biomass for each pool 
BM_df_1 <- BM_df %>%
  rowwise() %>%
  mutate(Biomass = sum(c_across(-c(Inselberg, Pool)), na.rm = TRUE)) %>%
  ungroup() %>% 
  dplyr::select(Inselberg, Pool, Biomass)


# Calculate median biomass and mean log biomass for each inselberg
Biomass_gamma <- BM_df_1 %>%
  dplyr::select(Inselberg, Biomass) %>% 
  mutate(Inselberg = factor(Inselberg, levels = unique(Inselberg))) %>%
  group_by(Inselberg) %>%
  summarise(Biomass_median = median(Biomass, na.rm = TRUE),
            Meanlogbiomass = mean(log(Biomass+1), na.rm = TRUE))

#-------------------------------------------------------------------------------

# SEPPARATING ACTIVE AND PASSIVE DISPERSERS  
# Load data
AP_df <- read_csv("Data-raw/active-passive-split-taxa.csv")

BM_df <- read_csv("Data-raw/pool-community-biomass-data.csv")

# Define labels and taxa
labels_cols <- c("Inselberg", "Pool")
all_taxa <- setdiff(names(BM_df), labels_cols)

# Classify taxa
active_cols <- intersect(all_taxa, AP_df$Active)
passive_cols <- intersect(all_taxa, AP_df$Passive)

# Check for unclassified taxa
unclassified <- setdiff(all_taxa, c(active_cols, passive_cols))
unclassified

# Subset data
BM_df_Active <- BM_df %>% dplyr::select(all_of(labels_cols), all_of(active_cols))
BM_df_Passive <- BM_df %>% dplyr::select(all_of(labels_cols), all_of(passive_cols))

#-------------------------------------------------------------------------------

# BIOMASS FOR ACTIVE DISPERSERS 

# pools to remove - App14 has a missing depth, Kam11 is too deep and Kor13 is eutrophicated (most likely). 
remove_pools <- c("App14", "Kam11", "Kor13")
BM_df_Active <- BM_df_Active %>% filter(!Pool %in% remove_pools) 

# Calculating total community biomass for each pool 
BM_df_Active_1 <- BM_df_Active %>%
  rowwise() %>%
  mutate(Biomass = sum(c_across(-c(Inselberg, Pool)), na.rm = TRUE)) %>%
  ungroup() %>% 
  dplyr::select(Inselberg, Pool, Biomass)


# Calculate median biomass and mean log biomass for each inselberg
Biomass_gamma_Active <- BM_df_Active_1 %>%
  dplyr::select(Inselberg, Biomass) %>% 
  mutate(Inselberg = factor(Inselberg, levels = unique(Inselberg))) %>%
  group_by(Inselberg) %>%
  summarise(Biomass_median = median(Biomass, na.rm = TRUE),
            Meanlogbiomass = mean(log(Biomass+1), na.rm = TRUE))

#-------------------------------------------------------------------------------

# BIOMASS FOR PASSIVE DISPERSERS

# pools to remove - App14 has a missing depth, Kam11 is too deep and Kor13 is eutrophicated (most likely). 
remove_pools <- c("App14", "Kam11", "Kor13")
BM_df_Passive <- BM_df_Passive %>% filter(!Pool %in% remove_pools) 

# Calculating total community biomass for each pool 
BM_df_Passive_1 <- BM_df_Passive %>%
  rowwise() %>%
  mutate(Biomass = sum(c_across(-c(Inselberg, Pool)), na.rm = TRUE)) %>%
  ungroup() %>% 
  dplyr::select(Inselberg, Pool, Biomass)


# Calculate median biomass and mean log biomass for each inselberg
Biomass_gamma_Passive <- BM_df_Passive_1 %>%
  dplyr::select(Inselberg, Biomass) %>% 
  mutate(Inselberg = factor(Inselberg, levels = unique(Inselberg))) %>%
  group_by(Inselberg) %>%
  summarise(Biomass_median = median(Biomass, na.rm = TRUE),
            Meanlogbiomass = mean(log(Biomass+1), na.rm = TRUE))

#-------------------------------------------------------------------------------
# ALPHA RICHNESS 
#-------------------------------------------------------------------------------

library(iNEXT)
library(dplyr)

ab_df <- read_csv("Data-raw/pool-community-abundance-data.csv")

ab_df %>% ncol() # All taxa = 436 (excluding Inselberg and Pool column)
ab_df %>% nrow() # 241 rows 

# pools to remove - App14 has a missing depth, Kam11 is too deep and Kor13 is most likely too eutrophicated. 
remove_pools <- c("App14", "Kam11", "Kor13")

ab_df <- ab_df %>% 
  filter(!Pool %in% remove_pools) %>% 
  mutate(across(-c(Inselberg, Pool), round)) # Some values are decimals because of the way abundance was calculated 

# Get inselberg names as a vector 
row_df <- unique(ab_df$Inselberg)

# Create community datasets for each inselberg 
community <- lapply(row_df, function(x) {
  ab_df %>% 
    filter(Inselberg == x) %>%     
    select(-Inselberg) %>% 
    column_to_rownames("Pool")
})

# Calculate Chao alpha richness for each pool in the inselbergs 
alpha_richness <- lapply(seq_along(community), function(i) {
  inselberg_name <- row_df[i]
  comm <- community[[i]]
  pool_names <- rownames(comm)
  
  # Identify zero-sum and non-zero-sum pools
  row_sums <- rowSums(comm)
  non_zero_comm <- comm[row_sums > 0, , drop = FALSE]
  zero_pools <- names(row_sums[row_sums == 0])
  
  # If there's at least one non-zero pool, run iNEXT
  if (nrow(non_zero_comm) > 0) {
    alpha <- iNEXT(
      as.list(as.data.frame(t(non_zero_comm))),
      datatype = "abundance"
    )
    
    out <- alpha$AsyEst %>%
      filter(Diversity == "Species richness") %>% 
      select(Assemblage, Estimator, `s.e.`)
    
      } else {
    out <- tibble(Assemblage = character(), Estimator = numeric(), `s.e.` = numeric())
  }
  
  # Add back zero-pools with 0 values
  if (length(zero_pools) > 0) {
    zero_df <- tibble(
      Assemblage = zero_pools,
      Estimator = 0,
      `s.e.` = 0
    )
    out <- bind_rows(out, zero_df)
  }
  
  out %>%
    mutate(Inselberg = inselberg_name) %>%
    relocate(Inselberg) %>%
    arrange(as.numeric(gsub("\\D+", "", Assemblage)))
})

alpha_df <- bind_rows(alpha_richness) 

alpha_df <- alpha_df %>%
  rename(
    Pool = Assemblage,
    Alpha = Estimator,
    ASE = 's.e.')

#-------------------------------------------------------------------------------
# SEPPARATING ACTIVE AND PASSIVE DISPERSERS  

# Load data
AP_df <- read_csv("Data-raw/active-passive-split-taxa.csv")
ab_df <- read_csv("Data-raw/pool-community-abundance-data.csv")

# Define labels and taxa
labels_cols <- c("Inselberg", "Pool")
all_taxa <- setdiff(names(ab_df), labels_cols)

# Classify taxa
active_cols <- intersect(all_taxa, AP_df$Active)
passive_cols <- intersect(all_taxa, AP_df$Passive)

# Check for unclassified taxa
unclassified <- setdiff(all_taxa, c(active_cols, passive_cols))
unclassified

# Subset data
ab_df_Active <- ab_df %>% dplyr::select(all_of(labels_cols), all_of(active_cols))
ab_df_Passive <- ab_df %>% dplyr::select(all_of(labels_cols), all_of(passive_cols))

#-------------------------------------------------------------------------------
# ALPHA RICHNESS FOR ACTIVE DISPERSERS 

# pools to remove - App14 has a missing depth, Kam11 is too deep and Kor13 is most likely too eutrophicated. 
remove_pools <- c("App14", "Kam11", "Kor13")

ab_df_Active <- ab_df_Active %>% 
  filter(!Pool %in% remove_pools) %>% 
  mutate(across(-c(Inselberg, Pool), round)) # Some values are decimals because of the way abundance was calculated 

# Get inselberg names as a vector 
row_df <- unique(ab_df_Active$Inselberg)

# Create community datasets for each inselberg 
community_Active <- lapply(row_df, function(x) {
  ab_df_Active %>% 
    filter(Inselberg == x) %>%     
    select(-Inselberg) %>% 
    column_to_rownames("Pool")
})


# Calculate Chao alpha richness for each pool in the inselbergs 
alpha_richness_Active <- lapply(seq_along(community_Active), function(i) {
  inselberg_name <- row_df[i]
  comm <- community_Active[[i]]
  pool_names <- rownames(comm)
  
  # Identify zero-sum and non-zero-sum pools
  row_sums <- rowSums(comm)
  non_zero_comm <- comm[row_sums > 0, , drop = FALSE]
  zero_pools <- names(row_sums[row_sums == 0])
  
  # If there's at least one non-zero pool, run iNEXT
  if (nrow(non_zero_comm) > 0) {
    alpha <- iNEXT(
      as.list(as.data.frame(t(non_zero_comm))),
      datatype = "abundance"
    )
    
    out <- alpha$AsyEst %>%
      filter(Diversity == "Species richness") %>% 
      select(Assemblage, Estimator, `s.e.`)
  } else {
    out <- tibble(Assemblage = character(), Estimator = numeric(), `s.e.` = numeric())
  }
  
  # Add back zero-pools with 0 values
  if (length(zero_pools) > 0) {
    zero_df <- tibble(
      Assemblage = zero_pools,
      Estimator = 0,
      `s.e.` = 0
    )
    out <- bind_rows(out, zero_df)
  }
  
  out %>%
    mutate(Inselberg = inselberg_name) %>%
    relocate(Inselberg) %>%
    arrange(as.numeric(gsub("\\D+", "", Assemblage)))
})

alpha_df_Active <- bind_rows(alpha_richness_Active) 

alpha_df_Active <- alpha_df_Active %>%
  rename(
    Pool = Assemblage,
    Alpha = Estimator,
    ASE = 's.e.')

# Lots of warning messages - some sites have only one species - estimation not robust 

#-------------------------------------------------------------------------------
# ALPHA RICHNESS FOR PASSIVE DISPERSERS 

# pools to remove - App14 has a missing depth, Kam11 is too deep and Kor13 is most likely too eutrophicated. 
remove_pools <- c("App14", "Kam11", "Kor13")

ab_df_Passive <- ab_df_Passive %>% 
  filter(!Pool %in% remove_pools) %>% 
  mutate(across(-c(Inselberg, Pool), round)) # Some values are decimals because of the way abundance was calculated 

# Get inselberg names as a vector 
row_df <- unique(ab_df_Passive$Inselberg)

# Create community datasets for each inselberg 
community_Passive <- lapply(row_df, function(x) {
  ab_df_Passive %>% 
    filter(Inselberg == x) %>%     
    select(-Inselberg) %>% 
    column_to_rownames("Pool")
})


# Calculate Chao alpha richness for each pool in the inselbergs 
alpha_richness_Passive <- lapply(seq_along(community_Passive), function(i) {
  inselberg_name <- row_df[i]
  comm <- community_Passive[[i]]
  pool_names <- rownames(comm)
  
  # Identify zero-sum and non-zero-sum pools
  row_sums <- rowSums(comm)
  non_zero_comm <- comm[row_sums > 0, , drop = FALSE]
  zero_pools <- names(row_sums[row_sums == 0])
  
  # If there's at least one non-zero pool, run iNEXT
  if (nrow(non_zero_comm) > 0) {
    alpha <- iNEXT(
      as.list(as.data.frame(t(non_zero_comm))),
      datatype = "abundance"
    )
    
    out <- alpha$AsyEst %>%
      filter(Diversity == "Species richness") %>%
      select(Assemblage, Estimator, `s.e.`)
  } else {
    out <- tibble(Assemblage = character(), Estimator = numeric(), `s.e.` = numeric())
  }
  
  # Add back zero-pools with 0 values
  if (length(zero_pools) > 0) {
    zero_df <- tibble(
      Assemblage = zero_pools,
      Estimator = 0,
      `s.e.` = 0
    )
    out <- bind_rows(out, zero_df)
  }
  
  out %>%
    mutate(Inselberg = inselberg_name) %>%
    relocate(Inselberg) %>%
    arrange(as.numeric(gsub("\\D+", "", Assemblage)))
})

alpha_df_Passive <- bind_rows(alpha_richness_Passive) 

alpha_df_Passive <- alpha_df_Passive %>%
  rename(
    Pool = Assemblage,
    Alpha = Estimator,
    ASE = 's.e.')


#-------------------------------------------------------------------------------
# GAMMA DIVERSITY
#-------------------------------------------------------------------------------

library(iNEXT)
library(dplyr)

ab_df <- read_csv("Data-raw/pool-community-abundance-data.csv")

# pools to remove - App14 has a missing depth, Kam11 is too deep and Kor13 is eutrophicated (most likely). 
remove_pools <- c("App14", "Kam11", "Kor13")

ab_df <- ab_df %>% 
  filter(!Pool %in% remove_pools) %>% 
  dplyr::select(-Pool) %>%
  mutate(across(-Inselberg, round)) # Some values are decimals because of the way abundance was calculated 

row_df <- unique(ab_df$Inselberg)

# Create separate communities for each inselberg 
community <- lapply(row_df, function(x) {
  ab_df %>% 
    dplyr::filter(Inselberg == x) %>%     
    # tidyr::drop_na() %>%                  
    dplyr::select(-Inselberg) 
})

# Calculating the chao gamma div for each inselberg 
gamma <- lapply(seq_along(community), function(i) {
  inselberg_name <- row_df[i]
  loc <- community[[i]]
  loc[loc > 0] <- 1
  loc_NEW <- c(nrow(loc), colSums(loc))
  chao_gamma <- ChaoRichness(loc_NEW, datatype = "incidence_freq", conf = 0.95)
  
  data.frame(
    Inselberg = inselberg_name,
    Estimator = chao_gamma$Estimator[1],  # Gamma richness
    s.e. = chao_gamma$Est_s.e.[1]
  )
})

gamma_df <- bind_rows(gamma)

gamma_df <- gamma_df %>%
  rename(
    Gamma = Estimator,
    GSE = 's.e.')

#-------------------------------------------------------------------------------
# SEPPARATING ACTIVE AND PASSIVE DISPERSERS  

# Load data
AP_df <- read_csv("Data-raw/active-passive-split-taxa.csv")
ab_df <- read_csv("Data-raw/pool-community-abundance-data.csv")

# Define labels and taxa
labels_cols <- c("Inselberg", "Pool")
all_taxa <- setdiff(names(ab_df), labels_cols)

# Classify taxa
active_cols <- intersect(all_taxa, AP_df$Active)
passive_cols <- intersect(all_taxa, AP_df$Passive)

# Check for unclassified taxa
unclassified <- setdiff(all_taxa, c(active_cols, passive_cols))
unclassified

# Subset data
ab_df_Active <- ab_df %>% dplyr::select(all_of(labels_cols), all_of(active_cols))
ab_df_Passive <- ab_df %>% dplyr::select(all_of(labels_cols), all_of(passive_cols))

#-------------------------------------------------------------------------------

# GAMMA DIVERSITY FOR ACTIVE DISPERSERS 

# pools to remove - App14 has a missing depth, Kam11 is too deep and Kor13 is eutrophicated (most likely). 
remove_pools <- c("App14", "Kam11", "Kor13")

ab_df_Active <- ab_df_Active %>% 
  filter(!Pool %in% remove_pools) %>% 
  dplyr::select(-Pool) %>%
  mutate(across(-Inselberg, round)) # Some values are decimals because of the way abundance was calculated 

row_df <- unique(ab_df_Active$Inselberg)

# Create separate communities for each inselberg 
community_Active <- lapply(row_df, function(x) {
  ab_df_Active %>% 
    dplyr::filter(Inselberg == x) %>%     
    dplyr::select(-Inselberg) 
})

# Calculating the chao gamma div for each inselberg 
gamma_Active <- lapply(seq_along(community_Active), function(i) {
  inselberg_name <- row_df[i]
  loc <- community_Active[[i]]
  loc[loc > 0] <- 1
  loc_NEW <- c(nrow(loc), colSums(loc))
  chao_gamma <- ChaoRichness(loc_NEW, datatype = "incidence_freq", conf = 0.95)
  
  data.frame(
    Inselberg = inselberg_name,
    Estimator = chao_gamma$Estimator[1],  # Gamma richness
    s.e. = chao_gamma$Est_s.e.[1]
  )
})

gamma_df_Active <- bind_rows(gamma_Active)

gamma_df_Active <- gamma_df_Active %>%
  rename(
    Gamma = Estimator,
    GSE = 's.e.')

#-------------------------------------------------------------------------------

# GAMMA DIVERSITY FOR PASSIVE DISPERSERS 

# pools to remove - App14 has a missing depth, Kam11 is too deep and Kor13 is eutrophicated (most likely). 
remove_pools <- c("App14", "Kam11", "Kor13")

ab_df_Passive <- ab_df_Passive %>% 
  filter(!Pool %in% remove_pools) %>% 
  dplyr::select(-Pool) %>%
  mutate(across(-Inselberg, round)) # Some values are decimals because of the way abundance was calculated 

row_df <- unique(ab_df_Passive$Inselberg)

# Create separate communities for each inselberg 
community_Passive <- lapply(row_df, function(x) {
  ab_df_Passive %>% 
    dplyr::filter(Inselberg == x) %>%     
    dplyr::select(-Inselberg) 
})

# Calculating the chao gamma div for each inselberg 
gamma_Passive <- lapply(seq_along(community_Passive), function(i) {
  inselberg_name <- row_df[i]
  loc <- community_Passive[[i]]
  loc[loc > 0] <- 1
  loc_NEW <- c(nrow(loc), colSums(loc))
  chao_gamma <- ChaoRichness(loc_NEW, datatype = "incidence_freq", conf = 0.95)
  
  data.frame(
    Inselberg = inselberg_name,
    Estimator = chao_gamma$Estimator[1],  # Gamma richness
    s.e. = chao_gamma$Est_s.e.[1]
  )
})

gamma_df_Passive <- bind_rows(gamma_Passive)

gamma_df_Passive <- gamma_df_Passive %>%
  rename(
    Gamma = Estimator,
    GSE = 's.e.')

#-------------------------------------------------------------------------------
# DEPTH
#-------------------------------------------------------------------------------

Raw_df <- read_csv("Data-raw/raw-data-alpha.csv")

# pools to remove - App14 has a missing depth, Kam11 is too deep and Kor13 is eutrophicated (most likely). 
remove_pools <- c("App14", "Kam11", "Kor13")

depth_df <- Raw_df %>% 
  dplyr::select(c(Inselberg, Pool, Depth)) %>% 
  filter(!Pool %in% remove_pools) 

# Calculate median depth and mean log depth for each inselberg
depth_gamma <- depth_df %>%
  dplyr::select(Inselberg, Depth) %>% 
  mutate(Inselberg = factor(Inselberg, levels = unique(Inselberg))) %>%
  group_by(Inselberg) %>%
  summarise(Depth_median = median(Depth, na.rm = TRUE),
            Meanlogdepth = mean(log(Depth), na.rm = TRUE))

#-------------------------------------------------------------------------------
# PROMINENCE
#-------------------------------------------------------------------------------

Raw_data <- read_csv("Data-raw/raw-data-gamma.csv")

prom_df <- Raw_data %>% dplyr::select(Inselberg, Prominence)

#-------------------------------------------------------------------------------
# BIOCLIM variables (CLIMATE PC1 and PC2) 
#-------------------------------------------------------------------------------

library(sp)
library(raster)
library(geodata)
library(writexl)

Raw_df_1 <- read_csv("Data-raw/raw-data-gamma.csv")

# Loading coordinates of Inselbergs 
coords <- Raw_df_1 %>% 
  dplyr::select(c(Longitude, Latitude)) %>% 
  mutate(
    Longitude = as.numeric(Longitude),
    Latitude = as.numeric(Latitude)
  )

points <- SpatialPoints(coords, proj4string=CRS("+init=epsg:4326"))
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

#Download the bioclim data (it downloads as .tiff files)
bio <- worldclim_global(var="bio", res=2.5,path = "C:/Users/apathakb/Downloads/BEF_NEW", lon=coords$Longitude, lat=coords$Latitude) 

# Load the downloaded raster files 
raster_files <- list.files("C:/Users/apathakb/Downloads/BEF_NEW/climate/wc2.1_2.5m", pattern = "*.tif", full.names = TRUE)

# Read the raster files into a stack
raster_stack <- stack(raster_files)

# Create SpatialPoints object for your coordinates
coordinates(coords) <- ~Longitude+Latitude
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")

# Extract the values for your coordinates
extracted_df <- extract(raster_stack, coords)

# Arranging it in ascending order - bio1 to bio19 
current_colnames <- colnames(extracted_df)
sorted_colnames <- current_colnames[order(as.numeric(gsub(".*_(\\d+)$", "\\1", current_colnames)))]
extracted_df <- extracted_df[, sorted_colnames]

# Assign correct bioclim variable names and rownames (inselbergs)
bioclim <- data.frame(extracted_df)
colnames(bioclim)=c("Annual Tmean","Mean diurnal range","Isothermality","T seasonality","Tmax warmest month","Tmin coldest month","T annual range","Tmean wettest quarter","Tmean driest quarter","Tmean warmest quarter","Tmean coldest quarter","Annual P","P wettest month","P driest month","P seasonality","P wettest quarter","P driest quarter","P warmest quarter","P coldest quarter")
rownames(bioclim) <- c("IVC","SPN","KAM","SWE","FRA","USA","KOR","SWA","NWA","MAL")
print(bioclim)

# write to a .csv file for use in later analyses
readr::write_csv(x = bioclim, file = "Data-raw/bioclim-data.csv")


#-------------------------------------------------------------------------------

# PCA - 
library(vegan)
library(ggfortify)
library(ggplot2)
library(ggrepel)
library(factoextra)

pca1 <- prcomp(bioclim, scale = TRUE, center = TRUE)
summary(pca1)

# use the factorextra package to visualise the PCA

Inselberg <- rownames(bioclim)

p1 <- 
  factoextra::fviz_pca_biplot(pca1, repel = TRUE, 
                              label = c("var","ind"), col.var="black", alpha.var=0.2,labelsize = 3) +
  geom_point(aes(shape = NULL, colour = Inselberg), size = 3) +
  scale_colour_manual(values = wesanderson::wes_palette(name = "Darjeeling1", n = 10, type = "continuous")) +
  labs(title = NULL, x="PC1 (48.9%)", y="PC2 (24.3%)") +
  theme_classic() +
  theme(legend.position = "none")
plot(p1)

# Export figure
ggsave("fig_s1.1.pdf", p1, units = "cm", width = 13, height = 11)

# run a pca to extract the scores using rda from vegan
pca2 <- vegan::rda(bioclim, scale=TRUE)
summary(pca2)

# extract the scores
pca2_scores <- scores(pca2)

# convert to a data.frame
pca2_scores <- as.data.frame(pca2_scores$sites)

# remove the row.names
pca2_scores$Inselberg <- row.names(pca2_scores)

# convert to factor and rearrange
pca2_scores <- 
  pca2_scores %>% 
  dplyr::mutate(Inselberg = factor(Inselberg)) %>% 
  dplyr::select(Inselberg, PC1, PC2)

# remove the row.names
row.names(pca2_scores) <- NULL

# write to a .csv file for use in later analyses
readr::write_csv(x = pca2_scores, file = "Data-raw/bioclim-pca-scores.csv")

#-------------------------------------------------------------------------------
# MAKING NEW RAW DATA (Subsets) 
#-------------------------------------------------------------------------------

# ALPHA SCALE

alpha_all_raw <- alpha_df %>% 
  left_join(BM_df_1, by = c("Inselberg", "Pool")) %>% 
  left_join(gamma_df, by = "Inselberg") %>% 
  left_join(depth_df, by = c("Inselberg", "Pool")) %>% 
  left_join(pca2_scores, by = "Inselberg")

readr::write_csv(alpha_all_raw, "Data-raw/alpha-all-raw.csv")

alpha_active_raw <- alpha_df_Active %>% 
  left_join(BM_df_Active_1, by = c("Inselberg", "Pool")) %>% 
  left_join(gamma_df_Active, by = "Inselberg") %>% 
  left_join(depth_df, by = c("Inselberg", "Pool")) %>% 
  left_join(pca2_scores, by = "Inselberg")

readr::write_csv(alpha_active_raw, "Data-raw/alpha-active-raw.csv")

alpha_passive_raw <- alpha_df_Passive %>% 
  left_join(BM_df_Passive_1, by = c("Inselberg", "Pool")) %>% 
  left_join(gamma_df_Passive, by = "Inselberg") %>% 
  left_join(depth_df, by = c("Inselberg", "Pool")) %>% 
  left_join(pca2_scores, by = "Inselberg")

readr::write_csv(alpha_passive_raw, "Data-raw/alpha-passive-raw.csv")

#-------------------------------------------------------------------------------

# GAMMA SCALE 

gamma_all_raw <- gamma_df %>% 
  left_join(Biomass_gamma, by = "Inselberg") %>% 
  left_join(depth_gamma, by = "Inselberg") %>% 
  left_join(pca2_scores, by = "Inselberg") %>% 
  left_join(prom_df, by = "Inselberg")

readr::write_csv(gamma_all_raw, "Data-raw/gamma-all-raw.csv")

gamma_active_raw <- gamma_df_Active %>% 
  left_join(Biomass_gamma_Active, by = "Inselberg") %>% 
  left_join(depth_gamma, by = "Inselberg") %>% 
  left_join(pca2_scores, by = "Inselberg") %>% 
  left_join(prom_df, by = "Inselberg")

readr::write_csv(gamma_active_raw, "Data-raw/gamma-active-raw.csv")

gamma_passive_raw <- gamma_df_Passive %>% 
  left_join(Biomass_gamma_Passive, by = "Inselberg") %>% 
  left_join(depth_gamma, by = "Inselberg") %>% 
  left_join(pca2_scores, by = "Inselberg") %>% 
  left_join(prom_df, by = "Inselberg")

readr::write_csv(gamma_passive_raw, "Data-raw/gamma-passive-raw.csv")

#-------------------------------------------------------------------------------
# END
#-------------------------------------------------------------------------------


