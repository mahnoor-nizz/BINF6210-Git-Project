###############################################################
# Biodiversity Patterns in Cardinalidae (BOLD Data)
# Analysis of Barcode Index Number (BIN) diversity patterns across geographic and elevational gradients
#
# Mahnoor Nizamani
# 2025-10-09
# Data Source: BOLD Systems (www.boldsystems.org)
###############################################################


#----- Packages Used ------

library(dplyr)      # Wickham et al. (2023). dplyr: A Grammar of Data Manipulation
library(ggplot2)    # Wickham (2016). ggplot2: Elegant Graphics for Data Analysis
library(vegan)      # Oksanen et al. (2022). vegan: Community Ecology Package
library(tidyr)      # Wickham et al. (2023). tidyr: Tidy Messy Data
library(scales)     # Wickham & Seidel (2022). scales: Scale Functions for Visualization
library(readr)      # Wickham et al. (2023). readr: Read Rectangular Text Data

#----- Load Data ------

# Data acquired from BOLD Systems API
# Code I used to acquire the data (data download performed on Oct 08, 2025):
# Cardinalidae.df <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Cardinalidae&format=tsv")

Cardinalidae.df <- read_tsv(file = "../data/Cardinalidae bold data.txt")

#----- Comments -----

# I downloaded the Data onto my laptop from the BOLD website at first but I did not realize that when you click to view the data frame it only shows you the first 50 column and that the rest are on the next page. I then used the code above to acquire the data and saw that it had the longitude and latitude data but no country data(on the first page ;-;), so I spent way too much time trying to merge the two data frames together. They had different column names and different number of rows. I had to find a common key (column that uniquely identified specimens & matched both data sets) to complete the merge. It was a pain but I learned a lot from it. After successfully merging the two data frames, I noticed that the merged data frame had DOUBLE the variables, that's when I realized what a silly mistake I had made haha.

# I wasted a lot of time on that so when I had issues with statistical tests that I wanted to run later on, I had to just scrap them as I ran out of time before the due date.

#----- Initial Steps -----

# Check initial dataset dimensions
cat("Rows:", nrow(Cardinalidae.df), "\n")
# Rows: 597
cat("Columns:", ncol(Cardinalidae.df), "\n\n")
# Columns: 80

# Summary of key variables
cat("Unique BINs:", n_distinct(Cardinalidae.df$bin_uri, na.rm = TRUE), "\n")
# Unique BINs: 58
cat("Unique species:", n_distinct(Cardinalidae.df$species_taxID, na.rm = TRUE), "\n")
# Unique species: 39
cat("Unique countries:", n_distinct(Cardinalidae.df$country, na.rm = TRUE), "\n")
# Unique countries: 16
cat("Records with elevation data:", sum(!is.na(Cardinalidae.df$elev)), "\n\n")
# Records with elevation data: 170

# Noticed there were a lot of NAs in some columns so I began by removing columns that were mostly NAS

na_percentage <- colMeans(is.na(Cardinalidae.df))
columns_to_keep <- na_percentage <= 0.90
Cardinalidae.df <- Cardinalidae.df[, columns_to_keep]

cat("Columns after NA filter:", ncol(Cardinalidae.df), "\n")
# Columns after NA filter: 49

# Removed some columns that I didn't need to use later on (not necessary but I wanted to be able to look through the data frame more easily)
Cardinalidae.df <- Cardinalidae.df[, -c(4, 20, 21, 22, 24, 29, 36, 44)]

cat("Columns after manual removal:", ncol(Cardinalidae.df), "\n")
# Columns after manual removal: 41


################################################################
# Question 1: Dpes Sampling Completeness and Bin Richness Differs Among Countries?
#
# Rationale: Different countries may have varying sampling effort, which affects
# our ability to estimate true BIN diversity/richness. We use Chao1 estimator to assess
# sampling completeness.
#
# Method: Species accumulation curves and Chao1 richness estimator. (Chao, 1984)
################################################################

#----- Created presence/absence matrix for BINs by country -----

# Each row = site (sampling event), columns = BINs, NAs filtered.
pa_matrix <- Cardinalidae.df %>%
  filter(!is.na(bin_uri), !is.na(country)) %>%
  mutate(site = row_number(), presence = 1) %>%
  pivot_wider(
    id_cols = c(site, country),
    names_from = bin_uri,
    values_from = presence,
    values_fill = 0
  )


cat("Sites:", nrow(pa_matrix), "\n")
# Sites: 467
cat("BINs:", ncol(pa_matrix) - 2, "\n")  # Subtract site and country columns.
# BINs: 55


# Extracted numeric matrix (specpool would not run when encountering any non-numeric value) and country vector (component needed for analysis).
pa_numeric <- pa_matrix %>% select(-site, -country)
country_vector <- pa_matrix$country


# Verify no NAs in matrix
cat("NAs in PA matrix:", sum(is.na(pa_numeric)), "\n\n")
# NAs in PA matrix: 0


#----- Estimated richness per country using Chao1 -----


richness_per_country <- specpool(pa_numeric, country_vector)
richness_per_country$country <- rownames(richness_per_country)

# Richness estimates
print(richness_per_country %>%
  select(country, Species, chao, chao.se) %>%
  arrange(desc(Species)))

#----- Calculated completeness and prepared for plotting -----

# Completeness = Observed / Estimated richness
richness_per_country <- richness_per_country %>%
  mutate(
    completeness = Species / chao,
    completeness_scaled = completeness * max(Species),
    country = factor(country, levels = country[order(-Species)])
  )

# Completeness statistics data summary
cat("Mean completeness:", round(mean(richness_per_country$completeness, na.rm = TRUE), 3), "\n")
# Mean completeness: 0.856
cat("Range:", round(min(richness_per_country$completeness, na.rm = TRUE), 3), "-",
  round(max(richness_per_country$completeness, na.rm = TRUE), 3), "\n\n")
# Range: 0.591 - 1

#----- Visualize observed vs estimated richness with completeness -----

ggplot(richness_per_country, aes(x = country)) +
  geom_bar(aes(y = Species), stat = "identity", fill = "skyblue") +
  geom_point(aes(y = chao), color = "red", size = 3) +
  geom_errorbar(
    data = filter(richness_per_country, chao.se > 0),
    aes(ymin = chao - chao.se, ymax = chao + chao.se),
    width = 0.2, color = "red"
  ) +
  geom_point(aes(y = completeness_scaled), color = "darkgreen", size = 2) +
  geom_line(aes(y = completeness_scaled, group = 1), color = "darkgreen", linewidth = 1) +
  labs(
    x = "Country",
    y = "BIN richness",
    title = "Sampling Completeness and BIN Richness per Country",
    subtitle = "Blue = Observed | Red = Chao1 estimate | Green = Completeness (scaled)",
    caption = "Error bars represent Chao1 standard error"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  ) 
################################################################
# Question 2: Does BIN Composition Differs Between Regions?
#
# Rationale: Geographic isolation and environmental differences may lead to
# distinct BIN assemblies across countries.
#
# Method: Non-metric Multidimensional Scaling (NMDS) with Jaccard dissimilarity
# (Jaccard, 1912; Kruskal, 1964)
################################################################

# NMDS Ordination with Visualization

# rather than making a new data frame with countries as the rows and associated bin_uri as column values. easier to just use the presence absence matrix as before. Convert to proper matrix (not table)
# Rows = countries, columns = BINs, values = presence (1) or absence (0)
pa_matrix_country <- as.matrix(table(Cardinalidae.df$bin_uri, Cardinalidae.df$country))
pa_matrix_country[pa_matrix_country > 0] <- 1

# Country-BIN matrix data
cat("Countries:", nrow(t(pa_matrix_country)), "\n")
# Countries: 16
cat("BINs:", ncol(t(pa_matrix_country)), "\n")
# BINs: 58

# Transpose so countries are rows (required for vegan functions)
pa_matrix_t <- t(pa_matrix_country)
class(pa_matrix_t) <- "matrix"  # Force matrix class

#----- Run NMDS Ordination -----

# Jaccard distance appropriate for presence/absence data
# k=2 for 2D visualization
# trymax=100 to find stable solution

set.seed(123)  # For reproducibility
# Now NMDS works
nmds <- metaMDS(pa_matrix_t, distance = "jaccard", k = 2, trymax = 100)

# Check stress
nmds$stress  # Should be < 0.1
# [1] 0.02944058

# Extract scores
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$Country <- rownames(nmds_scores)

#----- Plot NMDS Scores -----

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, label = Country)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(hjust = -0.2, vjust = 0.4, size = 3, fontface = "bold") +
  labs(
    title = "Geographic Variation in BIN Composition",
    subtitle = paste0("NMDS ordination (Jaccard dissimilarity) | Stress = ",
      round(nmds$stress, 3)),
    x = "NMDS Axis 1",
    y = "NMDS Axis 2",
    caption = "Countries closer together have more similar BIN compositions"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )

####Start of Improvement 1

nmds_scores <- nmds_scores %>%
  left_join(richness_per_country %>%
              select(country, Species, completeness),
            by = c("Country" = "country"))

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, label = Country, color = completeness)) +
  geom_point(size = 2) +
  geom_text(hjust = -0.2, vjust = 0.4, size = 3) +
  scale_color_gradient(
    low = "red",
    high = "blue",
    name = "Sampling\ncompleteness"
  ) +
  labs(
    title = "Geographic Variation in BIN Composition and Sampling completeness",
    subtitle = paste0("NMDS ordination (Jaccard dissimilarity index) | Stress = ",
                      round(nmds$stress, 2)),
    x = "Major differences in BIN composition",
    y = "Secondary differences in BIN composition",
    caption = "Closer points = Higher BIN  composition similarity"
  ) +
  coord_cartesian(xlim = c(min(nmds_scores$NMDS1), max(nmds_scores$NMDS1) + 0.1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5),
    legend.position = "right"
  )

####End of Improvement 1


################################################################
# Question 3: Does Elevation Influence BIN Diversity?
#
# Rationale: Elevation gradients represent major environmental transitions that
# may structure species distributions (Körner, 2007).
#
# Method: Binning elevation data and calculating BIN richness per bin
################################################################


#----- Filter and Prepare Elevation Data -----

Elev.Bin.df <- Cardinalidae.df %>%
  filter(!is.na(elev), !is.na(bin_uri)) %>%
  select(bin_uri, elev)


# Elevation data
cat("Records with elevation:", nrow(Elev.Bin.df), "\n")
# Records with elevation: 162
cat("Elevation range:", round(min(Elev.Bin.df$elev)), "-",
  round(max(Elev.Bin.df$elev)), "m\n")
# Elevation range: 8 - 3650 m

#----- Create Elevation Bins -----
# Using 250m bins as a biologically meaningful scale

Elev.Bin.df <- Elev.Bin.df %>%
  mutate(elev_bin = cut(elev,
    breaks = seq(0, ceiling(max(elev) / 250) * 250, by = 250),
    include.lowest = TRUE,
    dig.lab = 10))  # Prevents scientific notation

#----- Calculate BIN Richness per Elevation Bin -----
elev_counts <- Elev.Bin.df %>%
  group_by(elev_bin) %>%
  summarise(n_BINS = n_distinct(bin_uri), .groups = "drop") %>%
  mutate(elev_label = gsub(",", "-", elev_bin),  # Replace comma with dash
    elev_label = factor(elev_label, levels = gsub(",", "-", levels(elev_bin))))  # Keep order


#----- Visualize Elevational Pattern -----
ggplot(elev_counts, aes(x = elev_label, y = n_BINS, group = 1)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(size = 2, color = "darkgreen") +
  labs(x = "Elevation (m)",
    y = "BIN richness",
    title = "Elevational Diversity Gradient in Cardinalidae",
    subtitle = "BIN richness across 250m elevation bins",
    caption = "Line shows trend; bars show observed richness") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10))


######Start of Improvement 2

ggplot(elev_counts, aes(x = elev_label, y = n_BINS)) +
  geom_bar(stat = "identity", fill = "orange") +       
  geom_smooth(aes(group = 1), method = "loess", se = TRUE, 
              color = "blue", linewidth = 0.9) +
  labs(
    x = "Elevation (meters)",
    y = "BIN richness",
    title = "BIN Diversity along an Elevational Gradient in Cardinalidae",
    subtitle = "Trend of BIN richness across 250m elevation bins",
    caption = "Blue line = LOESS trend. Bars = Observed richness"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, hjust=0.5),
    plot.subtitle = element_text(size = 10, hjust=0.5)
  )

######End of Improvement 2


######Start of Improvement 3
install.packages('rnaturalearth')
install.packages("rnaturalearthdata")
install.packages("sf")
install.packages("viridis")

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)


world <- ne_countries(scale = "medium", returnclass = "sf")

world_richness <- world %>%
  left_join(richness_per_country, by = c("name" = "country"))

ggplot(world_richness) +
  geom_sf(aes(fill = Species), color = "black") +   
  scale_fill_viridis_c(option = "viridis", na.value = "gray90") +
  labs(
    title = "Worldwide distribution of BIN richness",
    fill = "BIN richness"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust=0.5)
  )

######End of Improvement 3


################################################################
# END OF SCRIPT
################################################################
