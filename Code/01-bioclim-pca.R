#'
#' @title Dimension reduction of climate variables
#' 
#' @description Use a principal components analysis to reduce the dimensions
#' in the BIOCLIM dataset
#'

# load relevant libraries
library(vegan)
library(ggfortify)
library(ggplot2)
library(ggrepel)
library(factoextra)

# load the bioclim data
bioclim <- readr::read_csv("Data/bioclim-data.csv")

# convert to data.frame because some pca packages need row.names
bioclim <- as.data.frame(bioclim)

# remove the location columns
bioclim <- bioclim[ ,-1] 

# set the rownames
rownames(bioclim) <- c("IVC","SPN","KAM","SWE","FRA","USA","KOR","SWA","NWA","MAL")
Inselberg <- rownames(bioclim)

# run a pca on these climatic variable data
pca1 <- prcomp(bioclim, scale = TRUE, center = TRUE)
summary(pca1)

# use the factorextra package to visualise the PCA
p1 <- 
  factoextra::fviz_pca_biplot(pca1, repel = TRUE, 
                            label = c("var","ind"), col.var="black", alpha.var=0.2,labelsize = 3) +
  geom_point(aes(shape = NULL, colour = Inselberg), size = 3) +
  scale_colour_manual(values = wesanderson::wes_palette(name = "Darjeeling1", n = 10, type = "continuous")) +
  labs(title = NULL, x="PC1 (48.9%)", y="PC2 (24.3%)") +
  theme_classic() +
  theme(legend.position = "none")
plot(p1)

# export the figure as fig. s1.1
ggsave("Figures/fig_s1.1.pdf", p1, 
       units = "cm", width = 13, height = 11)

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
  pca2_scores |>
  dplyr::mutate(Inselberg = factor(Inselberg)) |>
  dplyr::select(Inselberg, PC1, PC2)

# remove the row.names
row.names(pca2_scores) <- NULL

# write to a .csv file for use in later analyses
readr::write_csv(x = pca2_scores, file = "Data/bioclim-pca-scores.csv")

### END
