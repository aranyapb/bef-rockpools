#'
#' @title Make the species accumulation curves
#' 
#' @description Use the community abundance data to examine the 
#' species accumulation curves for the appendix
#'

# load relevant libraries
library(dplyr)
library(readr)
library(vegan)
library(ggplot2)
library(ggpubr)

# load the relevant data
sac_dat <- readr::read_csv2("Data/pool-community-data.csv")

# remove the columbia inselberg
unique(sac_dat$Inselberg)

sac_dat <- 
  sac_dat |>
  dplyr::filter(Inselberg != "COL")

# split the data by location
sac_list <- split(sac_dat, sac_dat$Inselberg)

# function to plot the species accumulation curve
SACplot <- function(data, location) {
  
  # get the community without the identifiers
  sp <- data.frame(sapply(data[,-c(1:2)], as.numeric))
  
  # round abundances
  sp <- round(sp, digits = 0)
  
  # retain columns with at least one taxa
  sp_pres <- sp[, colSums(sp) > 0]
  
  # perform the rarefaction with 10000 permutations
  sac_curve <- specaccum(sp_pres, method = "rarefaction", permutations = 10000)
  
  # get data.frame for the accumulation curve
  sac_curve_df <- data.frame(cbind(sac_curve[["sites"]],sac_curve[["richness"]],sac_curve[["sd"]])) 
  
  # rename the columns
  colnames(sac_curve_df) <- c("sites","richness","sd")
  
  # make sd ribbon intervals
  sac_curve_df[["low_sd"]] <- sac_curve_df[["richness"]]-sac_curve_df[["sd"]]
  sac_curve_df[["high_sd"]] <- sac_curve_df[["richness"]]+sac_curve_df[["sd"]]
  
  p1 <- 
    ggplot(data = sac_curve_df)+
    geom_line(mapping = aes(x = sites, y = richness), size = 0.65) +
    geom_ribbon(mapping = aes(x = sites, ymin = low_sd, ymax = high_sd), alpha = 0.1) +
    labs(x = "", y = "") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1,
             label = location, parse = TRUE) +
    theme_classic()
  
  p1
  
}

# for each inselberg, make an accumlation curve
figs2.1 <- vector("list", length = length(sac_list))
for(i in 1:length(sac_list)) {
  
  # get an SAC plot for each inselberg
  figs2.1[[i]] <- SACplot(data = sac_list[[i]], location = sac_list[[i]]$Inselberg[1])
  
}

# check the figures
plot(figs2.1[[1]])

# appendix: fig. s2.1
figs2.1 <- 
  ggarrange(plotlist = figs2.1, ncol = 4, nrow = 3,
            font.label = list(size = 11, face = "plain"),
            common.legend = TRUE, legend = "bottom")

figs2.1 <- annotate_figure(figs2.1,
                           left = text_grob("Rarefied richness", rot = 90, vjust = 1, size = 14),
                           bottom = text_grob("Number of Rock pools", vjust = -1, size = 14))

plot(figs2.1)

# export the figure
ggsave("Figures/fig_s2.1.pdf", figs2.1, units = "cm", width = 20, height = 13)

### END
