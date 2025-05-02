#'
#' @title Slope variation analysis and Compare BEF slopes between active and passive dispersers
#' 
#' @description Fit linear models to show the interaction between inselberg
#' and alpha diversity on biomass.
#'

# load relevant libraries
library(dplyr)
library(ggplot2)
library(wesanderson)

# load plotting theme
source("Code/helper-plotting-theme.R")

# get a list of data files
files <- list.files("Data/")

# which groups
taxa <- c("active", "passive")

# apply over the two taxa
slope_list <- 
  
  lapply(taxa, function(x) {
    
  # load the alpha scale active dispersers
  alph_dat <- readRDS(paste0("Data/", files[ grepl(pattern = paste0(x, "-taxa-alpha.rds"), files)] ))
  head(alph_dat)
  
  # fit a linear model and extract the slope of the standardised coefficients
  alph_dat <- 
    alph_dat |>
    dplyr::mutate(logalpha_std = scale(logalpha)[,1],
                  logbiomass_std = scale(logbiomass)[,1],
                  logdepth_std = scale(logdepth)[,1])
  
  # fit a linear model and extract the slope of the standardised coefficients
  lm_alph <- lm(logbiomass_std ~ logdepth_std + Inselberg + logalpha_std:Inselberg, data = alph_dat)
  plot(alph_dat$logbiomass_std, predict(lm_alph))
  abline(0, 1)
  
  # extract the passive disperser BEF slopes
  alph_b <- data.frame( summary(lm_alph)$coefficients )
  
  # make the row.names into a column
  alph_b <- dplyr::bind_cols(data.frame(taxa = x, term = row.names(alph_b)), alph_b)
  
  # remove row.names
  row.names(alph_b) <- NULL
  
  # rename the columns
  names(alph_b) <- c("taxa", "term", paste0(x, "_Est") , paste0(x, "_SE"), "t-val", "P-val")
  
  # subset the relevant columns
  alph_b <- dplyr::select(alph_b, term, paste0(x, "_Est"), paste0(x, "_SE"))
  
  # filter the correct terms
  alph_b <- dplyr::filter(alph_b, grepl(pattern = ":logalpha_std", term) )
  
  alph_b
  
} )

# combine into a data.frame
slope_df <- dplyr::full_join(slope_list[[1]], slope_list[[2]], by = "term")

# make an inselberg column
slope_df$Inselberg <- substr(slope_df$term, start = 10, stop = 12 )
slope_df$Inselberg <- factor(slope_df$Inselberg)

# Get slopes and SE 
std <- function(x) sd(x)/sqrt(length(x))

# for acive dispersers
mean(slope_df$active_Est)
std(slope_df$active_Est)

# for passive dispersers
mean(slope_df$passive_Est)
std(slope_df$passive_Est)

# test for differences in slope using a wilcox-test
wilcox.test(slope_df$active_Est, slope_df$passive_Est, paired= TRUE)

# use a one-sample t-test
t.test(x = (slope_df$active_Est - slope_df$passive_Est))

### END


