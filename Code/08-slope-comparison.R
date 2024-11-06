#'
#' @title Compare BEF slopes between active and passive dispersers
#' 
#' @description Estimate the BEF slope for active and passive dispersers
#' and compare them to produce fig. s9
#'

# load relevant libraries
library(dplyr)
library(ggplot2)
library(wesanderson)

# load plotting theme
source("Code/helper-plotting-theme.R")

# get a list of data files
files <- list.files("Data/")

#'
#' @title Slope variation analysis
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

# make the slope plot
p1 <- 
  ggplot(data = slope_df,
       mapping = aes(x = passive_Est, y = active_Est, colour = Inselberg)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_vline(xintercept = 0, size = 0.25) +
  geom_point() +
  geom_errorbarh(mapping = aes(xmin = passive_Est - passive_SE,
                               xmax = passive_Est + passive_SE), height = 0, 
                 show.legend = FALSE) +
  geom_errorbar(mapping = aes(ymin = active_Est - active_SE,
                               ymax = active_Est + active_SE), width = 0,
                show.legend = FALSE) +
  scale_colour_manual(values = wes_palette(name = "Darjeeling1", n = 10, type = "continuous")) +
  ylab("BEF slope (±SE) active dispersers") +
  xlab("BEF slope (±SE) passive dispersers") +
  theme_meta()+
  theme(legend.key = element_blank())

plot(p1)
  
# export the figure
ggsave("Figures/fig_s9.pdf", p1, 
       units = "cm", width = 12, height = 9.5)

# Get slopes and SE 
std <- function(x) sd(x)/sqrt(length(x))
mean(slope_df$active_Est)
std(slope_df$active_Est)
mean(slope_df$passive_Est)
std(slope_df$passive_Est)

# test for differences in slope using a wilcox-test
wilcox.test(slope_df$active_Est, slope_df$passive_Est, paired= TRUE)

# use a one-sample t-test
t.test(x = (slope_df$active_Est - slope_df$passive_Est))

### END


