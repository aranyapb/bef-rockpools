#'
#' @title Clean the raw data
#' 
#' @description Load the raw data, clean the data and output cleaned versions
#' for analysis
#'

# load relevant libraries
library(dplyr)
library(readr)

# get a list of the data files
files <- list.files("Data/")

# make a vector of file names for the different taxa
tax_names <- c("all", "active", "passive")

# gamma-scale

# get the gamma datasets
gam_files <- files[grepl(pattern = "gamma-", x = files)]

# order the gam files by the gam names
gam_files <- gam_files[sapply(tax_names, function(x) which(grepl(pattern = x, x = gam_files)))]

# read in each file and run the cleaning functions
for(i in 1:length(gam_files)) {
  
  # gamma scale
  dat_gam <- read_delim(paste0("Data/", gam_files[i]), delim = ";")
  
  # add ln-transformed variables
  dat_gam <- 
    dat_gam |>
    rename(Inselberg = `Inselberg `) |>
    mutate(loggamma = log(Gamma),
           logprom = log(Prominence),
           Inselberg = factor(Inselberg))
  
  # select the relevant columns
  dat_gam <-
    dat_gam |>
    select(Inselberg, PC1, PC2, Gamma, GSE, loggamma, Prominence, logprom, Slope, SSE) |>
    rename(Gamma_SE = GSE, Slope_SE = SSE)
  
  # save as a .rds object
  saveRDS(object = dat_gam, 
          file = paste0("Data/0", i, "-", tax_names[i], "-taxa-gamma.rds" ))
  
}

# alpha-scale

# get the gamma datasets
alph_files <- files[grepl(pattern = "alpha-", x = files)]

# order the gam files by the gam names
alph_files <- alph_files[sapply(tax_names, function(x) which(grepl(pattern = x, x = alph_files)))]

# set-up the datapoints that need to be removed
alpha_rem <- c(114, 102, 112)

# read in each file and run the cleaning functions
for(i in 1:length(alph_files)) {
  
  # alpha scale 
  dat_alph <- read_delim(paste0("Data/", alph_files[i]), delim = ";")
  
  # remove a datapoint without depth information
  dat_alph <- dat_alph[-alpha_rem[i], ]
  
  # add ln-transformed variables
  dat_alph <- 
    dat_alph |>
    mutate(logalpha = log(Alpha),
           loggamma = log(Gamma),
           logbiomass = log(Biomass),
           logdepth = log(Depth),
           Inselberg = factor(Inselberg))
  names(dat_alph)  
  
  # reorder and rename the columns
  dat_alph <- 
    dat_alph |>
    select(Inselberg, Pool, Depth, logdepth, PC1, PC2, 
           Gamma, GSE, loggamma, Alpha, ASE, logalpha,
           Biomass, logbiomass) |>
    rename(Gamma_SE = GSE, Alpha_SE = ASE)
  
  # save as a .rds object
  saveRDS(object = dat_alph, 
          file = paste0("Data/0", i, "-", tax_names[i], "-taxa-alpha.rds" ))
  
}

### END
