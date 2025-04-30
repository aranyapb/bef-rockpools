
# function to prep the environment
# taxa: "all", "active", "passive"
# load relevant libraries
library(lme4)
library(dagitty)
library(splines)
library(DHARMa)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(piecewiseSEM)

# load relevant functions
source("Code/helper-plotting-theme.R")
source("Code/helper-plot-functions.R")
source("Code/helper-ci-lmm.R")
source("Code/helper-compare-lmms.R")
source("Code/helper-dfbetas-lmm.R")
source("Code/helper-get-dags.R")
source("Code/helper-fisher-test.R")

# get a list of data files
files <- list.files("Data/")

# load the alpha-scale data
dat <- readRDS(paste0("Data/", files[ grepl(pattern = paste0(taxa, "-taxa-alpha.rds"), files)] ))

# set-up frequently used axis labels
alpha_div <- expression("ln("~alpha~"-diversity )")
gamma_div <- expression("ln("~gamma~"-diversity )")
ln_biomass <- "ln( biomass ) (mg)"
ln_depth <- "ln (depth) (cm)"
