
# function to prep the environment
# taxa: "all", "active", "passive" (or "")
# load relevant libraries
library(lme4)
library(dagitty)
library(splines)
library(DHARMa)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(piecewiseSEM)

# get a list of the helper functions
function_list <- list.files(here::here("Code/"))
function_list <- function_list[grepl(pattern = "helper-", function_list)]

# load the helper functions
for (i in function_list){
  print(i)
  source(here::here(file.path("Code", i)))
}

# get a list of data files
files <- list.files(here::here("Data/"))

# get the relevant file based on the taxa
file_sel <- files[grepl(pattern = paste0(taxa, "-data-clean.rds"), files)]

# load the alpha-scale data
if (taxa != "") {
  dat <- readRDS(here::here(file.path("Data", file_sel)))
}


# set-up frequently used axis labels
alpha_div <- expression("ln("~alpha~"-diversity )")
gamma_div <- expression("ln("~gamma~"-diversity )")
ln_biomass <- "ln( biomass ) (mg)"
ln_depth <- "ln (depth) (cm)"
