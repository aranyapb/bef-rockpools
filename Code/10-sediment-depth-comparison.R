#'
#' @title Supplementary figure examining pool sediment depth
#' 
#' @description Examine the relationship between pool size and sediment
#' depth to see if this can help us explain the results
#'

# load the relevant libraries
library(dplyr)
library(ggplot2)
library(readr)
library(wesanderson)

# load plotting theme
source("Code/helper-plotting-theme.R")

# load the sediment data
sed_dat <- readr::read_csv2("Data/pool-sediment-data.csv")
head(sed_dat)

# convert the lndepth and lnsed variables to numeric variables
sed_dat <- 
  sed_dat |>
  dplyr::mutate(Lndepth = as.numeric(Lndepth),
                LnSedPerVRatio = as.numeric(LnSedPerVRatio))

# rename the Location column
sed_dat <- dplyr::rename(sed_dat, Inselberg = Location)

# check the inselberg names
unique(sed_dat$Inselberg)
sed_dat <- 
  sed_dat |>
  dplyr::mutate(Inselberg = ifelse(Inselberg == "SOZ", "SWA",
                                   ifelse(Inselberg == "NOZ", "NWA", Inselberg)))

# convert inselberg to a factor
sed_dat <- dplyr::mutate(sed_dat, Inselberg = factor(Inselberg))

# sort out the coloours because not all inselbergs are present here
inselbergs <- c("IVC","SPN","KAM","SWE","FRA","USA","KOR","SWA","NWA","MAL")
inselbergs <- sort(inselbergs)
cols <- wesanderson::wes_palette(name = "Darjeeling1", n = 10, type = "continuous")

# which inselbergs are present
cols <- cols[which(inselbergs %in% unique(sed_dat$Inselberg))]

# perform a correlation test on these data
sed_cor <- cor.test(sed_dat$Lndepth, sed_dat$LnSedPerVRatio, method="pearson") # r = -0.73; P<0.001
sed_est <- round(sed_cor$estimate, 2)
  
# plot fig. s10
p1 <- 
  ggplot(data = sed_dat, mapping = aes(x = Lndepth, y = LnSedPerVRatio, colour = Inselberg)) +
  geom_point(shape = 1, stroke = 0.5, alpha = 0.8, size = 1.75) +
  labs(x = "Ln( depth ) (cm)", y = "Ln( sediment/volume ) (% cover/m³)")+
  scale_colour_manual(values = cols) +
  guides(size = "none", 
         colour = guide_legend(override.aes = list(size = 1.5, 
                                                   alpha = 1, shape = 1, stroke = 1),
                                                   ncol = 1)) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("r = ", sed_est, " *")) +
  theme_meta() +
  theme(legend.position = "right")+
  theme(legend.key = element_blank())

plot(p1)

# export the figure
ggsave("Figures/fig_s10.pdf", p1, 
       units = "cm", width = 12, height = 9.5)

### END
