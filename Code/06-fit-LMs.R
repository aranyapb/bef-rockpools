#'
#' @title Slope variation analysis
#' 
#' @description Fit linear models to show the interaction between inselberg
#' and alpha diversity on biomass.
#'

# load relevant libraries
library(lme4)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggbeeswarm)
library(car)

# load plotting theme
source("Code/helper-plotting-theme.R")
source("Code/helper-select-taxa.R")

# get a list of data files
files <- list.files("Data/")

# choose the set of taxa: all, active, passive
taxa <- select_taxa() 

# load the gamma-scale data
gam_dat <- readRDS(paste0("Data/", files[ grepl(pattern = paste0(taxa, "-taxa-gamma.rds"), files)] ))
print(gam_dat)

# check for missing values
any(is.na(gam_dat))

# check the structure of the data
str(gam_dat)

# load the alpha-scale data
alph_dat <- readRDS(paste0("Data/", files[ grepl(pattern = paste0(taxa, "-taxa-alpha.rds"), files)] ))
head(alph_dat)

# check for missing values
any(is.na(alph_dat))

# check the structure of the data
str(alph_dat)

# set-up frequently used axis labels
alpha_div <- expression("ln("~alpha~"-diversity )")
gamma_div <- expression("ln("~gamma~"-diversity )")
ln_biomass <- "ln( biomass ) (mg)"
ln_depth <- "ln (depth) (cm)"

# check summary statistics of alph_dat
summary(alph_dat)

# local-scale models

# fit local models to choose the appropriate model structure

# no effect of individual inselberg
local1 <- lm(logbiomass ~ logdepth + logalpha:Inselberg, data = alph_dat)

# run an analysis of variance
Anova(local1, type="3") # ALL SIGNIFICANT
summary(local1) # R2 = 17.7%

# random intercept for inselberg
local2 <- lmer(logbiomass ~ logdepth + logalpha:Inselberg + (1|Inselberg), data = alph_dat)
Anova(local2, type="3" ) # ALL SIGNIFICANT (depth = 0.04)
summary(local2)

# fixed factor for inselberg
local3 <- lm(logbiomass ~ logdepth + Inselberg + logalpha:Inselberg, data = alph_dat)
Anova(local3, type="3") # DEPTH NOT SIGNIFICANT
summary(local3) # R2 = 25.6%

# fit two different null models: fixed and random effect for inselberg
null1 <- lm(logbiomass ~ Inselberg, data = alph_dat)
null2 <- lmer(logbiomass ~ (1|Inselberg), data = alph_dat) 

# compare the models using AIC: local3 has the lowest AIC
AIC(local1, local2, local3)

# compare the these models to the two null models: local3 has lower AIC than nulls
AIC(local1, local2, local3, null1, null2)

# check assumptions of local3

# can we reject the null hypothesis of normality?
shapiro.test( resid(local3) ) # OK

# can we reject the null hypothesis of equal variance
ncvTest(local3) # OK

# are we justified in modelling the interaction between inselberg and logalpha?
local4 <- lm(logbiomass ~ logdepth + logalpha + Inselberg, data = alph_dat) 

# comppare this model without logalpha interacting with Inselberg to the interaction model
AIC(local3, local4)


# main text: fig. 4

# conditional effect plot holding depth constant
alph_pred <- 
  alph_dat |>
  dplyr::mutate(logdepth = mean(alph_dat$logdepth))

# # predict new logbiomass values but with constant logdepth, se.fit=T for SE
x <- dplyr::as_tibble( predict(local3, alph_pred, interval = "confidence") )
alph_pred$logbiomass <- x$fit
alph_pred$CI_low <- x$lwr
alph_pred$CI_upp <- x$upr

# get correct labels
if(taxa == "all") {
  title <- "All taxa"
  fig <- "a"
  ylab <- ln_biomass
  xlab <- " "
} else if(taxa == "active") {
  title <- "Active dispersers"
  fig <- "b"
  ylab <- NA
  xlab <- alpha_div
} else if(taxa == "passive") {
  title <- "Passive dispersers"
  fig <- "c"
  ylab <- NA
  xlab <- " "
}

p1 <- 
  ggplot() +
  geom_quasirandom(data = alph_dat, 
                   aes(x = logalpha, y = logbiomass, colour = Inselberg, size = Alpha_SE), 
                   width = 0.025, shape = 1, stroke = 0.5, alpha = 0.5) +
  geom_line(data = alph_pred, aes(x = logalpha, y = logbiomass, colour = Inselberg),
            show.legend = FALSE) +
  geom_ribbon(data = alph_pred, 
              aes(x = logalpha, ymax = CI_upp, ymin = CI_low, fill = Inselberg), 
              alpha = 0.1, show.legend = FALSE) +
  ggtitle(title) +
  ylab(if(!is.na(ylab)){ylab}else{NULL}) +
  xlab(xlab) +
  scale_colour_manual(values = wes_palette(name = "Darjeeling1", n = 10, type = "continuous")) +
  scale_fill_manual(values = wes_palette(name = "Darjeeling1", n = 10, type = "continuous")) +
  guides(size = "none", colour = guide_legend(override.aes = list(size = 2.5, alpha = 1, shape = 1, stroke = 1),
                                              nrow = 1))+
  scale_size_continuous(range = c(1, 2.5)) +
  theme_meta() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 11))
plot(p1)

# save the figure as a .rds file
saveRDS(object = p1, file = paste0("Figures/fig_6", fig, ".rds"))


# predict the slope

# all predictors normally distributed
cor_PC1 <- cor.test(gam_dat$Slope, gam_dat$PC1, method="pearson") # r = 0.29; P=0.42
cor_PC2 <- cor.test(gam_dat$Slope, gam_dat$PC2, method="pearson") # r = -0.50; P=0.14
cor_logprom <- cor.test(gam_dat$Slope, gam_dat$logprom, method="pearson") # r = -0.48; P=0.16

# get a correlation data.frame for adding the correlations to the plots
cor_list <- list(cor_PC1, cor_PC2, cor_logprom)

# figure 7 correlation plots

# get correct labels
if(taxa == "all") {
  title <- c("All taxa", NA, NA)
  fig <- c("a", "d", "g")
  ylab <- c(" ", "Local BEF slope (+-SE)", " ")
  xlab <- rep(" ", 3)
} else if(taxa == "active") {
  title <- c("Active dispersers", NA, NA)
  fig <- c("b", "e", "h")
  ylab <- rep(NA, 3)
  xlab <- c("PC1", "PC2", "ln( prominence ) (m)")
} else if(taxa == "passive") {
  title <- c("Passive dispersers", NA, NA)
  fig <- c("c", "f", "i")
  ylab <- rep(NA, 3)
  xlab <- rep(" ", 3)
}

# get the standard errors of the slopes
gam_dat$SE_low <- gam_dat$Slope - gam_dat$Slope_SE
gam_dat$SE_upp <- gam_dat$Slope + gam_dat$Slope_SE

# set the x variables
x_vars <- c("PC1", "PC2", "logprom")

# loop over the three different varialbes
for(i in 1:length(x_vars)) {
  
  p <- 
    ggplot(data = gam_dat) +
    geom_errorbar(aes_string(x = x_vars[i], ymin = "SE_low", ymax = "SE_upp", colour = "Inselberg"),
                  width = 0, alpha = 0.85, show.legend = FALSE) +
    geom_point(aes_string(x = x_vars[i], y = "Slope", colour = "Inselberg"), 
               shape = 16, stroke = 0.5, alpha = 0.9,
               size = 2.5) +
    ggtitle(if(!is.na(title[i])){title[i]}else{NULL}) +
    ylab(if(!is.na(ylab[i])){ylab[i]}else{NULL}) +
    xlab(xlab[i]) +
    scale_colour_manual(values = wes_palette(name = "Darjeeling1", n = 10, type = "continuous")) +
    scale_fill_manual(values = wes_palette(name = "Darjeeling1", n = 10, type = "continuous")) +
    guides(size = "none", colour = guide_legend(override.aes = list(size = 2.5, alpha = 1, shape = 16),
                                                nrow = 1)) +
    annotate(geom = "text", 
             label = paste0("r = ", round(cor_list[[i]]$estimate, 2), " ; ",
                            "P = ", round(cor_list[[i]]$p.value, 2)),
             x = -Inf, y = -Inf, vjust = -1, hjust = -0.2, size = 3.5) +
    theme_meta() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 11))
  
  saveRDS(object = p, file = paste0("Figures/fig_7", fig[i], ".rds"))
  
}

### END

