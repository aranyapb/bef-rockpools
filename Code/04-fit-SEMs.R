#'
#' @title SEM analysis
#' 
#' @description Fit structural equation models to the different sets of taxa i.e.
#' all taxa, active dispersers and passive dispersers
#'

# load relevant libraries
library(lme4)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggbeeswarm)
library(piecewiseSEM)
library(car)

# load plotting theme
source("Code/helper-plotting-theme.R")
source("Code/helper-plot-functions.R")
source("Code/helper-ci-lmm.R")
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


# model with direct effect of gamma on biomass

# fit the structural equation model using piecewiseSEM

# fit each model individually to check assumptions and convergence: all converge

# model 1
lm1 <- lm(loggamma ~ PC1 + PC2, data = gam_dat)

# test for residual normality
shapiro.test(resid(lm1)) 

# test for homogeneity of variance
ncvTest(lm1)

# plot the predicted values versus the residuals
plot(predict(lm1), resid(lm1))

# model 2: https://stats.stackexchange.com/questions/615331/should-i-pool-multiple-observations-from-the-same-experimental-unit-or-use-mixe
lm2 <- lmer(logalpha ~ logdepth + loggamma + (1|Inselberg), data = alph_dat)

# add observation random effect
lm2a <- lmer(logalpha ~ (1|obs), data = alph_dat, REML = TRUE)
alph_dat$obs = factor(1:nrow(alph_dat))

# test for residual normality
shapiro.test(resid(lm2)) 

# test for homogeneity of variance
ncvTest(lm(resid(lm2)~predict(lm2)))

# plot the predicted values versus the residuals
plot(predict(lm2), resid(lm2))

# model 3: https://stats.stackexchange.com/questions/615331/should-i-pool-multiple-observations-from-the-same-experimental-unit-or-use-mixe
lm3 <- lmer(logbiomass ~ loggamma + PC1 + PC2 + logalpha + logdepth + (1|Inselberg), data = alph_dat)

# test for residual normality
shapiro.test(resid(lm3)) 

# test for homogeneity of variance
ncvTest(lm(resid(lm3)~predict(lm3)))

# plot the predicted values versus the residuals
plot(predict(lm3), resid(lm3))

# combine models with logalpha into a piecewise SEM
SEM <- psem(
  lm(loggamma ~ PC1 + PC2, data = gam_dat),
  lmer(logalpha ~ logdepth + loggamma + (1|Inselberg), data = alph_dat), 
  lmer(logbiomass ~ loggamma + PC1 + PC2 + logalpha + logdepth + (1|Inselberg), data = alph_dat)
  )
SEM_sum <- summary(SEM, conditioning = T, standardize = "scale") 

# model convergence warning can be ignored since individual models 
# do converge: https://github.com/jslefche/piecewiseSEM/issues/242

# check the output 
print(SEM_sum)

# plot the SEMs
plot(SEM,show = "std",digits=2)


# conditional effect plots: lm3

# get summary object of the fitted model
lm3_sum <- summary(lm3)
print(lm3_sum)

# extracting the coefficients and giving them a name
for(i in 1:length(lm3@beta)) {
  assign(c("a", paste0("b", 1:(length(lm3@beta)-1) ))[i], lm3@beta[i])
}

# predict the logbiomass with the model formula manually
alph_dat$logbiomass_pred <- 
  with(alph_dat,
       a + (b1*loggamma) + (b2*PC1) + (b3*PC2) + (b4*logalpha) + (b5*logdepth))

# appendix: gamma -> BM
# fig. s1.2

# make predictions for each inselberg
min_pred <- ( min(alph_dat$loggamma) - 0.1)
max_pred <- ( max(alph_dat$loggamma) + 0.1)
gam_pred <- 
  expand.grid(loggamma = seq(min_pred, max_pred, 0.05),
              logalpha = mean(alph_dat$logalpha),
              logdepth = mean(alph_dat$logdepth),
              PC1 = mean(alph_dat$PC1),
              PC2 = mean(alph_dat$PC2))
str(gam_pred)

# use the model to get the mean prediction
gam_pred$logbiomass <- 
  with(gam_pred,
       a + (b1*loggamma) + (b2*PC1) + (b3*PC2) + (b4*logalpha) + (b5*logdepth) )

# calculate the 95% confidence interval
CI_df <- ci_lmm(obs_data = alph_dat, pred_data = gam_pred, obs_resp = "logbiomass", pred = "logbiomass_pred")

# calculate the confidence intervals
gam_pred$CI_low <- gam_pred$logbiomass + CI_df[["CI_low"]]
gam_pred$CI_upp <- gam_pred$logbiomass + CI_df[["CI_upp"]]

p1 <- 
  cond_eff_plot(data_raw = alph_dat, data_pred = gam_pred, 
                x_var = "loggamma", y_var = "logbiomass", 
                xlab = gamma_div, ylab = ln_biomass, size_var = "Gamma_SE",
                taxa = taxa, labels = c("a", "b", "c") )

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_s1.2", p1[[2]], ".rds"))


# model with total effect of gamma on biomass

# simply remove alpha from gamma - biomass model

# combine models with logalpha into a piecewise SEM
SEM <- psem(
  lm(loggamma ~ PC1 + PC2, data = gam_dat),
  lmer(logalpha ~ logdepth + loggamma + (1|Inselberg), data = alph_dat), 
  lmer(logbiomass ~ loggamma + PC1 + PC2 + logdepth + (1|Inselberg), data = alph_dat)
)

SEM_sum <- summary(SEM, conditioning = T, standardize = "scale") 

# model convergence warning can be ignored since individual models 
# do converge: https://github.com/jslefche/piecewiseSEM/issues/242

# check the output 
print(SEM_sum)

# plot the SEMs
plot(SEM,show = "std",digits=2)

# conditional effect plots: lm3

# model 3: https://stats.stackexchange.com/questions/615331/should-i-pool-multiple-observations-from-the-same-experimental-unit-or-use-mixe
lm3 <- lmer(logbiomass ~ loggamma + PC1 + PC2 + logdepth + (1|Inselberg), data = alph_dat)

# get summary object of the fitted model
lm3_sum <- summary(lm3)
print(lm3_sum)

# extracting the coefficients and giving them a name
for(i in 1:length(lm3@beta)) {
  assign(c("a", paste0("b", 1:(length(lm3@beta)-1) ))[i], lm3@beta[i])
}

# predict the logbiomass with the model formula manually
alph_dat$logbiomass_pred <- 
  with(alph_dat,
       a + (b1*loggamma) + (b2*PC1) + (b3*PC2) + (b4*logdepth))

# main text: PC1 -> BM
# fig. 5

# make predictions for each inselberg
min_pred <- ( min(alph_dat$loggamma) - 0.1)
max_pred <- ( max(alph_dat$loggamma) + 0.1)
gam_pred <- 
  expand.grid(loggamma = seq(min_pred, max_pred, 0.05),
              logdepth = mean(alph_dat$logdepth),
              PC1 = mean(alph_dat$PC1),
              PC2 = mean(alph_dat$PC2))
str(gam_pred)

# use the model to get the mean prediction
gam_pred$logbiomass <- 
  with(gam_pred,
       a + (b1*loggamma) + (b2*PC1) + (b3*PC2) + (b4*logdepth) )

# calculate the 95% confidence interval
CI_df <- ci_lmm(obs_data = alph_dat, pred_data = gam_pred, obs_resp = "logbiomass", pred = "logbiomass_pred")

# calculate the confidence intervals
gam_pred$CI_low <- gam_pred$logbiomass + CI_df[["CI_low"]]
gam_pred$CI_upp <- gam_pred$logbiomass + CI_df[["CI_upp"]]

p1 <- 
  cond_eff_plot(data_raw = alph_dat, data_pred = gam_pred, 
                x_var = "loggamma", y_var = "logbiomass", size_var = "Gamma_SE",
                xlab = gamma_div, ylab = ln_biomass, 
                taxa = taxa, labels = c("a", "b", "c") )

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_5", p1[[2]], ".rds"))

# appendix: PC1 -> BM
# figure s1.4

# make predictions for each inselberg
min_pred <- ( min(alph_dat$PC1) - 0.1)
max_pred <- ( max(alph_dat$PC1) + 0.1)
PC1_pred <- 
  expand.grid(PC1 = seq(min_pred, max_pred, 0.05),
              logalpha = mean(alph_dat$logalpha),
              logdepth = mean(alph_dat$logdepth),
              loggamma = mean(alph_dat$loggamma),
              PC2 = mean(alph_dat$PC2))
print(PC1_pred)

# use the model to get the mean prediction
PC1_pred$logbiomass <- 
  with(PC1_pred,
       a + (b1*loggamma) + (b2*PC1) + (b3*PC2) + (b4*logalpha) + (b5*logdepth) )

# calculate the 95% confidence interval
CI_df <- ci_lmm(obs_data = alph_dat, pred_data = PC1_pred, obs_resp = "logbiomass", pred = "logbiomass_pred")

# calculate the confidence intervals
PC1_pred$CI_low <- PC1_pred$logbiomass + CI_df[["CI_low"]]
PC1_pred$CI_upp <- PC1_pred$logbiomass + CI_df[["CI_upp"]]

# plot pc1-bm relationship for supplementary
p1 <- 
  cond_eff_plot(data_raw = alph_dat, data_pred = PC1_pred, 
                x_var = "PC1", y_var = "logbiomass", 
                xlab = "PC1", ylab = ln_biomass, 
                taxa = taxa, labels = c("a", "b", "c") )

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_s1.4", p1[[2]], ".rds"))


# appendix: PC2 - BM
# figure s1.4

# make predictions for each inselberg
min_pred <- ( min(alph_dat$PC2) - 0.1)
max_pred <- ( max(alph_dat$PC2) + 0.1)
PC2_pred <- 
  expand.grid(PC2 = seq(min_pred, max_pred, 0.05),
              logalpha = mean(alph_dat$logalpha),
              logdepth = mean(alph_dat$logdepth),
              loggamma = mean(alph_dat$loggamma),
              PC1 = mean(alph_dat$PC1))
print(PC2_pred)

# use the model to get the mean prediction
PC2_pred$logbiomass <- 
  with(PC2_pred,
       a + (b1*loggamma) + (b2*PC1) + (b3*PC2) + (b4*logalpha) + (b5*logdepth) )

# calculate the 95% confidence interval
CI_df <- ci_lmm(obs_data = alph_dat, pred_data = PC2_pred, obs_resp = "logbiomass", pred = "logbiomass_pred")

# calculate the confidence intervals
PC2_pred$CI_low <- PC2_pred$logbiomass + CI_df[["CI_low"]]
PC2_pred$CI_upp <- PC2_pred$logbiomass + CI_df[["CI_upp"]]

# plot pc2-bm relationship for supplementary
p1 <- 
  cond_eff_plot(data_raw = alph_dat, data_pred = PC2_pred, 
                x_var = "PC2", y_var = "logbiomass", 
                xlab = "PC2", ylab = ln_biomass, 
                taxa = taxa, labels = c("d", "e", "f") )

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_s1.4", p1[[2]], ".rds"))


# appendix: depth - BM
# figure s1.3

# make predictions for each inselberg
min_pred <- ( min(alph_dat$Depth) - 0.1)
max_pred <- ( max(alph_dat$Depth) + 0.1)
depth_pred <- 
  expand.grid(logdepth = seq(0.60, 4.6, 0.05),
              logalpha = mean(alph_dat$logalpha),
              PC1 = mean(alph_dat$PC1),
              loggamma = mean(alph_dat$loggamma),
              PC2 = mean(alph_dat$PC2))

# use the model to get the mean prediction
depth_pred$logbiomass <- 
  with(depth_pred,
       a + (b1*loggamma) + (b2*PC1) + (b3*PC2) + (b4*logalpha) + (b5*logdepth) )

# calculate the 95% confidence interval
CI_df <- ci_lmm(obs_data = alph_dat, pred_data = depth_pred, obs_resp = "logbiomass", pred = "logbiomass_pred")

# calculate the confidence intervals
depth_pred$CI_low <- depth_pred$logbiomass + CI_df[["CI_low"]]
depth_pred$CI_upp <- depth_pred$logbiomass + CI_df[["CI_upp"]]

# plot depth-bm relationship for supplementary
p1 <- 
  cond_eff_plot(data_raw = alph_dat, data_pred = depth_pred, 
                x_var = "logdepth", y_var = "logbiomass", 
                xlab = ln_depth, ylab = ln_biomass, 
                taxa = taxa, labels = c("a", "b", "c") )

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_s1.3", p1[[2]], ".rds"))


# conditional effect plots: lm2

# get summary object of the fitted model
lm2_sum <- summary(lm2)
print(lm2_sum)

# extracting the coefficients and giving them a name
for(i in 1:length(lm2@beta)) {
  assign(c("a", paste0("b", 1:(length(lm2@beta)-1) ))[i], lm2@beta[i])
}

# predict logalpha using the model
alph_dat$logalpha_pred <- 
  with(alph_dat, 
       a + (b1*logdepth) + (b2*loggamma))


# appendix: gamma -> alpha
# fig. s1.6

# make predictions for each inselberg
min_pred <- ( min(alph_dat$loggamma) - 0.1)
max_pred <- ( max(alph_dat$loggamma) + 0.1)
gamma_pred <- 
  expand.grid(loggamma = seq(min_pred, max_pred, 0.05),
              logdepth = mean(alph_dat$logdepth))

# use the model to get the mean prediction
gamma_pred$logalpha <- 
  with(gamma_pred,
       a + (b1*logdepth) + (b2*loggamma))

# calculate the 95% confidence interval
CI_df <- ci_lmm(obs_data = alph_dat, pred_data = gamma_pred, obs_resp = "logalpha", pred = "logalpha_pred")

# calculate the confidence intervals
gamma_pred$CI_low <- gamma_pred$logalpha + CI_df[["CI_low"]]
gamma_pred$CI_upp <- gamma_pred$logalpha + CI_df[["CI_upp"]]

# plot the gamma alpha relationship
p1 <- 
  cond_eff_plot(data_raw = alph_dat, data_pred = gamma_pred, 
                x_var = "loggamma", y_var = "logalpha", 
                xlab = gamma_div, ylab = alpha_div,
                taxa = taxa, labels = c("a", "b", "c") )

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_s1.6", p1[[2]], ".rds"))


# appendix: depth - alpha
# fig. s1.3

# make predictions for each inselberg
min_pred <- ( min(alph_dat$logdepth) - 0.1)
max_pred <- ( max(alph_dat$logdepth) + 0.1)
alpha_pred <- 
  expand.grid(logdepth = seq(min_pred, max_pred, 0.05),
              loggamma = mean(alph_dat$loggamma))

# use the model to get the mean prediction
alpha_pred$logalpha <- 
  with(alpha_pred,
       a + (b1*logdepth) + (b2*loggamma))

# calculate the 95% confidence interval
CI_df <- ci_lmm(obs_data = alph_dat, pred_data = alpha_pred, obs_resp = "logalpha", pred = "logalpha_pred")

# calculate the confidence intervals
alpha_pred$CI_low <- alpha_pred$logalpha + CI_df[["CI_low"]]
alpha_pred$CI_upp <- alpha_pred$logalpha + CI_df[["CI_upp"]]

# plot the gamma alpha relationship
p1 <- 
  cond_eff_plot(data_raw = alph_dat, data_pred = alpha_pred, 
                x_var = "logdepth", y_var = "logalpha", size_var = "Alpha_SE",
                xlab = ln_depth, ylab = alpha_div,
                taxa = taxa, labels = c("d", "e", "f") )

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_s1.3", p1[[2]], ".rds"))


# conditional effect plots: lm1

# get summary object of the fitted model
lm1_sum <- summary(lm1)
print(lm1_sum)

# extract model coefficients
lm1_beta <- coef(lm1)

# extracting the coefficients and giving them a name
for(i in 1:length(lm1_beta)) {
  assign(c("a", paste0("b", 1:(length(lm1_beta)-1) ))[i], lm1_beta[i])
}

# predict loggamma using the model
gam_dat$loggamma_pred <- 
  with(gam_dat, 
       a + (b1*PC1) + (b2*PC2))

# appendix: PC1 - gamma
# fig. s1.5

# make predictions for each inselberg
min_pred <- ( min(gam_dat$PC1) - 0.1)
max_pred <- ( max(gam_dat$PC1) + 0.1)
gam_pred <- 
  expand.grid(PC1 = seq(min_pred, max_pred, 0.05),
              PC2 = mean(gam_dat$PC2))

# use the model to get the mean prediction
gam_pred$loggamma <- 
  with(gam_pred,
       a + (b1*PC1) + (b2*PC2))

# calculate the 95% confidence interval
CI_df <- ci_lmm(obs_data = gam_dat, pred_data = gam_pred, obs_resp = "loggamma", pred = "loggamma_pred")

# calculate the confidence intervals
gam_pred$CI_low <- gam_pred$loggamma + CI_df[["CI_low"]]
gam_pred$CI_upp <- gam_pred$loggamma + CI_df[["CI_upp"]]

# plot the PC1 - gamma relationship
p1 <- 
  cond_eff_plot(data_raw = gam_dat, data_pred = gam_pred, 
                x_var = "PC1", y_var = "loggamma", size_var = "Gamma_SE",
                xlab = "PC1", ylab = gamma_div,
                taxa = taxa, labels = c("a", "b", "c") )

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_s1.5", p1[[2]], ".rds"))

# appendix: PC2 - gamma
# fig. s1.5

# make predictions for each inselberg
min_pred <- ( min(gam_dat$PC2) - 0.1)
max_pred <- ( max(gam_dat$PC2) + 0.1)
gam_pred <- 
  expand.grid(PC2 = seq(min_pred, max_pred, 0.05),
              PC1 = mean(gam_dat$PC1))

# use the model to get the mean prediction
gam_pred$loggamma <- 
  with(gam_pred,
       a + (b1*PC1) + (b2*PC2))

# calculate the 95% confidence interval
CI_df <- ci_lmm(obs_data = gam_dat, pred_data = gam_pred, obs_resp = "loggamma", pred = "loggamma_pred")

# calculate the confidence intervals
gam_pred$CI_low <- gam_pred$loggamma + CI_df[["CI_low"]]
gam_pred$CI_upp <- gam_pred$loggamma + CI_df[["CI_upp"]]

# plot the PC2 - gamma relationship
p1 <- 
  cond_eff_plot(data_raw = gam_dat, data_pred = gam_pred, 
                x_var = "PC2", y_var = "loggamma", size_var = "Gamma_SE",
                xlab = "PC2", ylab = gamma_div,
                taxa = taxa, labels = c("d", "e", "f") )

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_s1.5", p1[[2]], ".rds"))

### END

