#'
#' @title SEM analysis
#' 
#' @description Fit structural equation models to the different sets of taxa i.e.
#' all taxa, active dispersers and passive dispersers
#'

# load relevant libraries
library(lme4)
library(dagitty)
library(splines)
library(DHARMa)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggbeeswarm)
library(piecewiseSEM)
library(car)
library(ppcor)

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

# convert to a dat
dat <- alph_dat

# check for missing values
any(is.na(dat))

# check the structure of the data
str(dat)

# set-up frequently used axis labels
alpha_div <- expression("ln("~alpha~"-diversity )")
gamma_div <- expression("ln("~gamma~"-diversity )")
ln_biomass <- "ln( biomass ) (mg)"
ln_depth <- "ln (depth) (cm)"

# specify the causal hypothesis

# causal hypothesis 1:

# create the model using dagitty
dag1 <- dagitty('
dag {
  bb="0,0,1,1"
  Alpha [pos="0.473,0.490"]
  Biomass [pos="0.600,0.270"]
  Depth [pos="0.603,0.489"]
  Gamma [pos="0.340,0.273"]
  PC1 [pos="0.469,0.136"]
  PC2 [pos="0.469,0.341"]
  Alpha -> Biomass
  Depth -> Alpha
  Depth -> Biomass
  Gamma -> Alpha
  Gamma -> Biomass
  PC1 -> Biomass
  PC1 -> Gamma
  PC2 -> Biomass
  PC2 -> Gamma
}
')

# plot the model
plot(dag1)

# check the implied conditional independencies
impliedConditionalIndependencies(dag1)

# test the causal hypothesis via the conditional independencies

# 1. test unconditional independence: PC1 ⊥ PC2
test1 <- cor.test(dat$PC1, dat$PC2)

# 2. test conditional independence: PC1 ⊥ Alpha | Gamma
test2 <- pcor.test(dat$PC1, dat$Alpha, dat$Gamma, method = "pearson")

# 3. test unconditional independence: PC1 ⊥ Depth
test3 <- cor.test(dat$PC1, dat$Depth)

# 4. test conditional independence: PC2 ⊥ Alpha | Gamma
test4 <- pcor.test(dat$PC2, dat$Alpha, dat$Gamma, method = "pearson")

# 5. test unconditional independence: PC2 ⊥ Depth
test5 <- cor.test(dat$PC2, dat$Depth)

# 6. test unconditional independence: Gamma ⊥ Depth
test6 <- cor.test(dat$Gamma, dat$Depth)

# pull results into a list
test_list <- list(
  test1,
  test2,
  test3,
  test4,
  test5,
  test6
)

# extract the p-value
test_p <- sapply(test_list, function(x) x[["p.value"]])

# run a bonferroni correction
test_p_adj <- p.adjust(test_p, method = "bonferroni")
round(test_p_adj, 6)

# apply Fisher's combined probability test manually
fisher_stat <- -2 * sum(log(test_p))

# degrees of freedom = 2 * number of p-values
df <- 2 * length(test_p)

# get p-value from chi-squared distribution
fisher_pval <- pchisq(fisher_stat, df = df, lower.tail = FALSE)
round(fisher_pval, 6)

# the d-separation tests failed
# connection between PC2 and alpha is strong (correlation: 0.30)
# theory suggests that PC2 could affect alpha (i.e. climate could affect alpha)

# revise the causal hypothesis to account for the failed tests

# causal hypothesis 2:

# create the model using dagitty
dag2 <- dagitty('
dag {
  bb="0,0,1,1"
  Alpha [pos="0.473,0.490"]
  Biomass [pos="0.600,0.270"]
  Depth [pos="0.603,0.489"]
  Gamma [pos="0.340,0.273"]
  PC1 [pos="0.469,0.136"]
  PC2 [pos="0.469,0.341"]
  Alpha -> Biomass
  Depth -> Alpha
  Depth -> Biomass
  Gamma -> Alpha
  Gamma -> Biomass
  PC1 -> Biomass
  PC1 -> Gamma
  PC2 -> Alpha
  PC2 -> Biomass
  PC2 -> Gamma
}
')

# plot the model
plot(dag2)

# check the implied conditional independencies
impliedConditionalIndependencies(dag2)

# test the revised causal hypothesis via the conditional independencies

# 1. test unconditional independence: PC1 ⊥ PC2
test1 <- cor.test(dat$PC1, dat$PC2)

# 2. test conditional independence: PC1 ⊥ Alpha | Gamma, PC2
test2 <- pcor.test(dat$PC1, dat$Alpha, dat[,c("Gamma", "PC2")], method = "pearson")

# 3. test unconditional independence: PC1 ⊥ Depth
test3 <- cor.test(dat$PC1, dat$Depth)

# 4. test unconditional independence: PC2 ⊥ Depth
test4 <- cor.test(dat$PC2, dat$Depth)

# 5. test unconditional independence: Gamma ⊥ Depth
test5 <- cor.test(dat$Gamma, dat$Depth)

# pull results into a list
test_list <- list(
  test1,
  test2,
  test3,
  test4,
  test5
)

# extract the p-value
test_p <- sapply(test_list, function(x) x[["p.value"]])
round(test_p, 6)

# run a bonferroni correction
test_p_adj <- p.adjust(test_p, method = "bonferroni")
round(test_p_adj, 6)

# apply Fisher's combined probability test manually
fisher_stat <- -2 * sum(log(test_p))

# degrees of freedom = 2 * number of p-values
df <- 2 * length(test_p)

# get p-value from chi-squared distribution
fisher_pval <- pchisq(fisher_stat, df = df, lower.tail = FALSE)
round(fisher_pval, 6)

# analysis

# hypothesis 1:

# total effect of gamma on biomass

# adjustment set: PC1, PC2
adjustmentSets(x = dag2, exposure = "Gamma", outcome = "Biomass")

# h1_dat
h1_dat <-
  dat |>
  dplyr::select(logbiomass, loggamma, PC1, PC2, Inselberg, Gamma_SE) |>
  dplyr::mutate(logbiomass = scale(logbiomass)[, 1],
                loggamma = scale(loggamma)[, 1],
                PC1 = scale(PC1)[, 1],
                PC2 = scale(PC2)[, 1])

# fit the relevant linear mixed model without any splines
lm1 <- lmer(logbiomass ~ loggamma + PC1 + PC2 + (1 | Inselberg), data = h1_dat)

# check model assumptions using DHARMa

# simulate residuals
sim_res <- simulateResiduals(fittedModel = lm1, n = 1000)

# qq plot residuals
plotQQunif(sim_res)

# residuals versus predicted values
plotResiduals(sim_res)

# residual versus predicted values suggests problems with non-linearity
# run a model comparison using basis splines

# fit the relevant linear mixed models with splines

# create list of model formulas
model_formulas <- list(
  m1 = logbiomass ~ loggamma + PC1 + PC2 + (1 | Inselberg),
  m2 = logbiomass ~ bs(loggamma, degree = 2) + PC1 + PC2 + (1 | Inselberg),
  m3 = logbiomass ~ loggamma + bs(PC1, degree = 2) + PC2 + (1 | Inselberg),
  m4 = logbiomass ~ loggamma + PC1 + bs(PC2, degree = 2) + (1 | Inselberg),
  m5 = logbiomass ~ loggamma + bs(PC1, degree = 2) + bs(PC2, degree = 2) + (1 | Inselberg),
  m6 = logbiomass ~ bs(loggamma, degree = 2) + bs(PC1, degree = 2) + bs(PC2, degree = 2) + (1 | Inselberg)
)

# fit models using a loop
model_fits <- list()
model_checks <- list()
for (name in names(model_formulas)) {
  # fit the model and add to output list
  model_fits[[name]] <- lmer(
    formula = model_formulas[[name]],
    data = h1_dat,
    REML = FALSE
  )
  # check model assumptions using DHARMa
  sim_res <- simulateResiduals(fittedModel = model_fits[[name]] , n = 1000)
  # create the plot
  model_checks[[name]] <- sim_res
}

# compare models by AIC
model_aics <- sapply(model_fits, AIC)
model_comparison <- data.frame(
  Model = names(model_aics ),
  AIC = as.numeric(model_aics )
)

# order the models
model_comparison <- model_comparison[order(model_comparison$AIC), ]
model_comparison

# check the graphs to see which models meet relevant assumptions
sapply(model_checks, plot)

# m4 and m6 met the assumptions
# m4 and m6 have similar AIC values
# use m4 as it is the simpler model

# fit the m4 model to test hypothesis 1
lm1 <- lmer(logbiomass ~ loggamma + PC1 + bs(PC2, degree = 2) + (1 | Inselberg), data = h1_dat)

# plot the model fit
plot(predict(lm1), h1_dat$logbiomass)
abline(0, 1)

# check assumptions using DHARMa

# simulate residuals
lm1_sim_res <- simulateResiduals(fittedModel = lm1, n = 2000)

# create diagnostic plots

# qq plot residuals
plotQQunif(lm1_sim_res)

# residuals versus predicted values
plotResiduals(lm1_sim_res)

# test for overdispersion
testDispersion(lm1_sim_res)

# use a likelihood ratio test to test hypothesis 1
# loggamma parameter is non-significant
# no total effect of loggamma on logbiomass
drop1(lm1, test = "Chisq")

# plot fig 5

# get model predictions
lm1_pred <-
  predict_lmer_with_intervals(model = lm1, data = h1_dat,
                              predictors = c("loggamma", "PC1", "PC2"),
                              vary = "loggamma",
                              n_points = 100, nsim = 1000)

# plot the data with the model predictions
p1 <- 
  cond_eff_plot(data_raw = h1_dat, data_pred = lm1_pred, 
                x_var = "loggamma", y_var = "logbiomass", 
                xlab = gamma_div, ylab = ln_biomass, size_var = "Gamma_SE",
                taxa = taxa, labels = c("a", "b", "c") )
p1

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_5", p1[[2]], ".rds"))


# hypothesis 2:

# direct effect of alpha on biomass and variation across inselbergs

# adjustment set: PC1, PC2
adjustmentSets(x = dag2, exposure = "Alpha", outcome = "Biomass", effect = "direct")

# h2_dat
h2_dat <-
  dat |>
  dplyr::select(logbiomass, loggamma, logalpha, logdepth, PC2, Inselberg, Alpha_SE) |>
  dplyr::mutate(logbiomass = scale(logbiomass)[, 1],
                loggamma = scale(loggamma)[, 1],
                logalpha = scale(logalpha)[, 1],
                PC2 = scale(PC2)[, 1],
                logdepth = scale(logdepth)[, 1])

# fit the relevant linear mixed model without any splines
lm2 <- lmer(logbiomass ~ loggamma + logdepth + PC2 + logalpha + (1 | Inselberg), 
            data = h2_dat, REML = FALSE)

drop1(lm2, test = "Chisq")
summary(lm2)

lm2a <- lmer(logbiomass ~ loggamma + logdepth + PC2 + logalpha + (logalpha | Inselberg), 
             data = h2_dat, REML = FALSE)
summary(lm2a)

anova(lm2, lm2a)

# check model assumptions using DHARMa

# simulate residuals
sim_res <- simulateResiduals(fittedModel = lm2a, n = 1000)

# qq plot residuals
plotQQunif(sim_res)

# residuals versus predicted values
plotResiduals(sim_res)

# use a likelihood ratio test to test hypothesis 2
# interaction inselberg:logalpha is significant
drop1(lm2, test = "Chisq")

# get model predictions


# create a data.frame where logalpha and inselberg vary but the rest are held constant
lm2_pred <- 
  h2_dat |>
  dplyr::select(Inselberg, logdepth, PC2, loggamma, logalpha) |>
  dplyr::mutate(logdepth = mean(alph_dat$logdepth),
                loggamma - mean(loggamma),
                PC2 = mean(PC2))

# # predict new logbiomass values but with constant logdepth, se.fit=T for SE
x <- dplyr::as_tibble(predict(lm2a, lm2_pred, interval = "confidence"))
lm2_pred$logbiomass <- x$fit
lm2_pred$CI_low <- x$lwr
lm2_pred$CI_upp <- x$upr

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
  geom_quasirandom(data = h2_dat, 
                   aes(x = logalpha, y = logbiomass, colour = Inselberg, size = Alpha_SE), 
                   width = 0.025, shape = 1, stroke = 0.5, alpha = 0.5) +
  geom_line(data = lm2_pred, aes(x = logalpha, y = logbiomass, colour = Inselberg),
            show.legend = FALSE) +
  geom_ribbon(data = lm2_pred, 
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
















# model 2: https://stats.stackexchange.com/questions/615331/should-i-pool-multiple-observations-from-the-same-experimental-unit-or-use-mixe
lm2 <- lmer(logalpha ~ logdepth + loggamma + (1|Inselberg), data = alph_dat)

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
p1[[1]]
p1[[2]]

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_s1.2", p1[[2]], ".rds"))

# model with total effect of gamma on biomass

# simply remove alpha from gamma - biomass model

# conditional effect plots: lm4

# model 4:
lm4 <- lmer(logbiomass ~ loggamma + PC1 + PC2 + (1|Inselberg), data = alph_dat)

# get summary object of the fitted model
lm4_sum <- summary(lm4)
print(lm4_sum)

# extracting the coefficients and giving them a name
for(i in 1:length(lm4@beta)) {
  assign(c("a", paste0("b", 1:(length(lm4@beta)-1) ))[i], lm4@beta[i])
}

# predict the logbiomass with the model formula manually
alph_dat$logbiomass_pred <- 
  with(alph_dat,
       a + (b1*loggamma) + (b2*PC1) + (b3*PC2))

# main text: PC1 -> BM
# fig. 5

# make predictions for each inselberg
min_pred <- ( min(alph_dat$loggamma) - 0.1)
max_pred <- ( max(alph_dat$loggamma) + 0.1)
gam_pred <- 
  expand.grid(loggamma = seq(min_pred, max_pred, 0.05),
              PC1 = mean(alph_dat$PC1),
              PC2 = mean(alph_dat$PC2))
str(gam_pred)

# use the model to get the mean prediction
gam_pred$logbiomass <- 
  with(gam_pred,
       a + (b1*loggamma) + (b2*PC1) + (b3*PC2))

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
p1

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

