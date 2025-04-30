
# test hypothesis 1

# prepare the environment
taxa = "all"
hypothesis = "1"
source("Code/helper-prep-env.R")

# specify the causal hypothesis via a dag

# causal hypothesis 1:

# create the model using dagitty
dag_h1a <- get_dag_h1a()

# plot the model
plot(dag_h1a)

# check the implied conditional independencies
impliedConditionalIndependencies(dag_h1a)

# test the causal hypothesis via the conditional independencies

# 1. test unconditional independence: Depth ⊥ PC1
test1 <- cor.test(dat$logdepth, dat$PC1)

# 2. test unconditional independence: Depth ⊥ PC2
test2 <- cor.test(dat$logdepth, dat$PC2)

# 3. test conditional independence: PC1 ⊥ PC2
test3 <- cor.test(dat$PC1, dat$PC2)

# extract the p-value
test_p <- sapply(list(test1, test2, test3), function(x) x[["p.value"]])

# run a bonferroni correction
test_p_adj <- p.adjust(test_p, method = "bonferroni")
round(test_p_adj, 6)

# run the fisher's combined test
fisher_combined_test(p_values = test_p)

# revise the causal hypothesis to account for the failed tests

# causal hypothesis 2:

# create the model using dagitty
dag_h1b <- get_dag_h1b()

# plot the model
plot(dag_h1b)

# check the implied conditional independencies
impliedConditionalIndependencies(dag_h1b)

# test the revised causal hypothesis via the conditional independencies

# 1. test unconditional independence: Depth ⊥ PC1
test1 <- cor.test(dat$logdepth, dat$PC1)

# 2. test conditional independence: PC1 ⊥ PC2
test2 <- cor.test(dat$PC1, dat$PC2)

# extract the p-value
test_p <- sapply(list(test1, test2), function(x) x[["p.value"]])

# run a bonferroni correction
test_p_adj <- p.adjust(test_p, method = "bonferroni")
round(test_p_adj, 6)

# run the fisher's combined test
fisher_combined_test(p_values = test_p)

# analysis

# hypothesis 1:

# total effect of gamma on biomass

# adjustment set: PC1, PC2
adjustmentSets(x = dag_h1b, exposure = "Gamma", outcome = "Biomass", effect = "total")

# h1_dat
h1_dat <-
  dat |>
  dplyr::select(logbiomass, loggamma, logdepth, PC1, PC2, Inselberg, Gamma_SE) |>
  dplyr::mutate(logbiomass = scale(logbiomass)[, 1],
                loggamma = scale(loggamma)[, 1],
                logdepth = scale(logdepth)[, 1],
                PC1 = scale(PC1)[, 1],
                PC2 = scale(PC2)[, 1])

# fit the relevant linear mixed models with and without splines in the covariates

# create list of model formulas
model_formulas <- list(
  m1 = logbiomass ~ loggamma + logdepth + PC1 + PC2 + (1 | Inselberg),
  m2 = logbiomass ~ loggamma + bs(logdepth, degree = 2) + PC1 + PC2 + (1 | Inselberg),
  m3 = logbiomass ~ loggamma + logdepth + bs(PC1, degree = 2)  + PC2 + (1 | Inselberg),
  m4 = logbiomass ~ loggamma + logdepth + PC1 + bs(PC2, degree = 2)  + (1 | Inselberg),
  m5 = logbiomass ~ loggamma + bs(logdepth, degree = 2)  + bs(PC1, degree = 2) + PC2 + (1 | Inselberg),
  m6 = logbiomass ~ loggamma + logdepth + bs(PC1, degree = 2) + bs(PC2, degree = 2) + (1 | Inselberg),
  m7 = logbiomass ~ loggamma + bs(logdepth, degree = 2) + PC1 + bs(PC2, degree = 2) + (1 | Inselberg),
  m8 = logbiomass ~ loggamma + bs(logdepth, degree = 2) + bs(PC1, degree = 2) + bs(PC2, degree = 2) + (1 | Inselberg),
  m9 = logbiomass ~ loggamma + bs(logdepth, degree = 3) + bs(PC1, degree = 3) + bs(PC2, degree = 3) + (1 | Inselberg)
)

# fit and compare the models
model_comparison <- 
  fit_and_compare_lmm_models(formulas = model_formulas, 
                             data = h1_dat, n_sim = 1000)

# check the model outputs
model_comparison$comparison

# m5 converges, meets the assumptions (from DHARMA) and has the lowest AIC

# fit the m5 model to test hypothesis 1
lm1 <- lmer(logbiomass ~ loggamma + bs(logdepth, degree = 2)  + bs(PC1, degree = 2) + PC2 + (1 | Inselberg),
            data = h1_dat,
            REML = FALSE)

# plot the model fit
dev.off()
plot(predict(lm1), h1_dat$logbiomass)
abline(0, 1)

# check the rsquared value
piecewiseSEM::rsquared(lm1)

# check assumptions using DHARMa

# simulate residuals
sim_res <- simulateResiduals(fittedModel = lm1, n = 2000)

# create diagnostic plots

# qq plot residuals
plotQQunif(sim_res)

# residuals versus predicted values
plotResiduals(sim_res)

# test for overdispersion
testDispersion(sim_res)

# hypothesis test

# use a likelihood ratio test to test hypothesis 1
# loggamma parameter is non-significant
# no total effect of loggamma on logbiomass
drop1(lm1, test = "Chisq")

# plot fig 5

# get model predictions
lm1_pred <-
  predict_lmer_with_intervals_merTools(model = lm1, 
                                       data = h1_dat, 
                                       predictors = c("loggamma", "PC1", "PC2", "logdepth"),
                                       vary = "loggamma", 
                                       n_points = 100, nsim = 1000,
                                       include_re = FALSE,
                                       group_levels = NULL,
                                       seed = 123)

# plot the data with the model predictions
p1 <- 
  cond_eff_plot(data_raw = h1_dat, data_pred = lm1_pred, 
                x_var = "loggamma", y_var = "logbiomass", 
                xlab = gamma_div, ylab = ln_biomass, size_var = "Gamma_SE",
                taxa = taxa, labels = c("a", "b", "c") )
p1

# save the figure as a .rds file
saveRDS(object = p1[[1]], file = paste0("Figures/fig_5", p1[[2]], ".rds"))





