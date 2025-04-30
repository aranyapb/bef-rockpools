
# test hypotheses 2 and 3

# load relevant libraries
library(lme4)
library(dagitty)
library(splines)
library(DHARMa)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(piecewiseSEM)

# load plotting theme
source("Code/helper-plotting-theme.R")
source("Code/helper-ci-lmm.R")
source("Code/helper-select-taxa.R")

# get a list of data files
files <- list.files("Data/")

# choose the set of taxa: all, active, passive
taxa <- select_taxa()

# load the alpha-scale data
dat <- readRDS(paste0("Data/", files[ grepl(pattern = paste0(taxa, "-taxa-alpha.rds"), files)] ))
head(dat)

# check for missing values
any(is.na(dat))

# check the structure of the data
str(dat)

# set-up frequently used axis labels
alpha_div <- expression("ln("~alpha~"-diversity )")
gamma_div <- expression("ln("~gamma~"-diversity )")
ln_biomass <- "ln( biomass ) (mg)"
ln_depth <- "ln (depth) (cm)"

## hypothesis 2

# specify the causal hypothesis

# causal hypothesis 1:

# create the model using dagitty
dag1 <- dagitty('
  dag {
  bb="0,0,1,1"
  Alpha [exposure,pos="0.492,0.246"]
  Biomass [outcome,pos="0.615,0.325"]
  Depth [pos="0.382,0.332"]
  Inselberg [pos="0.486,0.093"]
  Alpha -> Biomass
  Depth -> Alpha
  Depth -> Biomass
  Inselberg -> Alpha
  Inselberg -> Biomass
  Inselberg -> Depth
}')

# plot the model
plot(dag1)

# check the implied conditional independencies
impliedConditionalIndependencies(dag1)

# there are no conditional independencies but this causal structure is uncontroversial

# hypothesis 2:

# effect of alpha on biomass within inselbergs

# adjustment set: Depth, Inselberg
adjustmentSets(x = dag1, exposure = "Alpha", outcome = "Biomass", effect = "direct")

# h2_dat
h2_dat <-
  dat |>
  dplyr::select(logbiomass, logalpha, logdepth, Inselberg, Alpha_SE) |>
  dplyr::mutate(logbiomass = scale(logbiomass)[, 1],
                logalpha = scale(logalpha)[, 1],
                logdepth = scale(logdepth)[, 1])

# fit the relevant linear mixed model without any splines

# test hypothesis by comparing model without random slope to model with random slope

# build the models
model_formulas <- list(
  m1 = logbiomass ~ logalpha + logdepth + (1 | Inselberg),
  m2 = logbiomass ~ logalpha + logdepth + (1 + logalpha | Inselberg),
  m3 = logbiomass ~ logalpha + logdepth + (1 + logalpha + logdepth | Inselberg),
  m4 = logbiomass ~ logalpha + bs(logdepth, degree = 2) + (1 | Inselberg),
  m5 = logbiomass ~ logalpha + bs(logdepth, degree = 2) + (1 + logalpha | Inselberg)
)

# fit models using a loop
model_fits <- list()
model_checks <- list()

for (name in names(model_formulas)) {
  # fit the model and add to output list
  model_fits[[name]] <- lmer(
    formula = model_formulas[[name]],
    data = h2_dat,
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
  AIC = as.numeric(model_aics ),
  singular = sapply(model_fits, isSingular)
)

# order the models
model_comparison <- model_comparison[order(model_comparison$AIC), ]
model_comparison

# check the graphs to see which models meet relevant assumptions
for (i in seq_along(model_checks)) {
  plot(model_checks[[i]])
  title(main = paste("Model", i), line = 0.5)
}

# m4 and m5 do meet the assumptions although there is some deviation in m5
# however, this is probably justifiable
# indeed, the AIC of m5 is considerably better than m4

# fit the m4 and m5 models to perform the hypothesis test

# without random alpha slope: m4
lm2a <- lmer(logbiomass ~ logalpha + bs(logdepth, degree = 2) + (1 | Inselberg), 
             data = h2_dat,
             REML = FALSE)

# check model assumptions

# check for influential points
lm2a_dfbetas <- plot_lmer_dfbetas(lm2a, data = h2_dat)

# examine influential points
lm2a_inf <- apply(lm2a_dfbetas, 2, function(x) which(x == max(x)))
names(lm2a_inf) <- NULL

# print these data points
h2_dat[lm2a_inf, ]
dat[lm2a_inf, ]
summary(dat)

# DHARMa tests

# simulate residuals
sim_res <- simulateResiduals(fittedModel = lm2a, n = 1000)

# qq plot residuals
plotQQunif(sim_res)

# residuals versus predicted values
plotResiduals(sim_res)

# test for overdispersion
testDispersion(sim_res)

# check model fit
dev.off()
plot(fitted(lm2a), lm2a@frame[[1]])
abline(0, 1)

# check the r2 value
piecewiseSEM::rsquared(lm2a)

# with random alpha slope: m5
lm2b <- lmer(logbiomass ~ logalpha + bs(logdepth, degree = 2) + (1 + logalpha | Inselberg), 
             data = h2_dat,
             REML = FALSE)

# check model assumptions

# check for influential points
lm2b_dfbetas <- plot_lmer_dfbetas(lm2b, data = h2_dat)

# examine influential points
lm2b_inf <- apply(lm2b_dfbetas, 2, function(x) which(x == max(x)))
names(lm2b_inf) <- NULL

# print these data points
h2_dat[lm2b_inf, ]
dat[lm2b_inf, ]
summary(dat)

# DHARMa tests

# simulate residuals
sim_res <- simulateResiduals(fittedModel = lm2b, n = 1000)

# qq plot residuals
plotQQunif(sim_res)

# residuals versus predicted values
plotResiduals(sim_res)

# test for overdispersion
testDispersion(sim_res)

# check model fit
dev.off()
plot(fitted(lm2b), lm2b@frame[[1]])
abline(0, 1)

# check the r2 value
piecewiseSEM::rsquared(lm2b)

# use a likelihood ratio test to test hypothesis 2
# expect lm2b to fit the data better as it allows a varying slope
anova(lm2a, lm2b)

# plot the data
lm2_pred <-
  predict_lmer_with_intervals_merTools(model = lm2b, 
                                       data = h2_dat, 
                                       predictors = c("logalpha", "logdepth"),
                                       vary = "logalpha", 
                                       n_points = 100, nsim = 1000,
                                       include_re = TRUE,
                                       group_levels = list(Inselberg = as.character(unique(h2_dat$Inselberg))),
                                       seed = 123)

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
  geom_point(data = h2_dat, 
             aes(x = logalpha, y = logbiomass, colour = Inselberg, size = Alpha_SE), 
             shape = 1, stroke = 0.5, alpha = 0.5) +
  geom_line(data = lm2_pred, aes(x = logalpha, y = predicted_mean, colour = Inselberg),
            show.legend = FALSE) +
  geom_ribbon(data = lm2_pred, 
              aes(x = logalpha, ymax = upper, ymin = lower, fill = Inselberg), 
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


## hypothesis 3

# extract the random effects
lm2b_ranef <- ranef(lm2b)

# extract the slope offsets from lm2b
lm2b_slopes <- lm2b_ranef$Inselberg[, "logalpha"]

# calculate the slope based on the overall slope
lm2b_slopes <- lm2b@beta[2] + lm2b_slopes

# pull into a data.frame
lm2b_slopes <- dplyr::tibble(Inselberg = row.names(lm2b_ranef$Inselberg),
                             logalpha_slope = lm2b_slopes)

