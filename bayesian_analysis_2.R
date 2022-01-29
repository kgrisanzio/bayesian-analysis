# Katherine Grisanzio
# Bayesian statistical analysis

##---------Load packages-------

if (!require(pacman)) {install.packages("pacman"); require(pacman)}
p_load("ggstatsplot",
       "BayesFactor",
       "brms",
       "bayestestR",
       "MASS", 
       "ggplot2",
       "rstanarm",
       "sjPlot",
       "devtools")
if (!require("INLA")) {install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE); require("INLA")}  
if (!require("INLAutils")) {install_github("timcdlucas/INLAutils"); require("INLAutils")}
if (!require("brinla")) {install_github("julianfaraway/brinla"); require("brinla")}
if (!require("coefINLA")) {install_github("hesscl/coefINLA"); require("coefINLA")} 

##---------Load and process data-------

hcp_data <- read.csv("hcp_data.csv")

dim(hcp_data)
head(hcp_data)
str(hcp_data)

# Convert relevant vars to factors
hcp_data$PIN <- as.factor(hcp_data$PIN)
hcp_data$gender <- as.factor(hcp_data$gender)
levels(hcp_data$gender) # 1: male, 2: female

# Rename four variables for ease of interpretation
names(hcp_data)[names(hcp_data) == "ML1"] <- "general_anxiety"
names(hcp_data)[names(hcp_data) == "ML2"] <- "anger"
names(hcp_data)[names(hcp_data) == "ML4"] <- "evaluative_anxiety"
names(hcp_data)[names(hcp_data) == "ML3"] <- "sadness"


##---------GAMs with INLA-------

# DV: sadness, anger, evaluative anxiety, general anxiety
# IV: age

# Examine scatter plots of age x each DV
ggplot(hcp_data, aes(x = age, y = general_anxiety)) + geom_point() + ggtitle("General Anxiety Scatterplot")
ggplot(hcp_data, aes(x = age, y = anger)) + geom_point() + ggtitle("Anger Scatterplot")
ggplot(hcp_data, aes(x = age, y = evaluative_anxiety)) + geom_point() + ggtitle("Evaluative Anxiety Scatterplot")
ggplot(hcp_data, aes(x = age, y = sadness)) + geom_point() + ggtitle("Sadness Scatterplot")
# Most of these have no obvious shape examining by eye, except maybe evaluative anxiety which appears to increase ~linearly

# Model 1: age as the only predictor
fit_genanx1 <- inla(general_anxiety ~ f(age, model = "rw1"), data = hcp_data, control.predictor = list(compute = TRUE))

# We get an error running this model originally:
# Error in inla.check.location(location[[r]], term = gp$random.spec[[r]]$term,  : 
#                             Locations are too close for f(age, model="rw1", ...):  min(diff(sort(x)))/diff(range(x)) = 7.944e-06 < 1e-05
#                             You can fix this by some kind of binning, see ?inla.group
#                             If you want/need to bypass this check at your own risk, do
#                             > m = get("inla.models", inla.get.inlaEnv())
#                             > m$latent$rw1$min.diff = NULL
#                             > assign("inla.models", m, inla.get.inlaEnv())

# First let's gry the inla.group, as suggested:
fit_genanx1 <- inla(general_anxiety ~ f(inla.group(age), model = "rw1"), data = hcp_data, control.predictor = list(compute = TRUE))
summary(fit_genanx1)
bri.band.ggplot(fit_genanx1, name = "age", type = "random", xlab = "age") + ggtitle("General Anxiety Trajectory")    

# It runs now, but bri.band.ggplot now gives an error:
# Error in FUN(X[[i]], ...) : object 'f.lb' not found

# Let's try their other suggestion, which is to run these lines to bypass the problem at our own risk:
m = get("inla.models", inla.get.inlaEnv())
m$latent$rw1$min.diff = NULL
assign("inla.models", m, inla.get.inlaEnv())

# Re-run  model
fit_genanx1 <- inla(general_anxiety ~ f(age, model = "rw1"), data = hcp_data, control.predictor = list(compute = TRUE))
summary(fit_genanx1)
bri.band.ggplot(fit_genanx1, name = "age", type = "random", xlab = "age") + ggtitle("General Anxiety Trajectory")    

# This now works, but since we ignored a warning we probalby need to examine output closely.
# Back to our models - 

# Model 1: age as the only predictor and model = "rw1" (random walk order of 1)
fit_genanx1 <- inla(general_anxiety ~ f(age, model = "rw1"), data = hcp_data, control.predictor = list(compute = TRUE))
fit_ang1 <- inla(anger ~ f(age, model = "rw1"), data = hcp_data, control.predictor = list(compute = TRUE))
fit_evalanx1 <- inla(evaluative_anxiety ~ f(age, model = "rw1"), data = hcp_data, control.predictor = list(compute = TRUE))
fit_sad1 <- inla(sadness ~ f(age, model = "rw1"), data = hcp_data, control.predictor = list(compute = TRUE))

summary(fit_genanx1)
summary(fit_ang1)
summary(fit_evalanx1)
summary(fit_sad1)

bri.band.ggplot(fit_genanx1, name = "age", type = "random", xlab = "age") + ggtitle("General Anxiety Trajectory\n(Random Walk Order: 1)")    
bri.band.ggplot(fit_ang1, name = "age", type = "random", xlab = "age") + ggtitle("Anger Trajectory\n(Random Walk Order: 1)")    
bri.band.ggplot(fit_evalanx1, name = "age", type = "random", xlab = "age") + ggtitle("Evaluative Anxiety Trajectory\n(Random Walk Order: 1)")    
bri.band.ggplot(fit_sad1, name = "age", type = "random", xlab = "age") + ggtitle("Sadness Trajectory\n(Random Walk Order: 1)")    

# These all look too wiggly. General anxiety and anger appear mostly flat across age, evaluative anxiety
# appears to increase linearly with age, and sadness appears to increase mostly in the second part of our
# age range. Let's try model = "rw2" (random walk order of 2) next.

# Code from warning to bypass error (now for "rw2")
m = get("inla.models", inla.get.inlaEnv())
m$latent$rw2$min.diff = NULL
assign("inla.models", m, inla.get.inlaEnv())

# Model 2: age as the only predictor and model = "rw2" (random walk order of 2)
fit_genanx2 <- inla(general_anxiety ~ f(age, model = "rw2"), data = hcp_data, control.predictor = list(compute = TRUE))
fit_ang2 <- inla(anger ~ f(age, model = "rw2"), data = hcp_data, control.predictor = list(compute = TRUE))
fit_evalanx2 <- inla(evaluative_anxiety ~ f(age, model = "rw2"), data = hcp_data, control.predictor = list(compute = TRUE))
fit_sad2 <- inla(sadness ~ f(age, model = "rw2"), data = hcp_data, control.predictor = list(compute = TRUE))

summary(fit_genanx2)
summary(fit_ang2)
summary(fit_evalanx2)
summary(fit_sad2)

bri.band.ggplot(fit_genanx2, name = "age", type = "random", xlab = "age") + ggtitle("General Anxiety Trajectory\n(Random Walk Order: 2)")    
bri.band.ggplot(fit_ang2, name = "age", type = "random", xlab = "age") + ggtitle("Anger Trajectory\n(Random Walk Order: 2)")    
bri.band.ggplot(fit_evalanx2, name = "age", type = "random", xlab = "age") + ggtitle("Evaluative Anxiety Trajectory\n(Random Walk Order: 2)")    
bri.band.ggplot(fit_sad2, name = "age", type = "random", xlab = "age") + ggtitle("Sadness Trajectory\n(Random Walk Order: 2")    

# These look better and less wiggly. General anxiety and anger still appear to be mostly flat across age,
# evaluative anxiety definitely appears to increase linearly with age, and sadness still appears to increase during
# the second part of our age range.

# Fitted values
yhat_genanx <- fit_genanx2$summary.fitted.values$mean 
yhat_ang <- fit_ang2$summary.fitted.values$mean 
yhat_evalanx <- fit_evalanx2$summary.fitted.values$mean 
yhat_sad <- fit_sad2$summary.fitted.values$mean 

# R-squared values
round(cor(yhat_genanx, hcp_data$general_anxiety)^2, 3)   
round(cor(yhat_ang, hcp_data$anger)^2, 3)    
round(cor(yhat_evalanx, hcp_data$evaluative_anxiety)^2, 3)    
round(cor(yhat_sad, hcp_data$sadness)^2, 3)    

# Plot regression lines
ind <- order(hcp_data$age)

with(hcp_data, plot(age, general_anxiety, xlab = "Age", ylab = "General Anxiety", main = "General Anxiety Regression"))
lines(sort(hcp_data$age), yhat_genanx[ind], lwd = 2, col = "red")

with(hcp_data, plot(age, anger, xlab = "Age", ylab = "Anger", main = "Anger Regression"))
lines(sort(hcp_data$age), yhat_ang[ind], lwd = 2, col = "red")

with(hcp_data, plot(age, evaluative_anxiety, xlab = "Age", ylab = "Evaluative Anxiety", main = "Evaluative Anxiety Regression"))
lines(sort(hcp_data$age), yhat_evalanx[ind], lwd = 2, col = "red")

with(hcp_data, plot(age, sadness, xlab = "Age", ylab = "Sadness", main = "Sadness Regression"))
lines(sort(hcp_data$age), yhat_sad[ind], lwd = 2, col = "red")


# Model 3: age and gender as predictors and model = "rw2" (random walk order of 2)
fit_genanx3 <- inla(general_anxiety ~ gender + f(age, model = "rw2"), data = hcp_data, control.predictor = list(compute = TRUE))
fit_ang3 <- inla(anger ~ gender + f(age, model = "rw2"), data = hcp_data, control.predictor = list(compute = TRUE))
fit_evalanx3 <- inla(evaluative_anxiety ~ gender + f(age, model = "rw2"), data = hcp_data, control.predictor = list(compute = TRUE))
fit_sad3 <- inla(sadness ~ gender + f(age, model = "rw2"), data = hcp_data, control.predictor = list(compute = TRUE))

summary(fit_genanx3) # gender has an influence
summary(fit_ang3) # gender has no influence
summary(fit_evalanx3) # gender has an influence
summary(fit_sad3) # gender has an influence


##---------GAMs with BRMS - default priors-------

# Default priors

op <- options(mc.cores = 4) 
options()$mc.cores     

if (!file.exists("fit_genanx4.rda")) {
  fit_genanx4 <- brm(bf(general_anxiety ~ s(age)), data = hcp_data, family = gaussian(), cores = 4, seed = 42,
                     iter = 5000, warmup = 1000, thin = 10)
  save(fit_genanx4, file = "fit_genanx4.rda")
} else {
  load("fit_genanx4.rda")
}

if (!file.exists("fit_ang4.rda")) {
  fit_ang4 <- brm(bf(anger ~ s(age)), data = hcp_data, family = gaussian(), cores = 4, seed = 42,
                     iter = 5000, warmup = 1000, thin = 10)
  save(fit_ang4, file = "fit_ang4.rda")
} else {
  load("fit_ang4.rda")
}

if (!file.exists("fit_evalanx4.rda")) {
  fit_evalanx4 <- brm(bf(evaluative_anxiety ~ s(age)), data = hcp_data, family = gaussian(), cores = 4, seed = 42,
                  iter = 5000, warmup = 1000, thin = 10)
  save(fit_evalanx4, file = "fit_evalanx4.rda")
} else {
  load("fit_evalanx4.rda")
}

if (!file.exists("fit_sad4.rda")) {
  fit_sad4 <- brm(bf(sadness ~ s(age)), data = hcp_data, family = gaussian(), cores = 4, seed = 42,
                      iter = 5000, warmup = 1000, thin = 10)
  save(fit_sad4, file = "fit_sad4.rda")
} else {
  load("fit_sad4.rda")
}

summary(fit_genanx4)
summary(fit_ang4)
summary(fit_evalanx4)
summary(fit_sad4)
# sds(sage_1): variance parameter, which controls the wiggliness of the smooth. The larger the value, the more wiggly.
  # if credible interval doesn't include 0, it suggests a smooth is needed
# sage_1: fixed effect of spline

# Plots
plot(fit_genanx4)
plot(conditional_effects(fit_genanx4), points = TRUE)

plot(fit_ang4)
plot(conditional_effects(fit_ang4), points = TRUE)

plot(fit_evalanx4)
plot(conditional_effects(fit_evalanx4), points = TRUE)

plot(fit_sad4)
plot(conditional_effects(fit_sad4), points = TRUE)

# Posterior predictive checks
pp_check(fit_genanx4)
pp_check(fit_ang4)
pp_check(fit_evalanx4)
pp_check(fit_sad4)

# Let's take just sadness (which we expect to be non-linear) and evaluative anxiety (which we expect to be linear)
# and compare them to a linear fit using brms:

if (!file.exists("sadness_lin.rda")) {
  sadness_lin <- brm(sadness ~ age, data = hcp_data)
  save(sadness_lin, file = "sadness_lin.rda")
} else {
  load("sadness_lin.rda")
}

if (!file.exists("evalanx_lin.rda")) {
  evalanx_lin <- brm(evaluative_anxiety ~ age, data = hcp_data)
  save(evalanx_lin, file = "evalanx_lin.rda")
} else {
  load("evalanx_lin.rda")
}

summary(sadness_lin)
summary(evalanx_lin)

# Compare linear and spline model posterior predictive checks
pp_check(fit_sad4)
pp_check(sadness_lin) # these look the same

pp_check(fit_evalanx4)
pp_check(evalanx_lin) # these also look the same

# Compare linear and GAM model using leave one out cross validation
loo(sadness_lin, fit_sad4)
loo(evalanx_lin, fit_evalanx4)
# smaller LOOIC values indicate better model fit, so we see that a GAM fits better
# than the linear model for sadness, and the linear model fits better than the GAM
# for evaluative anxiety. This matches what we expect given the INLA results.


##---------GAMs with BRMS - setting priors-------

# Now let's change the priors
# We'll Start with sadness as the DV

get_prior(bf(sadness ~ s(age)), data = hcp_data, family = gaussian(), cores = 4, seed = 42,
    iter = 5000, warmup = 1000, thin = 10)


# Use N(0, 100) slope priors for age, keep intercept and sigma as defaults
priors1 <- c(prior(normal(0, 100), class = b, coef = "sage_1"))
priors1 # new priors

# Fit model again with new priors
if (!file.exists("fit_sad5.rda")) {
  fit_sad5 <- brm(bf(sadness ~ s(age)), data = hcp_data, family = gaussian(), cores = 4, prior = priors1, seed = 42,
                  iter = 5000, warmup = 1000, thin = 10)
  save(fit_sad5, file = "fit_sad5.rda")
} else {
  load("fit_sad5.rda")
}
prior_summary(fit_sad5) # double check new priors

# Now use N(0, 20) slope priors
priors2 <- c(prior(normal(0, 20), class = b, coef = "sage_1"))
priors2 # new priors

if (!file.exists("fit_sad6.rda")) {
  fit_sad6 <- brm(bf(sadness ~ s(age)), data = hcp_data, family = gaussian(), cores = 4, prior = priors2, seed = 42,
                  iter = 5000, warmup = 1000, thin = 10)
  save(fit_sad6, file = "fit_sad6.rda")
} else {
  load("fit_sad6.rda")
}
prior_summary(fit_sad6)

# And N(0, 10) slope priors
priors3 <- c(prior(normal(0, 10), class = b, coef = "sage_1"))
priors3 # new priors

if (!file.exists("fit_sad7.rda")) {
  fit_sad7 <- brm(bf(sadness ~ s(age)), data = hcp_data, family = gaussian(), cores = 4, prior = priors3, seed = 42,
                  iter = 5000, warmup = 1000, thin = 10)
  save(fit_sad7, file = "fit_sad7.rda")
} else {
  load("fit_sad7.rda")
}
prior_summary(fit_sad7)

# Finally, N(0, 1) slope priors
priors4 <- c(prior(normal(0, 1), class = b, coef = "sage_1"))
priors4 # new priors

if (!file.exists("fit_sad8.rda")) {
  fit_sad8 <- brm(bf(sadness ~ s(age)), data = hcp_data, family = gaussian(), cores = 4, prior = priors4, seed = 42,
                  iter = 5000, warmup = 1000, thin = 10)
  save(fit_sad8, file = "fit_sad8.rda")
} else {
  load("fit_sad8.rda")
}
prior_summary(fit_sad8)

# Compare models
parameters::model_parameters(fit_sad5) # N(0, 100) prior
parameters::model_parameters(fit_sad6) # N(0, 20) prior
parameters::model_parameters(fit_sad7) # N(0, 10) prior
parameters::model_parameters(fit_sad8) # N(0, 1) prior
# We see the median and 95% CI stays roughly the same for the first three, then we see shrinkage with 
# the more narrow N(0, 1) prior. Let's keep the N(0, 10) prior.

# MCMC diagnostics and posterior plots
plot(fit_sad7, variable = "^b_", regex = TRUE)
# I think I can see some wiggle to the plot, but overall looks good

# R-squared
bayes_R2(fit_sad7)

# 95% CIs
parameters::model_parameters(fit_sad7)

# Effects plot
plot_model(fit_sad7) 

# ROPE
rope_range(fit_sad7) # default
plot(rope(fit_sad7, ci = 0.95))


# Let's try changing the priors on another model that might have a stronger effect
# We'll use evaluative anxiety as the DV and include gender as a covariate

get_prior(bf(evaluative_anxiety ~ gender + s(age)), data = hcp_data, family = gaussian(), cores = 4, seed = 42,
          iter = 5000, warmup = 1000, thin = 10)

# Set priors
priors11 <- c(prior(normal(0, 20), class = b, coef = "sage_1"),
              prior(normal(0, 20), class = b, coef = "gender2"))
priors12 <- c(prior(normal(0, 10), class = b, coef = "sage_1"),
              prior(normal(0, 10), class = b, coef = "gender2"))
priors13 <- c(prior(normal(0, 1), class = b, coef = "sage_1"),
              prior(normal(0, 1), class = b, coef = "gender2"))

if (!file.exists("fit_evalanx5.rda")) {
  fit_evalanx5 <- brm(bf(evaluative_anxiety ~ gender + s(age)), data = hcp_data, family = gaussian(), cores = 4, prior = priors11, seed = 42,
                      iter = 5000, warmup = 1000, thin = 10)
  save(fit_evalanx5, file = "fit_evalanx5.rda")
} else {
  load("fit_evalanx5.rda")
}
prior_summary(fit_evalanx5)

if (!file.exists("fit_evalanx6.rda")) {
  fit_evalanx6 <- brm(bf(evaluative_anxiety ~ gender + s(age)), data = hcp_data, family = gaussian(), cores = 4, prior = priors12, seed = 42,
                      iter = 5000, warmup = 1000, thin = 10)
  save(fit_evalanx6, file = "fit_evalanx6.rda")
} else {
  load("fit_evalanx6.rda")
}
prior_summary(fit_evalanx6)

if (!file.exists("fit_evalanx7.rda")) {
  fit_evalanx7 <- brm(bf(evaluative_anxiety ~ gender + s(age)), data = hcp_data, family = gaussian(), cores = 4, prior = priors13, seed = 42,
                      iter = 5000, warmup = 1000, thin = 10)
  save(fit_evalanx7, file = "fit_evalanx7.rda")
} else {
  load("fit_evalanx7.rda")
}
prior_summary(fit_evalanx7)

# Compare models
parameters::model_parameters(fit_evalanx5) # N(0, 20) prior
parameters::model_parameters(fit_evalanx6) # N(0, 10) prior
parameters::model_parameters(fit_evalanx7) # N(0, 1) prior
# With the N(0, 1) prior we start to see shrinkage for age (but not gender)
# Let's examine this model more closely

# Evaluate model

plot(fit_evalanx7, variable = "^b_", regex = TRUE) # intercept and slope posteriors
bayes_R2(fit_evalanx7) # r-squared
parameters::model_parameters(fit_evalanx7) # posterior output

# Effects plots 
plot_model(fit_evalanx7, type = "pred", terms = "age") # clear age effect   
plot_model(fit_evalanx7, type = "pred", terms = "gender") # clear gender effect (1 = male, 2 = female)

# ROPE
rope_range(fit_evalanx7) # default
plot(rope(fit_evalanx7), ci = 0.95)

