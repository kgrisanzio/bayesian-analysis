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
       "sjPlot")


##---------Load and process data-------

hcp_data <- read.csv("hcp_data.csv")
morbid_data <- read.csv("morbid_data.csv")

dim(hcp_data)
head(hcp_data)
str(hcp_data)

dim(morbid_data)
head(morbid_data)
str(morbid_data)

# Convert relevant vars to factors
hcp_data$PIN <- as.factor(hcp_data$PIN)
hcp_data$gender <- as.factor(hcp_data$gender)
levels(hcp_data$gender) # 1: male, 2: female
morbid_data$SID <- as.factor(morbid_data$SID)
morbid_data$Choice <- as.factor(morbid_data$Choice)


##---------Quick T-test Frequentist vs. Bayesian-------

# Q: Is there a difference between males and females on internalizing symptoms?
# IV: "gender" (male: 1, female: 2)
# DV: "internal_scale_sch" (internalizing symptoms)

# Histogram to visualize DV
hist(hcp_data$internal_scale_sch, xlab = "Internalizing", main = "Histogram Internalizing")
# Right skewed DV

# Violin plots to visualize DV by gender
ggstatsplot::ggbetweenstats(data = hcp_data, x = gender, y = internal_scale_sch, messages = FALSE) 

# Frequentist t-test fit
t_freq <- t.test(internal_scale_sch ~ gender, data = hcp_data)
t_freq # not significant 

# Bayes factors version
t_BF <- ttestBF(formula = internal_scale_sch ~ gender, data = hcp_data)
t_BF
extractBF(t_BF)[,"bf"]
# The Bayes Factor is ~.24, so according to Jeffrey's (1961) conventions, there is "moderate
# evidence for H0"

# Extract the posterior
set.seed(123)
chains <- posterior(t_BF, iterations = 10000) # resample
head(chains)
summary(chains) # posterior mean, sd, and credible intervals
plot(chains[,2], trace = FALSE, main = "Posterior Mean Differences") # posterior distribution


##---------Full Bayesian Regression Workflow: sadness ~ age using rstanarm-------

# Q: Does age relate to sadness? 
# IV: "age"
# DV: "ML3" (continuous factor score generated from EFA)

# Histogram of DV
ggplot(hcp_data, aes(x = ML3)) + geom_histogram(color="black", fill="white", bins = 20)
# Our DV is continuous, right-skewed, and not bounded between two values.
# Candidate distributions:
  # Gamma distribution, but X>0 for gamma. We could add a constant to all values so all values are positive?
        # Is that allowed??
  # Gaussian as comparison

# Make new DV variable "ML3_pos" with values shifted to be all positive
hcp_data$ML3_pos <- hcp_data$ML3 + round(abs(min(hcp_data$ML3)),2)
range(hcp_data$ML3_pos) # all positive values now
ggplot(hcp_data, aes(x = ML3_pos)) + geom_histogram(color="black", fill="white", bins = 20)

# Scatterplot of IV/DV
ggplot(hcp_data, aes(x = age, y = ML3)) +
  geom_point() + 
  ggtitle("Age/Sadness Scatterplot")


# Frequentist fit first (for comparison)
fitlm <- glm(ML3 ~ age, data = hcp_data)
summary(fitlm)

fitlm_gamma <- glm(ML3_pos ~ age, data = hcp_data, family = Gamma)
summary(fitlm_gamma)
# We're just running these to compare to later, but wow does the relationship
# between our IV and DV really change from significantly positive to significantly
# negative when we use a normal vs. gamma distribution?!


# (Step 1) Specify the data model and priors

set.seed(123) 
fitBlm <- stan_glm(ML3 ~ age, data = hcp_data, refresh = 0) # gaussian() as default family setting
fitBlm_gamma <- stan_glm(ML3_pos ~ age, data = hcp_data, family = Gamma(link="inverse"), refresh = 0) # gamma
  # How do we know what link to use?
  # If we compare a gaussian to another distribution that needs transformed data (like gamma or beta), do we use the transformed
  # variable for the gaussian distribution (even though it doesn't need it), or the original variable?
prior_summary(fitBlm) 
prior_summary(fitBlm_gamma) # default priors in rstanarm (weakly informative)


# Posterior predictive checks
set.seed(123)
pp_check(fitBlm, nreps = 100) # PPC for Gaussian fit
pp_check(fitBlm_gamma, nreps = 100) # PPC for Gamma fit
# Gamma fits better than Gaussian


# (Step 2) Estimate the model parameters
print(fitBlm_gamma, digits = 3)


# (Step 3) Check sampling quality
summary(fitBlm_gamma, digits = 3) 
# Checks look good - Rhat (scale reduction factor) is smaller than 1.1 and 
# n_eff (number of effective sample size) is greater than 1000

plot(fitBlm_gamma, "trace", pars = "(Intercept)")
plot(fitBlm_gamma, "trace", pars = "age")     
plot(fitBlm_gamma, "trace", pars = "shape")     
# We don't see any pattern in the trace plots


# (Step 4) Summarize and interpret the results

# Extract and plot the posterior draws
postIntercept <- as.matrix(fitBlm_gamma, pars = "(Intercept)")     
postSlope <- as.matrix(fitBlm_gamma, pars = "age")
postShape <- as.matrix(fitBlm_gamma, pars = "shape")

# Histograms
op <- par(mfrow = c(3,1))
hist(postIntercept, freq = FALSE, main = "Posterior Intercept", xlab = bquote(beta[0]))
hist(postSlope, freq = FALSE, main = "Posterior Slope", xlab = bquote(beta[1]))
hist(postShape, freq = FALSE, main = "Posterior Shape", xlab = bquote(shape))
par(op)

# Kernel densities
op <- par(mfrow = c(3,1))                  
plot(density(postIntercept), main = "Posterior Intercept", xlab = bquote(beta[0]))
plot(density(postSlope), main = "Posterior Slope", xlab = bquote(beta[1]))
plot(density(postShape), main = "Posterior Shape", xlab = bquote(shape))
par(op)

summary(fitBlm_gamma, pars = "age", digits = 3) 
# We see here that 0 is not contained in the the 80% credible interval

# Plotting with 95% UI
plot(fitBlm_gamma, plotfun = "areas", pars = "(Intercept)", prob = 0.95)
plot(fitBlm_gamma, plotfun = "areas", pars = "age", prob = 0.95)
plot(fitBlm_gamma, plotfun = "areas", pars = "shape", prob = 0.95)

# Compare with lm() fit
describe_posterior(fitBlm_gamma, digits = 3)
confint(fitlm_gamma)
summary(fitlm_gamma)
# Very similar


# Plot regression lines
posterior1 <- as.matrix(fitBlm_gamma)
means1 <- colMeans(posterior1)
betas <- posterior1[sample(nrow(posterior1), 100), 1:2]       ## sample of 100 beta_0/beta_1 pairs
blues <- color_scheme_get("brightblue")
mod1p1 <- ggplot(hcp_data, aes(x = age, y = ML3_pos)) +
  geom_point(color = "gray30") +
  geom_abline(intercept = betas[, 1], slope = betas[, 2], color = blues[[2]], size = 0.15, alpha = 0.5) +
  geom_abline(intercept = means1[1], slope = means1[2], size = 1.25, color = blues[[6]]) 
plot(mod1p1)

# Effect plot
plot_model(fitBlm_gamma, type = "pred", terms = c("age"))
predict(fitBlm_gamma)
# Predicting with new observations
predict(fitBlm_gamma, newdata = data.frame(age = c(10, 14, 17))) 
# So we can expect someone aged 10, say, to have a sadness score of ~0.83


##---------Bayesian Regression: depression symptoms ~ age + gender using brms-------

op <- options(mc.cores = 4) 
options()$mc.cores            

# Q: Does age relate to depression symptoms? 
# IV: "age", "gender"
# DV: "dsm5_scale1" (count of depression symptoms)

# Histogram of DV
hist(hcp_data$dsm5_scale1, breaks = 10, xlab = "Depression counts", main = "Histogram Depression")

# Candidate distributions: 
  # Poisson, 0-inflated Poisson, negative-binomial, 0-inflated negative-binomial 
  # Gaussian just for comparison

# Model fit and priors
set.seed(123)
fit_gaus <- brm(dsm5_scale1 ~ age + gender, data = hcp_data) # Gaussian
prior_summary(fit_gaus)

fit_pois <- brm(dsm5_scale1 ~ age + gender, family = poisson, data = hcp_data) # poisson
prior_summary(fit_pois)

fit_0pois <- brm(dsm5_scale1 ~ age + gender, family = zero_inflated_poisson, data = hcp_data) # 0-inflated poisson
prior_summary(fit_0pois)

fit_nb <- brm(dsm5_scale1 ~ age + gender, family = negbinomial, data = hcp_data) # negative-binomial
prior_summary(fit_nb)

fit_0nb <- brm(dsm5_scale1 ~ age + gender, family = zero_inflated_negbinomial, data = hcp_data) # 0-inflated negative-binomial
prior_summary(fit_0nb)

if (!file.exists("fit_gaus.rda")) {
  fit_gaus <- brm(dsm5_scale1 ~ age + gender, data = hcp_data)
  save(fit_gaus, file = "fit_gaus.rda")
} else {
  load("fit_gaus.rda")
}


if (!file.exists("fit_pois.rda")) {
  fit_pois <- brm(dsm5_scale1 ~ age + gender, family = poisson, data = hcp_data)
  save(fit_pois, file = "fit_pois.rda")
} else {
  load("fit_pois.rda")
}

if (!file.exists("fit_0pois.rda")) {
  fit_0pois <- brm(dsm5_scale1 ~ age + gender, family = zero_inflated_poisson, data = hcp_data)
  save(fit_0pois, file = "fit_0pois.rda")
} else {
  load("fit_0pois.rda")
}

if (!file.exists("fit_nb.rda")) {
  fit_nb <- brm(dsm5_scale1 ~ age + gender, family = negbinomial, data = hcp_data)
  save(fit_nb, file = "fit_nb.rda")
} else {
  load("fit_nb.rda")
}

if (!file.exists("fit_0nb.rda")) {
  fit_0nb <- brm(dsm5_scale1 ~ age + gender, family = zero_inflated_negbinomial, data = hcp_data)
  save(fit_0nb, file = "fit_0nb.rda")
} else {
  load("fit_0nb.rda")
}


# Which distribution to proceed with?

set.seed(5)
brms::pp_check(fit_gaus, ndraws = 100) # clearly doesn't fit well
brms::pp_check(fit_gaus, ndraws = 11, type = "hist")

brms::pp_check(fit_pois, ndraws = 100) 
brms::pp_check(fit_pois, ndraws = 11, type = "hist")

brms::pp_check(fit_0pois, ndraws = 100) 
brms::pp_check(fit_0pois, ndraws = 11, type = "hist")

brms::pp_check(fit_nb, ndraws = 100) 
brms::pp_check(fit_nb, ndraws = 11, type = "hist")

brms::pp_check(fit_0nb, ndraws = 100) 
brms::pp_check(fit_0nb, ndraws = 11, type = "hist")


# R squared
bayes_R2(fit_gaus)
bayes_R2(fit_pois)
bayes_R2(fit_0pois)
bayes_R2(fit_nb)
bayes_R2(fit_0nb)

# Based on plots and R2 we can go with.. negative binomial?

# Check sampling quality
fit_nb # Rhat and ESS are within recommended ranges
plot(fit_nb, N = 5) # no pattern to trace plots


# Summarize and interpret results
describe_posterior(fit_nb)

conditional_effects(fit_nb, effects = "age")
conditional_effects(fit_nb, effects = "gender")

# For both of our IVs, the 95% credible interval contains 0, so we can say neither has a "significant" effect
# on our depression symptoms variable.

# Let's compare this really quick to a frequentist approach to double check
fitglm_dep <- glm.nb(dsm5_scale1 ~ age + gender, data = hcp_data)
summary(fitglm_dep)
confint(fitglm_dep)
# Similar and also not significant


##---------Bayesian Hierarchical Models-------

# Q: Does age relate to depression symptoms? 
# IV: "age", "gender"
# DV: "dsm5_scale1" (count of depression symptoms)

# Frequentist fit
M_viewtime <- lmer(mean_viewingRT ~ Choice + (1|SID), data = morbid_data); summary(M_viewtime)

# Visualize

# Histogram of DV
hist(morbid_data$mean_viewingRT, breaks = 10, xlab = "Viewing Time", main = "Histogram Viewing Time")
# Skewed reaction time 

ggline(morbid_data, x = "Choice", y = "mean_viewingRT", 
       add = c("mean_se", "jitter"), 
       ylab = "Viewing Time (Seconds)", xlab = "Choice")
# Subjects looked longer at negative images than at neutral and positive images

tspag <- ggplot(morbid_data, aes(x = Choice, y = mean_viewingRT, group = SID)) + xlab("Choice") + geom_line() +
  ylab("Mean Viewing Time") + ggtitle("Age by ")
tspag + aes(colour = SID) +  guides(colour=FALSE) 
# Difficult to tell pattern here with lots of data


# Bayesian fit

set.seed(123)
M_hier <- stan_lmer(mean_viewingRT ~ Choice + (1|SID), iter = 5000, data = morbid_data)

if (!file.exists("M_hier.rda")) {
  M_hier <- stan_lmer(mean_viewingRT ~ Choice + (1|SID), iter = 5000, data = morbid_data)
  save(M_hier, file = "M_hier.rda")
} else {
  load("M_hier.rda")
}

describe_posterior(M_hier, effects = "fixed")
describe_posterior(M_hier, effects = "random")
# Increased iterations due to warning

prior_summary(M_hier)

# Plot
plot_model(M_hier, type = "re", ci.lvl = 0.95, title = "Hierarchical Model") 


# Bayesian fit - gamma

set.seed(123)
M_hier_gamma <- stan_glmer(mean_viewingRT ~ Choice + (1|SID), family = Gamma(link = "inverse"), iter = 3000, data = morbid_data)

if (!file.exists("M_hier_gamma.rda")) {
  M_hier_gamma <- stan_glmer(mean_viewingRT ~ Choice + (1|SID), family = Gamma(link = "inverse"), iter = 3000, data = morbid_data)
  save(M_hier_gamma, file = "M_hier_gamma.rda")
} else {
  load("M_hier_gamma.rda")
}

describe_posterior(M_hier_gamma, effects = "fixed")
describe_posterior(M_hier_gamma, effects = "random")
# Increased iterations due to warning

prior_summary(M_hier_gamma)

# Plot
plot_model(M_hier_gamma, type = "re", ci.lvl = 0.95, title = "Hierarchical Model") 


