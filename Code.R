# estimate the gender wage gap on CPS data and explore heterogeneity

library(hdm)
library(glmnet)

# ------------ Double LASSO for estimation of gender wage gap

# (1) prepare data

# load data
data(cps2012)
str(cps2012)
head(cps2012)

# variable of interest
D <- cps2012$female	

# create predictors
# model.matrix converts a model formula to a actual matrix which we can then insert in the model params

# Build the high-dimensional control matrix W for Double LASSO:
# model.matrix() converts the formula into numeric regressors.
# -1 removes the intercept, because W should contain only covariates and not a constant term.
# (...)^2 includes all listed main effects plus all pairwise interactions,
# allowing flexible controls without estimating the final model yet.

# R formula syntax for interactions:
#
# x:z means: only the interaction between x and z
# Example:
#   x:z
# creates:
#   x:z
#
# x * z means: main effects plus interaction
# Example:
#   x * z
# is equivalent to:
#   x + z + x:z
#
# (a + b + c)^2 means:
# all main effects plus all pairwise interactions
# Example:
#   (a + b + c)^2
# is equivalent to:
#   a + b + c + a:b + a:c + b:c
#
# Important:
# In R formulas, ^2 does NOT mean "square each variable".
# It means "include interactions up to order 2".
#
# To literally square a variable, use I():
# Example:
#   I(exp1^2)
#
# In this code, exp2 and exp3 are already existing variables,
# so they are included as normal controls, not created by the formula.

W <- model.matrix(~-1 +(widowed + divorced + separated + nevermarried + hsd08 + hsd911 + hsg + cg + ad + mw + so +
                          we + exp1 + exp2 + exp3)^2, data = cps2012)


# Remove constant columns from W:
# apply(W, 2, var) computes the variance of each column.
# Columns with variance 0 contain no information, e.g. an interaction
# that is always 0 because two dummy categories are mutually exclusive.
# These columns can cause numerical problems and are not useful for LASSO,
# so we keep only columns with non-zero variance.
W <- W[, which(apply(W, 2, var) != 0)] # exclude all constant variables

# Demean / center each column of W:
# demean() subtracts the column mean from every observation.
# After this transformation, each control variable has mean 0.
#
# Example:
#   exp1_centered = exp1 - mean(exp1)
#
# This is useful for LASSO because the controls are centered around zero,
# the intercept is separated from the covariates, and numerical estimation
# is usually more stable.
demean <- function(x) { 
  x - mean(x) 
}

# apply(W, 2, FUN = demean) applies the demeaning function column-wise.
# The 2 means "apply over columns"; 1 would mean "apply over rows".
W <- apply(W, 2, FUN = demean)

# outcome
Y <- cps2012$lnw

# (2) double LASSO estimator

# ------------ Double LASSO / partialling-out estimator ------------
#
# Goal:
# Estimate the coefficient of D = female on Y = log wage,
# while flexibly controlling for many covariates W.
#
# The idea is based on partialling out:
#   1. Remove from Y the part explained by W.
#   2. Remove from D the part explained by W.
#   3. Regress the residualized Y on the residualized D.
#
# This targets the coefficient on D in a model like:
#   Y = alpha * D + g(W) + error
#
# Here alpha is the adjusted gender wage gap.

# Partial out W from Y:
# cv.glmnet(W, Y) runs LASSO regressions of Y on W for many lambda values.
# Since alpha = 1 by default, this is LASSO.
# The main purpose here is to use cross-validation to choose a good lambda
# penalty value; the default is 10-fold cross-validation.
fit.lasso.Y <- cv.glmnet(W, Y)

# Predict the part of Y that can be explained by W:
# predict(..., s = "lambda.min") uses the LASSO fit corresponding to the
# lambda value with the smallest cross-validation error.
# newx = W means we predict fitted values for the same observations.
# The output column is named "lambda.min", but the entries are predicted Y values,
# not lambda values
fitted.lasso.Y <- predict(fit.lasso.Y, newx = W, s = "lambda.min")

# Residualized outcome:
# Ytilde is the remaining part of log wage after removing the part
# predicted by the controls W.
Ytilde <- Y - fitted.lasso.Y


# Partial out W from D:
# cv.glmnet(W, D) runs LASSO regressions of D = female on W
# for many lambda values and chooses lambda by cross-validation.
# Since no family is specified, glmnet uses a Gaussian/linear model by default.
# Here this is used for residualizing D, not mainly for classification.
fit.lasso.D <- cv.glmnet(W, D)

# Predict the part of D that can be explained by W:
# Again, s = "lambda.min" uses the lambda with the smallest
# cross-validation error.
fitted.lasso.D <- predict(fit.lasso.D, newx = W, s = "lambda.min")

# Residualized treatment:
# Dtilde is the remaining part of female after removing the part
# predicted by the controls W.
Dtilde <- D - fitted.lasso.D


# Final partialling-out regression:
# Now we run the actual final regression for the Double LASSO estimator.
# We regress residualized wages on residualized gender.
# The coefficient on Dtilde estimates the adjusted gender wage gap.
fit.doubleLASSO <- lm(Ytilde ~ Dtilde)

# Extract the estimated coefficient on Dtilde.
beta1hat <- coef(fit.doubleLASSO)[2]

# Interpretation of the female coefficient:
#
# The outcome variable is log hourly wage:
#   Y = log(wage)
#
# Therefore, the coefficient on female is measured in log points.
# For example, an estimate of -0.28 means that women have about
# 0.28 lower log hourly wages than comparable men, conditional on controls.
#
# A rough interpretation is:
#   -0.28 ≈ -28%
#
# The exact percentage wage difference is obtained by exponentiating:
#   100 * (exp(beta) - 1)
#
# Example:
#   100 * (exp(-0.28) - 1) ≈ -24.4%
#
# So beta = -0.28 means that women earn about 24.4% less than comparable men,
# after adjusting for the selected/high-dimensional controls.
#
# Since this is observational data, this should be interpreted as an
# adjusted association unless strong causal assumptions are justified.


# (3) double LASSO estimator using hdm package

# ------------ Double LASSO estimator using the hdm package ------------

# Create one combined design matrix DW:
# DW contains both:
#   1. D = female, the variable of interest
#   2. W = all control variables and their pairwise interactions
#
# So conceptually:
#   DW = [D, W]
#
# This is different from the earlier W matrix, which contained only controls.
# Here female is included because rlassoEffects() expects the target variable
# and the controls together in one matrix.
DW <- model.matrix(
  ~ -1 + female + #this female col is the difference befor it was split to a different matrix
    (widowed + divorced + separated + nevermarried +
       hsd08 + hsd911 + hsg + cg + ad +
       mw + so + we +
       exp1 + exp2 + exp3)^2,
  data = cps2012
)

# Remove constant columns:
# Some interaction terms may be always 0 or otherwise constant.
# Such variables contain no useful information and can cause numerical issues.
DW <- DW[, which(apply(DW, 2, var) != 0)]

# Demean / center every column of DW:
# This subtracts the column mean from each variable.
# After this, every column has mean 0.
demean <- function(x) {
  x - mean(x)
}

DW <- apply(DW, 2, FUN = demean)

# Estimate the Double LASSO effect using hdm:
# rlassoEffects() estimates the effect of the variable specified by index,
# while treating all other columns as high-dimensional controls.
#
# index = 1 means:
#   estimate the effect of the first column of DW.
#
# Since DW was created with female first, index = 1 corresponds to female.
#
# Internally, rlassoEffects() does the partialling-out steps:
#   1. partial out W from Y
#   2. partial out W from D
#   3. regress residualized Y on residualized D
#
# Difference to the manual cv.glmnet version:
#   cv.glmnet() chooses lambda by cross-validation.
#   rlassoEffects() uses a theory-based/plugin penalty choice designed for inference.
#
# post = FALSE means:
#   use the LASSO estimates directly and do not refit post-LASSO OLS
#   on the selected variables.
fit.rlasso <- rlassoEffects(DW, Y, index = 1, post = FALSE)

# Print the estimated coefficient, standard error, t-statistic,
# p-value, and confidence interval for the effect of female on log wage.
summary(fit.rlasso)


# The estimate is very close to the manual Double LASSO estimate beta1hat,
# but not exactly identical because the two implementations choose the LASSO
# penalty differently:
#
#   manual version: cv.glmnet() chooses lambda by cross-validation
#   hdm version:    rlassoEffects() uses a theory-based/plugin penalty
#
# Therefore, the selected controls and residuals can differ slightly.
# Both estimates are around -0.28, meaning that female is associated with
# about 0.28 lower log hourly wages after adjusting for high-dimensional controls.
#
# Since the outcome is log wage, this corresponds approximately to a 28%
# lower wage, or more exactly:
#   100 * (exp(-0.28) - 1) ≈ -24.4%




# ------------ Recovering Heterogeneity in the Gender Wage Gap ------------

# The previous Double LASSO estimated one average adjusted gender wage gap:
#   effect of female on log wage.
#
# Here we ask a richer question:
#   Does the gender wage gap differ by marital status, education,
#   region, and experience?
#
# To study this, we include interactions between female and the characteristics.
# For example:
#   female:cg   = difference in the female wage gap for college graduates
#   female:ad   = difference in the female wage gap for advanced degrees
#   female:mw   = difference in the female wage gap in the Midwest
#   female:exp1 = how the female wage gap changes with experience
#
# The formula contains three parts:
#
#   1. female
#      Main female coefficient, i.e. the average/base gender wage gap.
#
#   2. female:(...)
#      Interactions between female and each listed characteristic.
#      In R formula syntax, ":" means interaction only.
#
#   3. (...)^2
#      Main effects and pairwise interactions among the control variables.
#      These are included as high-dimensional controls.
#
# So this setup allows the gender wage gap to be heterogeneous,
# while still flexibly controlling for many observed covariates.
X <- model.matrix(
  ~ -1 + female +
    female:(widowed + divorced + separated + nevermarried +
              hsd08 + hsd911 + hsg + cg + ad +
              mw + so + we +
              exp1 + exp2 + exp3) +
    (widowed + divorced + separated + nevermarried +
       hsd08 + hsd911 + hsg + cg + ad +
       mw + so + we +
       exp1 + exp2 + exp3)^2,
  data = cps2012
)

# Remove constant columns:
# Some interactions may be always 0 or otherwise constant.
# These contain no information and can cause numerical issues.
X <- X[, which(apply(X, 2, var) != 0)]

# Demean / center all columns:
# After this, each variable has mean 0.
demean <- function(x) {
  x - mean(x)
}

X <- apply(X, 2, FUN = demean)

# Select all female-related variables as target parameters:
# grep("female", colnames(X)) finds the column numbers of all variables
# whose names contain "female", e.g.:
#   female
#   female:widowed
#   female:cg
#   female:exp1
#
# These are the coefficients we want to estimate and interpret.
index.gender <- grep("female", colnames(X))

# Estimate the female-related coefficients using Double LASSO:
# rlassoEffects() estimates the selected target coefficients one by one,
# while using all other variables as high-dimensional controls.
#
# This is different from the previous step:
#   previous model: one average gender wage gap
#   this model: heterogeneous gender wage gaps by characteristics
fit <- rlassoEffects(X, Y, index = index.gender)

# Show coefficient estimates, standard errors, p-values, and confidence intervals.
summary(fit)

# Interpretation of significant heterogeneity effects:
#
# In the previous Double LASSO model, we estimated only one female effect:
#   Y = alpha * female + controls + error
#
# Therefore, the female coefficient from that model was interpreted as one
# average adjusted gender wage gap. In our case, this was around -0.28.
#
# In this heterogeneity model, we now include interactions between female
# and worker characteristics:
#   female + female:(marital status, education, region, experience)
#
# This means the female wage gap is allowed to differ across subgroups.
# The model is now of the form:
#   Y = theta0 * female
#       + theta1 * female:divorced
#       + theta2 * female:nevermarried
#       + theta3 * female:hsd911
#       + ...
#       + controls + error
#
# Therefore, the main female coefficient theta0 is NOT the overall average
# gender wage gap anymore. It is the baseline/reference female gap,
# i.e. the effect when the relevant female interaction terms are zero.
#
# A subgroup-specific gender wage gap is obtained by adding the relevant
# interaction coefficient(s):
#   subgroup gap = female coefficient + relevant female:subgroup interaction(s)
#
# So the earlier -0.28 and the new -0.1549 answer different questions:
#   -0.28   = one average adjusted gender wage gap from the simpler model
#   -0.1549 = baseline/reference female gap in the heterogeneity model
#
# Significant effects:
#
# female = -0.1549
# Baseline/reference female wage gap is negative.
# Women in the reference group earn about 0.155 lower log wages than comparable men.
# Exact percentage interpretation:
#   100 * (exp(-0.1549) - 1) ≈ -14.4%
#
# female:divorced = +0.1369
# The female wage gap is significantly less negative for divorced workers.
# Approximate subgroup gap:
#   -0.1549 + 0.1369 = -0.0180
# Exact percentage interpretation:
#   100 * (exp(-0.0180) - 1) ≈ -1.8%
#
# female:nevermarried = +0.1869
# The female wage gap is significantly less negative for never-married workers.
# Approximate subgroup gap:
#   -0.1549 + 0.1869 = 0.0319
# Exact percentage interpretation:
#   100 * (exp(0.0319) - 1) ≈ +3.2%
#
# female:hsd911 = -0.1193
# The female wage gap is significantly more negative for this education group.
# Approximate subgroup gap:
#   -0.1549 - 0.1193 = -0.2742
# Exact percentage interpretation:
#   100 * (exp(-0.2742) - 1) ≈ -24.0%
#
# female:exp2 = -0.1595 and female:exp3 = +0.0385
# The gender wage gap varies significantly and nonlinearly with experience.
# These polynomial experience interactions should be interpreted jointly,
# not separately as simple one-unit effects.
#
# Non-significant interaction terms are not interpreted as clear evidence
# of heterogeneity in the gender wage gap.
#
# Since the data are observational, these results should be interpreted as
# adjusted associations unless strong causal assumptions are justified.


# Plot the estimated female-related coefficients with 90% marginal confidence intervals.
# "Marginal" means intervals are shown coefficient by coefficient,
# not jointly adjusted for all coefficients at once.
plot(fit, level = 0.90)
