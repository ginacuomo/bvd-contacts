library(tidyverse)

# input
log_odds_baseline <- 0
beta_a <- 1
beta_b <- -1

# logit(p) = log(p / (1 - p)), probabilties to log odds
logit <- function(...) stats::qlogis(...)
# logistic(x) = 1 / (1 + exp(-x)), log odds to probabilties
logistic <- function(...) stats::plogis(...)

# Simulate data under a logistic regression model; make the covariates
# correlated so their associated regression coefficients are hard to estimate
# separately.
df <- tibble(a = c(rep(1, 40), rep(0, 40), rbernoulli(20)) %>% as.logical,
             b = c(rep(1, 40), rep(0, 40), rbernoulli(20)) %>% as.logical,
             log_odds = log_odds_baseline + if_else(a, beta_a, 0) + if_else(b, beta_b, 0),
             p = logistic(log_odds),
             outcome = rbernoulli(100, p))

# Fit the model, compare the estimated regression coefficients to the true values
model <- glm(data = df, outcome ~ a + b, family = "binomial")
summary(model)

# Report the estimated probability for a manually specified combination of a and
# b variables, using a standard error appropriate to the joint uncertainty in
# beta_a and beta_b. (We understand their combined effect well, though not their
# separate effects.)
predict.glm(model,
            newdata = tibble(a = TRUE,
                             b = TRUE),
            type = "response", # predictions on the response (outcome) scale, not the predictors scale
            se.fit = TRUE)
