## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bbreg)

## -----------------------------------------------------------------------------
fit <- bbreg(agreement ~ priming + eliciting, data = WT)

## -----------------------------------------------------------------------------
fit

## -----------------------------------------------------------------------------
fit_beta <- bbreg(agreement ~ priming + eliciting, data = WT, model = "beta")
fit_beta

## -----------------------------------------------------------------------------
fit_priming <- bbreg(agreement ~ priming + eliciting | priming, data = WT)
fit_priming

## -----------------------------------------------------------------------------
fit_priming_eliciting <- bbreg(agreement ~ priming + eliciting | priming + eliciting, data = WT)
fit_priming_eliciting

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
summary(fit_beta)

## -----------------------------------------------------------------------------
summary(fit_priming)

## -----------------------------------------------------------------------------
fit_cloglog <- bbreg(agreement ~ priming + eliciting, data = WT, link.mean = "cloglog")
fit_cloglog

## -----------------------------------------------------------------------------
fit_cloglog_sqrt <- bbreg(agreement ~ priming + eliciting, data = WT, 
                          link.mean = "cloglog", link.precision = "sqrt")
fit_cloglog_sqrt

## -----------------------------------------------------------------------------
fit_priming_sqrt <- bbreg(agreement ~ priming + eliciting| priming, data = WT, 
                          link.precision = "sqrt")
fit_priming_sqrt

## -----------------------------------------------------------------------------
fitted_means <- fitted(fit)
head(fitted_means)

## -----------------------------------------------------------------------------
fitted_variances <- fitted(fit, type = "variance")
head(fitted_variances)

## -----------------------------------------------------------------------------
new_data_example <- data.frame(priming = c(0,0,1), eliciting = c(0,1,1))

## -----------------------------------------------------------------------------
predict(fit_priming_eliciting, newdata = new_data_example)
predict(fit_priming_eliciting, newdata = new_data_example, type = "variance")

## -----------------------------------------------------------------------------
predict(fit, newdata = new_data_example)
predict(fit, newdata = new_data_example, type = "variance")

## -----------------------------------------------------------------------------
predict_without_newdata <- predict(fit)
identical(predict_without_newdata, fitted_means)

## -----------------------------------------------------------------------------
fit_envelope <- bbreg(agreement ~ priming + eliciting, envelope = 300, data = WT)

## -----------------------------------------------------------------------------
summary(fit_envelope)

## -----------------------------------------------------------------------------
plot(fit)

## -----------------------------------------------------------------------------
plot(fit_envelope, which = 2)

## -----------------------------------------------------------------------------
plot(fit, which = c(1,4), ask = FALSE)

