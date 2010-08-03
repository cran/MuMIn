`coeffs` <-
function (model) UseMethod("coeffs")

`coeffs.default` <-
function(model) model$coefficients

`coeffs.gls` <-
function (model) summary(model)$coefficients

`coeffs.lme` <-
function(model) model$coefficients$fixed

`coeffs.glmer` <-
`coeffs.lmer` <-
function(model) { 
	ret <- model@fixef
	names(ret) <- model@cnames$.fixed 
	return(ret)
}

`coeffs.mer` <-
function(model) model@fixef

`coeffs.spautolm` <-
function(model) model$fit$coefficients

