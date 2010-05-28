`tTable.default` <-
function(model, ...) {
	return(summary(model)$coefficients)
}

`tTable.gam` <-
function(model, ...) {
	cf <- model$coefficients
	se <- summary(model)$se
     return(cbind(`Estimate`=cf, `Std. Error` = se))
}

`tTable.glmer` <-
function(model, ...) {
	sm <- eval(expression(summary), environment(lmer))
	return (sm(model)@coefs)
	#return((lme4::summary(model))@coefs)
}

`tTable.gls` <-
function (model, ...) return(summary(model)$tTable)

`tTable.lme` <-
function(model, ...) {
	return(summary(model)$tTable)
}

`tTable.lmer` <-
function(model, ...) {
	sm <- eval(expression(summary), environment(lmer))
	return (sm(model)@coefs)
	#return((lme4::summary(model))@coefs)
}

`tTable.mer` <-
function(model, ...) {
	sm <- eval(expression(summary), environment(lmer))
	return (sm(model)@coefs)
	#return((lme4::summary(model))@coefs)
}

`tTable` <-
function (model, ...) {
	UseMethod("tTable")
}

`tTable.spautolm` <-
function(model, ...) {
return(summary(model)$Coef)
}
