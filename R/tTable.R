`tTable` <-
function (model, ...) 	UseMethod("tTable")

`tTable.default` <-
function(model, ...) return(summary(model)$coefficients)

`tTable.gam` <-
function(model, ...) {
	cf <- model$coefficients
	se <- summary(model)$se
     return(cbind(`Estimate`=cf, `Std. Error` = se))
}

`tTable.gls` <-
function (model, ...) return(summary(model)$tTable)

`tTable.lme` <-
function(model, ...) return(summary(model)$tTable)


# these are for old (buggy) versions of lme4
`tTable.mer` <-
`tTable.glmer` <-
`tTable.lmer` <-
function(model, ...) {
	sm <- eval(expression(summary), environment(lmer))
	return (sm(model)@coefs)
	#return((lme4::summary(model))@coefs)
}

`tTable.sarlm` <-
`tTable.spautolm` <-
function(model, ...) return(summary(model)$Coef)


`tTable.multinom` <-
function(model, ...) {
	ret <- as.data.frame(summary(model)[c("coefficients", "standard.errors")])
	colnames(ret) <- c("Estimate", "Std. Error")
	return(ret)
}

`tTable.coxph` <- function (model, ...)
	return(summary(model)$coefficients[,-2, drop=FALSE])
