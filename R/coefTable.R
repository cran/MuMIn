`tTable` <-
function (model, ...) 	{
	.Deprecated("coefTable")
	UseMethod("coefTable")
}

`coefTable` <-
function (model, dispersion = NULL, ...) 	UseMethod("coefTable")

`coefTable.default` <-
function(model, dispersion = NULL, ...)
return(summary(model, dispersion = dispersion)$coefficients)

`coefTable.gam` <-
function(model, dispersion = NULL, ...) {
	cf <- model$coefficients
	se <- summary(model, dispersion = dispersion)$se
    return(cbind(`Estimate` = cf, `Std. Error` = se))
}

`coefTable.glmmML` <- function(model, ...) {
    coef <- model$coefficients
    se <- model$coef.sd
    ret <- cbind(coef, se, coef / se, signif(1 - pchisq((coef/se)^2,
        1)))
    dimnames(ret) <- list(names(coef), c("Estimate", "Std. Error",
        "z", "Pr(>|z|)"))
  return(ret)
}

`coefTable.gls` <-
function (model, ...) return(summary(model)$tTable)

`coefTable.lme` <-
function(model, ...) return(summary(model)$tTable)

`coefTable.mer` <-
function(model, ...) {
	#sm <- eval(expression(summary), as.environment("package:lme4"))
	sm <- eval(expression(summary), asNamespace("lme4"))
	return (sm(model)@coefs)
	#return((summary(model))@coefs)
}

`coefTable.multinom` <-
function(model, ...) {
	ret <- do.call("cbind", summary(model)[c("coefficients", "standard.errors")])
	colnames(ret) <- c("Estimate", "Std. Error")
	return(ret)
}

`coefTable.sarlm` <-
`coefTable.spautolm` <-
function(model, ...) {
	cf <- coef(model)
	se <- sqrt(diag(summary(model)$resvar))[names(cf)]
	return(cbind(`Estimate`=cf, `Std. Error` = se))
}

`coefTable.survreg` <- function (model, ...) {
	return(cbind(
		`Estimate` = model$coefficients,
		`Std. Error` = sqrt(diag(model$var)))
	)
}

`coefTable.coxph` <- function (model, ...)
	return(summary(model)$coefficients[,-c(2L, 3L), drop = FALSE])
