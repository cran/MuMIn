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

`tTable.glmmML` <- function(model, ...) {
    coef <- model$coefficients
    se <- model$coef.sd
    ret <- cbind(coef, se, coef/se, signif(1 - pchisq((coef/se)^2,
        1)))
    dimnames(ret) <- list(names(coef), c("Estimate", "Std. Error",
        "z", "Pr(>|z|)"))
  return(ret)
}

`tTable.gls` <-
function (model, ...) return(summary(model)$tTable)

`tTable.lme` <-
function(model, ...) return(summary(model)$tTable)


`tTable.mer` <-
function(model, ...) {
	#sm <- eval(expression(summary), as.environment("package:lme4"))
	sm <- eval(expression(summary), asNamespace("lme4"))
	return (sm(model)@coefs)
	#return((summary(model))@coefs)
}

`tTable.multinom` <-
function(model, ...) {
	ret <- do.call("cbind", summary(model)[c("coefficients", "standard.errors")])
	colnames(ret) <- c("Estimate", "Std. Error")
	return(ret)
}

`tTable.sarlm` <-
`tTable.spautolm` <-
function(model, ...) {
	cf <- coef(model)
	se <- sqrt(diag(summary(model)$resvar))[names(cf)]
	return(cbind(`Estimate`=cf, `Std. Error` = se))
}

`tTable.survreg` <- function (model, ...) {
	return(cbind(
		`Estimate` = model$coefficients,
		`Std. Error` = sqrt(diag(model$var)))
	)
}

`tTable.coxph` <- function (model, ...)
	return(summary(model)$coefficients[,-c(2L, 3L), drop = FALSE])