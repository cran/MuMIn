# Extends: survival    #########################################################

# logLik for survival::coxph model
# https://stat.ethz.ch/pipermail/r-help/2006-December/122118.html
# originally by Charles C. Berry, mod. by KB: correction for the null model
`logLik.coxph` <- function(object,...) {
# Thx to Mathieu Basille:
    ret <-  object$loglik[length(object$loglik)]
	#ret <-  if(length(object$loglik) > 1)
	#	-1 * (object$loglik[1] - object$loglik[2]) else object$loglik
	#attr(ret,"nall") <-
	attr(ret, "nobs") <- object$n
    attr(ret, "df") <- if(is.null(object$coef)) 0 else sum(!is.na(object$coef))
	class(ret) <- "logLik"
	return(ret)
}

`logLik.survreg` <- function (object, ...)  {
	ret <- object$loglik[2L]
	attr(ret, "nobs") <- length(object$linear)
	attr(ret, "df") <- sum(object$df)
	class(ret) <- "logLik"
	return(ret)
}

`logLik.glmmML` <- function(object, ...) {
	ret <- -object$deviance / 2
	#ret <- df  - object$aic / 2
	n <- length(object$coefficients)
	attr(ret, "df") <- n + object$cluster.null.df - object$df.residual
	attr(ret, "nobs") <- n + object$cluster.null.df
	class(ret) <- "logLik"
	return(ret)
}

`logLik.glmmboot` <- function (object, ...) {
	ret <- object$logLik
	attr(ret, "nobs") <- object$n
	attr(ret, "df") <- object$n - object$df.residual
	class(ret) <- "logLik"
	return(ret)
}

