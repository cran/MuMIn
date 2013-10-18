# Extends: survival    #########################################################

# logLik for survival::coxph model
# https://stat.ethz.ch/pipermail/r-help/2006-December/122118.html
# originally by Charles C. Berry, mod. by KB: correction for the null model
`logLik.coxph` <- 
function(object,...) {
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

`logLik.survreg` <- 
function (object, ...)  {
	ret <- object$loglik[2L]
	attr(ret, "nobs") <- length(object$linear)
	attr(ret, "df") <- sum(object$df)
	class(ret) <- "logLik"
	return(ret)
}

`logLik.glmmML` <- 
function(object, ...) {
	ret <- -object$deviance / 2
	#ret <- df  - object$aic / 2
	n <- length(object$coefficients)
	attr(ret, "df") <- n + object$cluster.null.df - object$df.residual
	attr(ret, "nobs") <- n + object$cluster.null.df
	class(ret) <- "logLik"
	return(ret)
}

`logLik.glmmboot` <- 
function (object, ...) {
	ret <- object$logLik
	attr(ret, "nobs") <- object$n
	attr(ret, "df") <- object$n - object$df.residual
	class(ret) <- "logLik"
	return(ret)
}

`logLik.coxme` <-
function(object, type = c("penalized", "integrated"), ...) {
	type <- match.arg(type)
	i <- which(type == c("integrated", "penalized"))[1L]
	ret <- object$loglik[[i + 1L]]
	attr(ret, "df") <- object$df[i]
	attr(ret, "nobs") <- object$n[2L] # XXX: 1 or 2 ?
	class(ret) <- "logLik"
	ret
}

`logLik.lmekin` <-
function(object, ...) {
	ret <- object$loglik
	attr(ret, "nobs") <- object$n
	attr(ret, "df") <- length(object$coefficients$fixed) +
		length(object$coefficients$random) + 1L
	class(ret) <- "logLik"
	ret
}

`logLik.unmarkedFit` <- 
function(object, ...) {
	ret <- -object@negLogLike
	attr(ret, "df") <- length(object@opt$par)
	attr(ret, "nobs") <- #get("sampleSize", asNamespace("unmarked"))(object)
			unmarked::sampleSize(object)
	class(ret) <- "logLik"
	return(ret)
}

`logLik.splm` <- 
function (object, ...) {
	ret <- object$logLik
	#if(is.null(ret)) return(NA)
	if(is.null(ret)) ret <- NA_real_
	attr(ret, "nobs") <- length(resid(object))
	attr(ret, "df") <- length(object$coefficients) + length(object$errcomp) +
		length(object$arcoef) + 1L
	class(ret) <- "logLik"
	ret
}

`logLik.MCMCglmm` <-
function (object, ...)
	structure(NA, df = sum(object$Fixed$nfl, object$Random$nfl,
		object$Residual$nfl), nobs = object$Residual$nrl, 
		class = "logLik")


`logLik.gamm` <-
function (object, ...)
	logLik(object[[if(is.null(object$lme)) "mer" else "lme"]], ...)

	
`logLik.mark` <- 
function (object, adjust = TRUE, ...) {
	res <- -0.5 * object$results$lnl
	attr(res, "df") <- object$results[[if(!adjust && !is.null(object$results$npar.unadjusted))
		'npar.unadjusted' else 'npar']]
	attr(res, "nobs") <- object$results$n 
	class(res) <- "logLik"
	res
}

`logLik.logistf` <- 
function (object, ...) {
	res <- object$loglik[2L]
	attr(res, "nobs") <- object$n
	attr(res, "df") <- object$df + 1L
	class(res) <- "logLik"
	res
}