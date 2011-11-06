# Extends: survival    #########################################################


`logLik.glmmML` <- function(object, ...) {
  ll <- (-object$deviance)/2
  n <- length(object$coefficients)
  attr(ll, "df") <- n + object$cluster.null.df - object$df.residual
  attr(ll, "nobs") <- n + object$cluster.null.df
  class(ll) <- "logLik"
  ll
}


# logLik for survival::coxph model
# https://stat.ethz.ch/pipermail/r-help/2006-December/122118.html
# originally by Charles C. Berry, mod. by KB: correction for the null model
`logLik.coxph` <- function(object,...) {
# Thx to Mathieu Basille:
    y <-  object$loglik[length(object$loglik)]
	#y <-  if(length(object$loglik) > 1)
	#	-1 * (object$loglik[1] - object$loglik[2]) else object$loglik
    class(y) <- "logLik"
	#attr(y,"nall") <-
	attr(y, "nobs") <- object$n
    attr(y,'df') <- if(is.null(object$coef)) 0 else sum(!is.na(object$coef))
    return(y)
}

`logLik.survreg` <- function (object, ...)  {
	ret <- object$loglik[2L]
	attr(ret, "nobs") <- length(object$linear)
	attr(ret, "df") <- sum(object$df)
	ret
}

