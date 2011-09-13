# Extends: glmmML      #########################################################

`tTable.glmmML` <- function(model, ...) {
    coef <- model$coefficients
    se <- model$coef.sd
    ret <- cbind(coef, se, coef/se, signif(1 - pchisq((coef/se)^2,
        1)))
    dimnames(ret) <- list(names(coef), c("Estimate", "Std. Error",
        "z", "Pr(>|z|)"))
  return(ret)
}

`logLik.glmmML` <- function(object, ...) {
  ll <- (-object$deviance)/2
  n <- length(object$coefficients)
  attr(ll, "df") <- n + object$cluster.null.df - object$df.residual
  attr(ll, "nobs") <- n + object$cluster.null.df
  class(ll) <- "logLik"
  ll
}

`nobs.glmmML` <- function(object, ...) length(object$coefficients) +
	object$cluster.null.df

# This is merely to get rid of the annoying behaviour in summary.glmML.
# it does not do anything except for printing the model output.
`summary.glmmML` <- function(object, ...) identity(object)

# Extends: survival    #########################################################

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

`nobs.coxph` <- function(object, ...) object$n

# Extends: lme4        #########################################################

`nobs.mer` <- function(object, nall = TRUE, ...) {
	N <- object@dims[["n"]]
	p <- object@dims[["p"]]
	if (nall) return (N)
	REML <- object@dims[['REML']]
	N - REML * p
}

# Extends: nnet/spdep        #########################################################

`nobs.sarlm` <-
`nobs.spautolm` <-
`nobs.multinom` <-
function(object, ...) NROW(fitted(object))



# No longer needed
# Extends: nlme
`nobs.gls` <- function(object, nall = TRUE, ...) {
	p <- object$dims$p
	N <- object$dims$N
	if (nall) return (N)
	REML <- object$method == "REML"
	N - REML * p
}

`nobs.lme` <- function(object, nall = TRUE, ...) {
	N <- object$dims$N
	if (nall) return (N)
	p <- object$dims$ncol[object$dims$Q + 1]
	REML <- object$method == "REML"
	N - REML * p
}

# # p - the number of coefficients in the linear model.
# #N - the number of observations in the data,
# #Q - the number of grouping levels
# #ncol - the number of columns in the model matrix for each level of grouping from innermost to outermost
# #  (last two values are equal to the number of fixed effects and one).



#limited support for unmarked

logLik.unmarkedFit <- function(object, ...) {
  ll <- -object@negLogLike
  attr(ll, "df") <- length(object@opt$par)
  attr(ll, "nobs") <- unmarked::sampleSize(object)
  class(ll) <- "logLik"
  ll
}
#setMethod("logLik", "unmarkedFit", logLik.unmarkedFit)

formula.unmarkedFit <- function (x, ...) x@formula

getAllTerms.unmarkedFit <- function (x, ...)  {
  f <- formula(x)
  t1 <- getAllTerms(f[[2]])
  t2 <- getAllTerms(f[-2])
  structure(c(sprintf("psi(%s)",t1), sprintf("p(%s)",t2)), 
            intercept=c(attr(t1, "intercept"), attr(t2, "intercept")))
}

tTable.unmarkedFit <- function (model, ...) {
  do.call("rbind", lapply(model@estimates@estimates, function(y) {
    ret <- cbind(Estimate=y@estimates, SE=sqrt(diag(y@covMat)))
    rn <- rownames(ret)
    rn[rn=="(Intercept)"] <- "Int"
    rownames(ret) <- paste(y@short.name, "(", rn, ")", sep="")
    ret 
  }))
}

coefDf.unmarkedFit <- function(x) rep(NA, length(coef(x)))
nobs.unmarkedFit <- function(object, ...) unmarked::sampleSize(object)
