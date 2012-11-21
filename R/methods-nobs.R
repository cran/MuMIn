`nobs.glmmML` <- function(object, ...) length(object$coefficients) +
	object$cluster.null.df

# Extends: nlme --- No longer needed
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

# Extends: survival
`nobs.survreg` <- function (object, ...) {
	length(object$linear)
}

`nobs.coxph` <- function(object, ...) object$n

# Extends: lme4
`nobs.mer` <- function(object, nall = TRUE, ...) {
	N <- object@dims[["n"]]
	p <- object@dims[["p"]]
	if (nall) return (N)
	REML <- object@dims[['REML']]
	N - REML * p
}

# Extends: nnet/spdep
`nobs.sarlm` <-
`nobs.spautolm` <-
`nobs.multinom` <-
function(object, ...) NROW(fitted(object))
