`DIC` <-
function (object, ...) {
	if (length(list(...))) {
		lls <- sapply(list(object, ...), function(x) {
			c(extractDIC(x), attr(logLik(x), "df"))
		})
		val <- data.frame(df = lls[2L, ], DIC = lls[1L, ])
		Call <- match.call()
		row.names(val) <- make.unique(as.character(Call[-1L]))
		val
	} else extractDIC(object)
}

if(!exists("extractDIC", mode = "function")) {
	extractDIC <- function (fit, ...) UseMethod("extractDIC")
}

# from package 'arm'
`extractDIC.mer` <- function (fit, ...) {
	llik <- logLik(fit, REML = fit@dims["REML"])
    #dev <- fit@deviance["ML"]
    dev <- deviance(fit, REML = FALSE)
    c(2 * (dev + llik))[[1L]]
    #llik <- logLik(fit, fit@dims["REML"])
    #dev <- fit@deviance["ML"]
    ##n <- fit@dims["n"]
    #Dhat <- -2 * c(llik)
    #pD <- dev - Dhat
    #DIC <- dev + pD[[1]]
    #return(DIC)
}

`extractDIC.MCMCglmm` <- function (fit, ...) fit$DIC

`extractDIC.lme` <- function (fit, ...) .NotYetImplemented()


`formula.MCMCglmm` <-
function (x, ...) x$Fixed$formula

`nobs.MCMCglmm` <-
function (object, ...) object$Residual$nrl

`family.MCMCglmm` <-
function (object, ...) object$family

`logLik.MCMCglmm` <-
function (object, ...)
	structure(NA, df = object$Fixed$nfl + object$Random$nfl,
			  nobs = object$Residual$nrl, class = "logLik")

`coeffs.MCMCglmm` <-
function (model) summary(model)$solutions[, 1L]

`coefTable.MCMCglmm` <-
function (model, ...) {
	cf <- coeffs(model)
	.makeCoefTable(cf, se = rep(NA_real_, length.out = length(cf)))
}
