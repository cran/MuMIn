`DIC` <-
function (object, ...) {
	if (length(list(...))) {
		lls <- sapply(list(object, ...), function(x)
			c(x$DIC, x$Fixed$nfl + x$Random$nrl))
		val <- data.frame(df = lls[2L, ], DIC = lls[1L, ])
		Call <- match.call()
		Call$k <- Call$REML <- NULL
		row.names(val) <- as.character(Call[-1L])
		val
	} else object$DIC
}

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
