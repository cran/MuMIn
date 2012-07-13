`tTable` <-
function (model, ...) 	{
	.Deprecated("coefTable")
	coefTable(model, ...)
}

`coefTable` <-
function (model, ...) UseMethod("coefTable")

`print.coefTable` <-
function (x, ...)
stats::printCoefmat(x, has.Pvalue = FALSE)


.makeCoefTable <- function(x, se, df = NA_real_, coefNames = names(x)) {
	if(n <- length(x)) {
		xdefined <- !is.na(x)
		ndef <- sum(xdefined)
		if(ndef < n) {
			if(length(se) == ndef) {
				y <- rep(NA_real_, n); y[xdefined] <- se; se <- y
			}
			if(length(df) == ndef) {
				y <- rep(NA_real_, n); y[xdefined] <- df; df <- y
			}
		}
	}
	if(n && n != length(se)) stop("length(x) is not equal to length(se)")
	ret <- matrix(NA_real_, ncol = 3L, nrow = length(x),
		dimnames = list(coefNames, c("Estimate", "Std. Error", "df")))
	if(n) ret[, ] <- cbind(x, se, rep(if(is.null(df)) NA_real_ else df,
		length.out = n), deparse.level = 0L)
	class(ret) <- c("coefTable", "matrix")
	ret
}

`coefTable.default` <-
function(model, ...) {
	dfs <- tryCatch(df.residual(model), error = function(e) NA_real_)
	cf <- summary(model, ...)$coefficients
	.makeCoefTable(cf[, 1L], cf[, 2L], dfs)
}

`coefTable.coxph` <-
`coefTable.survreg` <-
`coefTable.lm` <-
function(model, ...)
	.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))), model$df.residual)

`coefTable.glmmML` <- function(model, ...)
	.makeCoefTable(model$coefficients, model$coef.sd)


`coefTable.gls` <-
function (model, ...)
	.makeCoefTable(coef(model), sqrt(diag(as.matrix(model$varBeta))),
		model$dims$N - model$dims$p)


`coefTable.lme` <-
function(model, adjustSigma = TRUE, ...) {
	se <- sqrt(diag(as.matrix(model$varFix)))
	if (adjustSigma && model$method == "ML")
		se <- se * sqrt(model$dims$N / (model$dims$N - length(se)))
	.makeCoefTable(fixef(model), se, model$fixDF[["X"]])
}

`coefTable.mer` <-
function(model, ...)
	#sm <- eval(expression(summary), asNamespace("lme4"))
	.makeCoefTable(fixef(model), vcov(model, ...)@factors$correlation@sd)

`coefTable.multinom` <-
function(model, ...) {
	s <- summary(model, ...)
	cf <- s$coefficients
	se <- s$standard.errors
	.makeCoefTable(if(is.vector(cf)) cf else cf[, 1L], 
		if(is.vector(se)) se else se[, 1L])
}

`coefTable.sarlm` <-
`coefTable.spautolm` <-
function(model, ...) {
	x <- coef(model)
	.makeCoefTable(x, sqrt(diag(summary(model, ...)$resvar))[names(x)])
}

`coefTable.coxme` <-
`coefTable.lmekin` <-
function(model, ...)  {
	# code from coxme:::print.coxme
	beta <- model$coefficients # for class coxme:
	if(is.list(beta) && !is.null(beta$fixed))
		beta <- beta$fixed # for class lmekin and older coxme
	nvar <- length(beta)
	if(nvar) {
		nfrail <- nrow(model$var) - nvar
		se <- sqrt(get("diag", getNamespace("Matrix"))(model$var)[nfrail + 1L:nvar])
	} else se <- NULL
	.makeCoefTable(beta, se)
}
