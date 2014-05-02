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

.makeCoefTable <- 
function(x, se, df = NA_real_, coefNames = names(x)) {
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
function (model, ...) {
	.makeCoefTable(coeffs(model), sqrt(diag(vcov(model, ...))))
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

`coefTable.rq` <- 
function(model, ...)
	.makeCoefTable(model$coefficients, rep(NA_real_, length(model$coefficients)))

	
`coefTable.zeroinfl` <-
function(model, ...)
	.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))

`coefTable.hurdle` <- 
function(model, ...) {
	cts <- summary(model)$coefficients
	ct <- do.call("rbind", unname(cts))
	cfnames <- paste(rep(names(cts), vapply(cts, nrow, 1L)), "_", rownames(ct),
		sep = "")
	.makeCoefTable(ct[, 1L], ct[, 2L], coefNames = cfnames)
	#.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))
}

`coefTable.aodql` <-
`coefTable.betareg` <- 
`coefTable.glimML` <-
`coefTable.unmarkedFit` <- 
function(model, ...)
	.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))

	
`coefTable.gee` <-
`coefTable.geeglm` <-
function(model, ..., type = c("naive", "robust")) {
	cf <- summary(model, ...)$coefficients
	type <- match.arg(type)
	j <- if(type == "naive") 2L else 4L
	.makeCoefTable(cf[, 1L], cf[, j], coefNames = rownames(cf))
}


`coefTable.geese` <-
function(model, ..., type = c("naive", "robust")) {
	cf <- summary(model, ...)$mean
	type <- match.arg(type)
	j <- if(type == "naive") 2L else 4L
	.makeCoefTable(cf[, 1L], cf[, j], coefNames = rownames(cf))
}

`coefTable.yagsResult` <-
function(model, ..., type = c("naive", "robust")) {
	type <- match.arg(type)
	vcv <- slot(model, if(type == "naive") "naive.parmvar" else "robust.parmvar")
	.makeCoefTable(model@coefficients, sqrt(diag(vcv)), coefNames = model@varnames)
}

`coefTable.splm` <- 
function (model, ...) {
	cf <- sapply(c("coefficients", "arcoef", "errcomp"), function(i)
		if(is.matrix(model[[i]])) model[[i]][, 1L] else model[[i]],
		simplify = FALSE)
	
	ncf <- sapply(cf, length)
	vcovlab <- c(coefficients = "vcov", arcoef = "vcov.arcoef", errcomp = "vcov.errcomp")
	se <- sqrt(unlist(lapply(names(vcovlab), function(i) {
		vcv2 <- diag(model[[vcovlab[i]]])
		c(vcv2, rep(NA_real_, ncf[[i]] - length(vcv2)))
	})))
	
	.makeCoefTable(unlist(cf, use.names = FALSE), se,
		coefNames = unlist(lapply(cf, names), use.names = FALSE))
}

`coefTable.MCMCglmm` <-
function (model, ...) {
	cf <- coeffs(model)
	.makeCoefTable(cf, se = rep(NA_real_, length.out = length(cf)))
}

`coefTable.gamm` <-
function (model, ...) coefTable.lm(model$gam, ...)


`coefTable.mark` <- 
function (model, orig.names = FALSE, ...) {
    dfs <- model$results[['n']] - model$results[['npar']]
    beta <- model$results[['beta']]
    .makeCoefTable(beta[, 1L], beta[, 2L], dfs,
		coefNames = if(orig.names) rownames(beta) else
			gsub("^([a-zA-Z]+):(.*)$", "\\1(\\2)", rownames(beta), perl = TRUE))
}

`coefTable.logistf` <-
function (model, ...)
.makeCoefTable(model$coefficients, sqrt(diag(model$var)))

`coefTable.aodml` <-
function (model, ...) {
	.makeCoefTable(coeffs(model), sqrt(diag(model$varparam)))
}

## XXX: fixed effects coefficients only
`coefTable.asreml` <- 
function (model, ...)  {
	.makeCoefTable(
		x = model$coefficients$fixed, 
		se = sqrt(model$vcoeff$fixed * model$sigma2) ## ?
		) 
}

`coefTable.cplm` <-
function (model, ...) 
.makeCoefTable(coef(model), sqrt(diag(vcov(model))),
			   model@df.residual)

`coefTable.cpglmm` <-
function (model, ...) 
.makeCoefTable(coeffs(model), sqrt(diag(vcov(model))))
