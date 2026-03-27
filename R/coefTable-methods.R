
`coefTable.default` <-
function(model, ...) {
	dfs <- tryCatch(df.residual(model), error = function(e) NA_real_)
	cf <- summary(model, ...)$coefficients
	.makeCoefTable(cf[, 1L], cf[, 2L], dfs, coefNames = rownames(cf))
}

`coefTable.lm` <-
`coefTable.betareg` <- 
function(model, ...)
	.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))), df.residual(model))

`coefTable.survreg` <- 
function(model, ...) {
.makeCoefTable(
	coeffs(model), 
	sqrt(diag(vcov(model, ...))), 
	NA,
	)
}

`coefTable.coxph` <-
function(model, ...) {
	.makeCoefTable(coef(model), if(all(is.na(model$var))) 
		rep(NA_real_, length(coef(model))) else sqrt(diag(model$var)), 
		model$df.residual)
}
	
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
	.makeCoefTable(nlme::fixef(model), se, model$fixDF[["X"]])
}

`coefTable.multinom` <- 
function (model, ...) {
	.makeCoefTable(coeffs(model), sqrt(diag(vcov(model, ...))))
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
	cfnames <- paste0(rep(names(cts), vapply(cts, nrow, 1L)), "_", rownames(ct))
	.makeCoefTable(ct[, 1L], ct[, 2L], coefNames = cfnames)
	#.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))
}

`coefTable.aodql` <-
`coefTable.glimML` <-
function(model, ...)
	.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))


coefTable.unmarkedFit <-
function (model, ...) {
    vcv <- umf_vcov(model, ...)
    cf <- coef(model)
    se <- if(is.null(vcv)) rep(NA_real_,length(cf)) else sqrt(diag(vcv))
    .makeCoefTable(cf, se)
}


`coefTable.gee` <-
`coefTable.geeglm` <-
function(model, ..., type = c("naive", "robust")) {
	cf <- summary(model, ...)$coefficients
	j <- if(match.arg(type) == "naive") 2L else 4L
	.makeCoefTable(cf[, 1L], cf[, j], coefNames = rownames(cf))
}

`coefTable.geem` <-
function(model, ..., type = c("naive", "robust")) {
	smr <- summary(model)
	.makeCoefTable(smr$beta, smr[[if(match.arg(type) == "naive")
								  "se.model" else "se.robust"]],
	               coefNames = smr$coefnames)
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
	.makeCoefTable(coeffs(model), sqrt(diag(vcov(model))))
	#.makeCoefTable(coeffs(model), sqrt(diag(model$varparam)))
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


`coefTable.maxlikeFit` <-
function (model, ...)
.makeCoefTable(model$Est[, 1L], model$Est[, 2L])


coefTable.bic.glm <-
function (model, ...) {
	.makeCoefTable(model$condpostmean, model$condpostsd, NA_integer_,
		dimnames(model$mle)[[2L]])
}
