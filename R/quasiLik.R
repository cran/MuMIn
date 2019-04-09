# Code based on 'compar.gee' from package 'ape'
## Comparative Analysis with GEEs
## compar.gee.R (2011-06-14)
## Copyright 2002-2010 Emmanuel Paradis
## https://svn.mpl.ird.fr/ape/dev/ape/R/compar.gee.R

##=============================================================================
## quasiLik 
##=============================================================================

`quasiLik` <- function (object, ...) UseMethod("quasiLik")

.qlik <- function(y, mu, fam, scale = 1) {
	switch(fam,
		gaussian = -0.5 * sum((y - mu)^2) / scale,
		binomial = sum(y * log(mu / (1 - mu)) + log(1 - mu)),
		#binomial.sqvar = sum(((2 * y - 1) * log(mu /(1 - mu))) - (y / mu) - ((1 - y)/(1 - mu))),
		poisson = sum((y * log(mu)) - mu),
		Gamma = -sum(y / mu + log(mu)) / scale,
		inverse.gaussian = sum(-(y / (2 * mu^2)) + (1 / mu)) / scale,
        negbin = , 
        negative.binomial = sum((y * log(mu)) - (2 * log(mu + 1))) / scale,
		cry(, "do not know how to calculate quasi-likelihood for family '%s'",
			 fam))
}

`print.quasiLik` <- function (x, digits = getOption("digits"), ...) {
    cat("'quasi Lik.' ", paste(format(c(x), digits = digits), collapse = ", "), 
        "\n", sep = "")
    invisible(x)
}



`quasiLik.geeglm` <-
`quasiLik.gee` <-
function(object, ...) {
	scale <- if(is.null(object$geese))
		object$scale else
			object$geese$gamma[[1L]]
	ret <- .qlik(object$y, object$fitted.values, family(object)$family, scale)
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(object$y)
	class(ret) <- "quasiLik"
	ret
}

`quasiLik.yagsResult` <- function(object, ...) {
	mu <- object@fitted.values
	ret <- .qlik(mu + object@residuals, mu, family(object)$family, object@phi)
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(mu)
	class(ret) <- "quasiLik"
	ret
}

`quasiLik.geem` <-
function(object, ...) {
	fam <- family(object)
	scale <- object$phi
	ret <- .qlik(object$y, fitted(object), if(inherits(fam, "family"))
				 fam$family else "custom", scale)
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(object$y)
	class(ret) <- "quasiLik"
	ret
}


quasiLik.wgee <-
function(object, ...) {
	dat <- match.fun(getOption('na.action'))(model.frame(object$model, data = object$data))
	bad <- attr(dat, "na.action")
	ret <- .qlik(if(!is.null(bad)) object$y[-bad] else object$y, object$mu_fit,
		family(object)$family, object$scale[[1L]])
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(object$mu_fit)
	class(ret) <- "quasiLik"
	ret
}

##=============================================================================
## QIC 
##=============================================================================

.qic2 <- function(y, mu, vbeta, mui, vbeta.naiv.i, fam, scale, typeR = FALSE) {
	#ql <- if(typeR) .qlik(y, mu, fam, scale) else .qlik(y, mui, fam, scale)
    ql <- .qlik(y, if(typeR) mu else mui, fam, scale)
	# XXX: should be typeR = TRUE for QICu???
	n <- length(y)
	AIinv <- solve(vbeta.naiv.i)
	tr <- sum(matmult(AIinv, vbeta, diag.only = TRUE)) 
	## tr <- sum(diag(AIinv %*% vbeta))
	#px <- length(mu)
	px <- dim(vbeta)[1L]
	## When all modelling specifications in GEE are correct tr = px.
	c(2 * (c(QIC = tr, QICu = px) - ql), n = n)
}

`getQIC` <- 
function(x, typeR = FALSE) UseMethod("getQIC")
	
`getQIC.default` <-
function(x, typeR = FALSE) .NotYetImplemented()

`getQIC.coxph` <- function(x, ...) {
	warnonce("getQIC.coxph", "QIC for 'coxph' is experimental")
	naive.var <- x[[ if (is.null(x$naive.var)) "var" else "naive.var" ]]
	# tr <- sum(diag(solve(naive.var) %*% x$var))
	tr <- sum(matmultdiag(solve(naive.var), x$var))
	ll <- x$loglik[2L]
	px <- dim(x$var)[1L]
	c(2 * (c(QIC = tr, QICu = px) - ll), n = length(x$y))
}

`getQIC.gee` <- 
function(x, typeR = FALSE) {
	if(x$model$corstr != "Independent")
		utils::capture.output(suppressMessages(xi <- update(x, corstr = "independence",
		silent = TRUE))) else
		xi <- x
	.qic2(x$y, x$fitted.values, x$robust.variance, 
		  xi$fitted.values, xi$naive.variance, family(x)$family,
		  scale = x$scale,
		  typeR = typeR)
}

`getQIC.geeglm` <- 
function(x, typeR = FALSE) {
	xi <- if(x$corstr != "independence")
		update(x, corstr = "independence") else x
	.qic2(x$y, x$fitted.values, x$geese$vbeta, 
		  xi$fitted.values, xi$geese$vbeta.naiv, family(x)$family,
		  scale = x$geese$gamma[[1L]],
		  typeR = typeR)
}

`getQIC.wgee` <- 
function(x, typeR = FALSE) {
	if(isTRUE(typeR)) warning("argument 'typeR' ignored.")
	qic <- getFrom("wgeesel", "QIC.gee")(x)
	c(QIC = qic[[1L]], QICu = qic[[2L]], n = length(x$mu_fit))
}

`getQIC.yagsResult` <- 
function(x, typeR = FALSE) {
	xi <- if(x@corstruct.tag != "independence")
		update(x, corstruct = "independence") else x
	.qic2(x@fitted.values + x@residuals, x@fitted.values, x@robust.parmvar, 
		  xi@fitted.values, xi@naive.parmvar, family(x)$family,
		  scale = x@phi,
		  typeR = typeR)
}

`getQIC.geem` <-
function(x, typeR = FALSE) {
	fam <- family(x)
	xi <- if(x$corr != "independence")
		update(x, corstr = "independence") else x
    .qic2(x$y, fitted(x), x$var, fitted(xi), xi$naiv.var,
		if(inherits(fam, "family")) fam$family else "custom",
		scale = x$phi,
        typeR = typeR)
}

`QIC` <- function (object, ..., typeR = FALSE) {
	if (!missing(...)) {
		res <- sapply(list(object, ...), getQIC, typeR = typeR)
		val <- as.data.frame(t(res[1L,, drop = FALSE]))
		colnames(val) <- c("QIC")
		Call <- match.call()
		Call$typeR <- NULL
		row.names(val) <- as.character(Call[-1L])
		val
	} else getQIC(object, typeR = typeR)[1L]
}

`QICu` <- function (object, ..., typeR = FALSE) {
	if (!missing(...)) {
		res <- sapply(list(object, ...), getQIC, typeR = typeR)
		val <- as.data.frame(t(res[2L,, drop = FALSE]))
		colnames(val) <- "QICu"
		Call <- match.call()
		Call$typeR <- NULL
		row.names(val) <- as.character(Call[-1L])
		val
	} else getQIC(object, typeR = typeR)[2L]
}
