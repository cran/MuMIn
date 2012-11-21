# Code based on 'compar.gee' from package 'ape'
## Comparative Analysis with GEEs
## compar.gee.R (2011-06-14)
## Copyright 2002-2010 Emmanuel Paradis
## https://svn.mpl.ird.fr/ape/dev/ape/R/compar.gee.R

##=============================================================================
## quasiLik 
##=============================================================================


`quasiLik` <- function (object, ...) UseMethod("quasiLik")

.qlik <- function(y, mu, fam) {
	ret <- switch(fam,
		   gaussian = -sum((y - mu)^2)/2,
		   binomial = sum(y * log(mu/(1 - mu)) + log(1 - mu)),
		   #binomial.sqvar = sum(((2 * y - 1) * log(mu /(1 - mu))) - (y / mu) - ((1 - y)/(1 - mu))),
		   poisson = sum(y * log(mu) - mu),
		   Gamma = -sum(y/mu + log(mu)),
		   inverse.gaussian = sum(-y/(2 * mu^2) + 1/mu),
		   stop("do not know how to calculate quasi-likelihood for family ",
				dQuote(fam))
		   )
	ret
}

`print.quasiLik` <- function (x, digits = getOption("digits"), ...) {
    cat("'quasi Lik.' ", paste(format(c(x), digits = digits), collapse = ", "), 
        "\n", sep = "")
    invisible(x)
}

`quasiLik.geeglm` <-
`quasiLik.gee` <-
function(object, ...) {
	ret <- .qlik(object$y, object$fitted.values, family(object)$family)
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(object$y)
	class(ret) <- "quasiLik"
	ret
}

`quasiLik.yagsResult` <- function(object, ...) {
	mu <- object@fitted.values
	ret <- .qlik(mu + object@residuals, mu, family(object)$family)
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(mu)
	class(ret) <- "quasiLik"
	ret
}


##=============================================================================
## QIC 
##=============================================================================
.qic <- function(mu, vbeta, i.vbeta.naiv, qlik) {
	AIinv <- solve(i.vbeta.naiv) # solve via indenity
	tr <- sum(diag(AIinv %*% vbeta))
	px <- length(mu) # number non-redunant columns in design matrix
	# QIC
	# ret <- -2 * qlik + 2 * tr
	ret <- 2 * (tr - qlik)
	QICu <- 2 * (px - qlik)    # Approximation assuming model structured correctly
	#attr(ret, "QICu") <- QICu
	c(ret, QICu)
}

`getQIC` <- function(x, typeR = FALSE) UseMethod("getQIC")
	
`getQIC.default` <-
function(x, typeR = FALSE) .NotYetImplemented()

	
`getQIC.gee` <- function(x, typeR = FALSE) {
	capture.output(suppressMessages(xi <- update(x, corstr = "independence",
		silent = TRUE)))
	qx <- if(typeR) x else xi
	ql <- .qlik(qx$y, qx$fitted.values, family(qx)$family)
	n <- length(x$y)
	# yags/yags.cc: p140 of Hardin and Hilbe   
	ql <- if(family(x)$family == "gaussian")
			(n * log(-2 * ql / n)) / -2	 else ql
	c(.qic(x$fitted.values, x$robust.variance, xi$naive.variance, ql), n)
}

`getQIC.geeglm` <- function(x, typeR = FALSE) {
	xi <- update(x, corstr = "independence")
	qx <- if(typeR) x else xi
	ql <- .qlik(qx$y, qx$fitted.values, family(qx)$family)
	n <- length(x$y)
	# yags/yags.cc: p140 of Hardin and Hilbe
	ql <- if(family(x)$family == "gaussian")
			(n * log(-2 * ql / n)) / -2	 else ql
	c(.qic(x$fitted.values, x$geese$vbeta, xi$geese$vbeta.naiv, ql), n)
}

`getQIC.yagsResult` <- function(x, typeR = FALSE) {
	xi <- update(x, corstruct = "independence")
	##
	#cl <- match.call(call = getCall(yags1), yags::yags)
	#cl[[1L]] <- as.name("model.frame.default")
	#cl$formula[[3L]] <- 1L
	#cl <- cl[c(TRUE, (names(cl)[-1L] %in% c("formula", "data", "subset")))]
	#y <- eval(cl, parent.frame())[, 1L]
	qx <- if(typeR) x else xi
	y <- qx@fitted.values + qx@residuals
	n <- length(y)
	ql <- .qlik(y, qx@fitted.values, family(qx)$family)
		# yags/yags.cc: p140 of Hardin and Hilbe
	ql <- if(family(x)$family == "gaussian")
			(n * log(-2 * ql / n)) / -2	 else ql
	c(.qic(x@fitted.values, x@robust.parmvar, xi@naive.parmvar, ql), n)
}

`QIC` <- function (object, ..., typeR = FALSE) {
	if (length(list(...))) {
		res <- sapply(list(object, ...), getQIC, typeR = typeR)
		#val <- data.frame(QIC = res[1L, ], QICu = res[2L, ])
		val <- as.data.frame(t(res[1L:2L, ]))
		colnames(val) <- c("QIC", "QICu")
		Call <- match.call()
		Call$typeR <- NULL
		row.names(val) <- as.character(Call[-1L])
		val
	} else getQIC(object, typeR = typeR)[1L]
}
