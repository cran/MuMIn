# Code based on 'compar.gee' from package 'ape'
## Comparative Analysis with GEEs
## compar.gee.R (2011-06-14)
## Copyright 2002-2010 Emmanuel Paradis
## https://svn.mpl.ird.fr/ape/dev/ape/R/compar.gee.R

##=============================================================================
## quasiLik 
##=============================================================================


`quasiLik` <- function (object, ...) UseMethod("quasiLik")

# TODO: add weights to families
# TODO: update calls to .qlik, add weights
.qlik <- function(y, mu, fam, w, scale = 1) {
    w <- w / (sum(w) / length(w))
	switch(fam$family,
		gaussian = -0.5 * sum(w * (y - mu)^2) / scale,
		binomial = sum(w * (y * log(mu / (1 - mu)) + log(1 - mu))),
		#binomial.sqvar = sum(((2 * y - 1) * log(mu /(1 - mu))) - (y / mu) - ((1 - y)/(1 - mu))),
		poisson = sum(w * ((y * log(mu)) - mu)),
		Gamma = -sum(w * (y / mu + log(mu))) / scale,
		inverse.gaussian =  sum(w * (mu - 0.5 * y) / mu^2) / scale,
        {
            # negative.binomial = sum((y * log(mu)) - (2 * log(mu + 1))) / scale,
            if(startsWith(tolower(fam$family), "negative binomial(")) {
                gt <- fam$getTheta
                th <- if(is.function(gt)) {
                    if(is.null(formals(gt)$trans)) gt() else gt(TRUE)
                } else environment(fam$aic)$.Theta
                a <- (1 / th) * mu
                ap1 <- a + 1
                sum(w * (lgamma(y + th) - lgamma(th) + y * log(a / ap1) +
                    th * log(1 / ap1)))
            } else
                cry(, "do not know how to calculate quasi-likelihood for family '%s'",
                    fam$family)
		})
}

`print.quasiLik` <-
function (x, digits = getOption("digits"), ...) {
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
	ret <- .qlik(object$y, object$fitted.values, family(object),
        1, # XXX should use weights(object) for 'geeglm', 'gee' gives no weights,
        scale)
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(object$y)
	class(ret) <- "quasiLik"
	ret
}

# XXX: check weights
`quasiLik.yagsResult` <- function(object, ...) {
	mu <- object@fitted.values
	ret <- .qlik(mu + object@residuals, mu, family(object),
         1, # 'yags' object gives no weights,
         object@phi)
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(mu)
	class(ret) <- "quasiLik"
	ret
}

`quasiLik.geem` <-
function(object, ...) {
	fam <- family(object)
	scale <- object$phi
	ret <- .qlik(object$y, fitted(object),
        if(inherits(fam, "family")) fam else list(family = "custom"),
        1, # object$weights,
        scale)
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
		family(object),
        1, # object$weight XXX" no -s,
        object$scale[[1L]])
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(object$mu_fit)
	class(ret) <- "quasiLik"
	ret
}

##=============================================================================
## QIC 
##=============================================================================

.qic2 <- function(y, mu, vbeta, mui, vbeta.naiv.i, fam, wts, scale, typeR = FALSE) {
    ql <- .qlik(y, if(typeR) mu else mui, fam, wts, scale)
	# XXX: should be typeR = TRUE for QICu???
	n <- length(y)
    
    invert <- if (is.matrix(vbeta.naiv.i) && 
        "MASS" %in% loadedNamespaces()) MASS::ginv else solve
    
	AIinv <- invert(vbeta.naiv.i)
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
	#warnonce("getQIC.coxph", "QIC for 'coxph' is experimental")
	warnonce("QIC for 'coxph' is experimental")
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
        
    y <- if(x$family$family == "binomial" && any(x$y > 1)) {
        cl <- getCall(x)
        cl <- cl[names(cl) %in% c("", "formula", "data", "subset", "na.action")]
        cl[[1L]] <- as.name("model.frame")
        mf <- eval.parent(cl)
        warning(gettextf("using response \"%s\" from the current %s.",
            names(mf)[1L], if(is.null(cl$data)) "environment" else sQuote(asChar(cl$data))))
        y <- mf[, 1L]
        y <- y[, 1L] / rowSums(y)
    } else x$y
    
    .qic2(y, x$fitted.values, x$robust.variance, 
		  xi$fitted.values, xi$naive.variance, family(x),
          1, # no weights
		  scale = x$scale,
		  typeR = typeR)
}

`getQIC.geeglm` <- 
function(x, typeR = FALSE) {
	xi <- if(x$corstr != "independence") {
        cl <- getCall(x)
        cl$corstr <- "independence"
        cl$zcor <- NULL
        eval.parent(cl)        
        } else x
	.qic2(x$y, x$fitted.values, x$geese$vbeta, 
		  xi$fitted.values, xi$geese$vbeta.naiv, family(x),
          1, # weights(x)
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
		  xi@fitted.values, xi@naive.parmvar,
          family(x),
          1, # no weights
		  scale = x@phi,
		  typeR = typeR)
}

`getQIC.geem` <-
function(x, typeR = FALSE) {
	fam <- family(x)
	xi <- if(x$corr != "independence")
		update(x, corstr = "independence") else x
    .qic2(x$y, fitted(x), x$var, fitted(xi), xi$naiv.var,
        if (inherits(fam, "family")) fam else list(family = "custom"),
        1, # x$weights
        scale = x$phi,
        typeR = typeR
    )
}

`QIC` <- function(object, ..., typeR = FALSE) {
    if (!missing(...)) {
        res <- sapply(list(object, ...), getQIC, typeR = typeR)
        val <- as.data.frame(t(res[1L, , drop = FALSE]))
        colnames(val) <- c("QIC")
        Call <- match.call()
        Call$typeR <- NULL
        row.names(val) <- as.character(Call[-1L])
        val
    } else
        getQIC(object, typeR = typeR)[1L]
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
