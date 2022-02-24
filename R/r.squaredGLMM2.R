
# general function
r2glmm <-
function(family, vfe, vre, vol, link, pmean, lambda, omega, n) {
    
    #message("lambda=", lambda)
    #message("omega=", omega)
    
    if(inherits(family, "family")) {
        link <- family$link
        family <- family$family
    }
    
    if(missing(vol) || !is.numeric(vol) || is.na(vol))
    vol <- switch(paste(family, link, sep = "."),
        gaussian.identity = vol,
        quasibinomial.logit =,
        binomial.logit = c(
            theoretical = 3.28986813369645 / n,
            delta = 1 / (n * pmean * (1 - pmean))
            ),
        quasibinomial.probit =,
        binomial.probit = c(
            theoretical = 1 / n,
            delta =
                6.2831853071795862 / n * pmean * (1 - pmean) *
                    exp((qnorm(pmean) / 1.4142135623730951)^2)^2
            ),
        quasibinomial.cloglog =,
		binomial.cloglog = c(
			theoretical = 1.6449340668482264 / n,  #  pi^2 / 6
			delta = pmean / n / log(1 - pmean)^2 / (1 - pmean)
			),
		Gamma.log =, poisson.log =, quasipoisson.log =, nbinom1.log = c(
			delta = omega / lambda,
			lognormal = log1p(omega / lambda),
			trigamma = trigamma(lambda / omega)
			),
        quasipoisson.sqrt =, nbinom1.sqrt =,
        "poisson.mu^0.5" =, poisson.sqrt = c(
            delta = 0.25 * omega
            ),
        nbinom2.log = {
            vdelta <- (1 / lambda) + (1 / omega)
            c(
            delta = vdelta,
            lognormal =  log1p(vdelta),
            trigamma = trigamma(1 / vdelta)
            )},
		Gamma.inverse =, # c( delta = 1 / nu / lambda^2 ),
        NotImplementedFamily =
            stop("not implemented yet for ", family, " and ", link),
        cry(sys.call(-1L), "do not know how to calculate variance for %s(%s)", family, dQuote(link))
    )
	
	vtot <- sum(vfe, vre)
	matrix(c(vfe, vtot) / (vtot + rep(vol, each = 2L)),
        ncol = 2L, byrow = TRUE, dimnames = list(names(vol), c("R2m", "R2c")))
}



`r.squaredGLMM` <-
function(object, null, ...) {
    warnonce("rsquaredGLMM",
	simpleWarning(paste0("'r.squaredGLMM' now calculates a revised statistic. See the help page.")))
    UseMethod("r.squaredGLMM")
}

`r.squaredGLMM.merMod` <-
function(object, null, envir = parent.frame(), pj2014 = FALSE, ...) {
    
    if(is.logical(envir)) { # backwards compatibility
        tmp <- envir
        if(!missing(pj2014)) envir <- pj2014
        pj2014 <- tmp
    }
    
	fam <- family(object)
    #varOL <- lambda <- omega <- NA
    fe <- .numfixef(object)
    ok <- !is.na(fe)
    fitted <- (model.matrix(object)[, ok, drop = FALSE] %*% fe[ok])[, 1L]
    varFE <- var(fitted)
	
   
	mmRE <- .remodmat(object)
    ##Note: Argument 'contrasts' can only be specified for fixed effects
	##contrasts.arg = eval(cl$contrasts, envir = environment(formula(object))))	
	vc <- .varcorr(object)
    if(!is.null(vc)) {
        for(i in seq.int(length(vc))) {
            a <- fixCoefNames(rownames(vc[[i]]))
            dimnames(vc[[i]]) <- list(a, a)
        }
        colnames(mmRE) <- fixCoefNames(colnames(mmRE))
        if(!all(unlist(sapply(vc, rownames), use.names = FALSE) %in% colnames(mmRE)))
            stop("RE term names do not match those in model matrix. \n",
                 "Have 'options(contrasts)' changed since the model was fitted?")
        varRE <- .varRESum(vc, mmRE) # == sum(as.numeric(VarCorr(fm)))
    } else {
        varRE <- 0
    }
    
    familyName <- fam$family
    if(startsWith(familyName, "Negative Binomial("))
        familyName <- "nbinom2"
    
    if(familyName %in% c("quasipoisson", "poisson", "nbinom1", "nbinom2",
        "binomial", "quasibinomial")) {
		if(missing(null) || !is.object(null)) null <- .nullFitRE(object, envir)
        fixefnull <- unname(.numfixef(null))
    }
    
    switch(familyName,
    gaussian =
        r2glmm(fam, varFE, varRE, vol = sigma2(object)^2),
    binomial =, quasibinomial = {
        vt <- .varRESum(.varcorr(null), mmRE)
		# XXX: inverse-link seems to give more reasonable value for non-logit
		# links, should inv-logit (plogis) be used here always?
        pmean <- fam$linkinv(fixefnull - 0.5 * vt * tanh(fixefnull * (1 + 2 * exp(-0.5 * vt)) / 6))
        r2glmm(fam, varFE, varRE, pmean = pmean, n = .binomial.sample.size(object))
    }, nbinom2 = {
        vt <- .varRESum(.varcorr(null), mmRE)
        lambda <- unname(exp(fixefnull + 0.5 * vt))
        theta <- sigma2(object)
        
        #koBrowseHere()

        r2glmm(familyName, varFE, varRE, lambda = lambda, omega = theta, link = fam$link)

    }, Gamma = {
        nu <- sigma2(object)^-2
        omega <- 1
        r2glmm(fam, varFE, varRE, lambda = nu, omega = omega)
    }, quasipoisson = , nbinom1 = {
        vt <- .varRESum(.varcorr(null), mmRE)
        lambda <- unname(exp(fixefnull + 0.5 * vt))
        omega <- sigma2(object)
        r2glmm(fam, varFE, varRE, lambda = lambda, omega = omega)
    }, poisson = {
        vt <- .varRESum(.varcorr(null), mmRE)
        lambda <- unname(exp(fixefnull + 0.5 * vt))
        omega <- 1
        rval <- r2glmm(fam, varFE, varRE, lambda = lambda, omega = omega)
        if(inherits(object, "merMod") && familyName == "poisson" && pj2014) {
            xo <- .OLREFit(object)
            vc <- .varcorr(xo)
            fe <- .numfixef(xo)
            ok <- !is.na(fe)
            fitted <- (model.matrix(xo)[, ok, drop = FALSE] %*% fe[ok])[, 1L]
            
            n <- nrow(mmRE)
            vname <- names(xo@flist)[sapply(xo@flist, nlevels) == n][1L]
            if(! vname %in% names(vc)) vname <- make.names(vname)
            stopifnot(vname %in% names(vc)) ### !!!
            varresid <- vc[[vname]][1L]
            rval <- rbind(pj2014 = r2glmm(fam, var(fitted),
                .varRESum(vc, mmRE) - varresid,
                vol = log1p(1 / exp(mean(fitted))) + varresid)[1L, ], rval)
        }
        rval
    }, r2glmm(fam, varFE, varRE))
}


`r.squaredGLMM.lme` <-
function(object, null, ...) r.squaredGLMM.merMod(object, null, ...)

`r.squaredGLMM.glmmTMB` <-
function(object, null, envir = parent.frame(), ...) {
    w <- c(fixef(object)$zi != 0L,
        !identical(object$modelInfo$allForm$dispformula, ~0, ignore.environment = TRUE))
    if(any(w)) warning("the effects of ",
            prettyEnumStr(c("zero-inflation", "dispersion model"), quote = FALSE), " are ignored")
    r.squaredGLMM.merMod(object, null, envir, ...)
}

`r.squaredGLMM.glmmadmb` <-
function(object, null, envir = parent.frame(), ...) {
	if(object$zeroInflation)
        warning("effects of zero-inflation are ignored")
		#stop("r.squaredGLMM cannot (yet) handle 'glmmADMB' object with zero-inflation")
	r.squaredGLMM.merMod(object, null, envir, ...)
}

`r.squaredGLMM.lm` <-
function(object, null, envir = parent.frame(), ...) {
	fam <- family(object)
    ok <- !is.na(coef(object))
    fitted <- (model.matrix(object)[, ok, drop = FALSE] %*% coef(object)[ok])[, 1L]
	delayedAssign("fixefnull",
		coef(if(missing(null) || !is.object(null))
			 .nullFitRE(object, envir) else null))
    varFE <- var(fitted)
    familyName <- fam$family
    if(substr(familyName, 1L, 17L) == "Negative Binomial")
        familyName <- "nbinom2"
    
    switch(familyName,
    gaussian =
        r2glmm(fam, varFE, 0, vol = sigma2(object)^2),
    binomial =, quasibinomial = {
        r2glmm(fam, varFE, 0, pmean = fam$linkinv(unname(fixefnull)),
            n = .binomial.sample.size(object))
    }, Gamma = {
        nu <- sigma2(object)^-2
        omega <- 1
        r2glmm(fam, varFE, 0, lambda = nu, omega = omega)
    }, nbinom2 = {
        r2glmm(familyName, varFE, 0, lambda = unname(exp(fixefnull)),
            omega = sigma2(object), link = fam$link)
    }, quasipoisson = , nbinom1 = {
        r2glmm(fam, varFE, 0, lambda = unname(exp(fixefnull)),
            omega = sigma2(object))
    }, poisson = {
        r2glmm(fam, varFE, 0, lambda = unname(exp(fixefnull)), omega = 1)
    }, r2glmm(fam, varFE, 0))
}


# TODO 
`r.squaredGLMM.glmmML` <-
function(object, null, ...) {
    .NotYetImplemented()
}

`r.squaredGLMM.cplm` <-
function(object, null, envir = parent.frame(), ...) {
	fam <- family(object)
    if(!fam$link %in% c("mu^0", "log"))
         stop("not implemented yet for ", fam$family, " and ", fam$link)
         
    fe <- .numfixef(object)
    ok <- !is.na(fe)
    fitted <- (model.matrix(object)[, ok, drop = FALSE] %*% fe[ok])[, 1L]
    varFE <- var(fitted)
	if(missing(null) || !is.object(null)) null <- .nullFitRE(object, envir)

	if(inherits(object, "cpglm")) {
		varRE <- vt <- 0
	} else {
		mmRE <- .remodmat(object)
		varRE <- .varRESum(.varcorr(object), mmRE) # == sum(as.numeric(VarCorr(fm)))
		vt <- .varRESum(.varcorr(null), mmRE)
	}
	mu <- unname(exp(.numfixef(null) + 0.5 * vt)) # the same as getting lambda
    phi <- object@phi # the dispersion parameter
    p <- object@p # the index parameter
    varO <- c(delta = phi * mu^(p - 2), lognormal = log1p(phi * mu^(p - 2)))
	r2glmm(NA, varFE, varRE, varO)
}

