## Helper functions:


# Consistent sigma
sigma2 <- function(object) UseMethod("sigma2")
sigma2.default <- function(object) sigma(object)
sigma2.glmmPQL <- function(object) {
    switch(family(object)$family,
        gaussian = , Gamma = object$sigma,
        object$sigma^2
    )
}
sigma2.glmmTMB <- function(object) {
    if(family(object)$family == "nbinom1") sigma(object) + 1 else sigma(object)
}
        
# VarCorr wrapper returning consistent format (list of named RE variances) 
.varcorr <- function(object, ...) UseMethod(".varcorr")
.varcorr.default <- function(object, ...) {
    unclass(VarCorr(object, ...))
}

# RE model matrix colnames for models with >1 random formulas are prefixed with
# the grouping factor name, e.g. :
# {~ 1 | X1, ~ 1 | X2} has model.matrix columns "X1.(Intercept)", "X2.(Intercept)"
# Need to rename VC matrix dimnames to match them.
.varcorr.lme <- function(object, ...) {
    reStruct <- object$modelStruct$reStruct
	rval <- lapply(reStruct, function(v, sig2) nlme::pdMatrix(v) * sig2, object$sigma^2)
    if ((m <- length(rval)) > 1L) {
		nm <- names(rval)
		for (i in seq.int(m)) {
			dn <- paste(nm[i], dimnames(rval[[i]])[[1L]], sep=".")
			dimnames(rval[[i]]) <- list(dn, dn)
		}
	}
    rval
}


# Note: currently ONLY FOR CONDITIONAL MODEL
.varcorr.glmmTMB <- function(object, ...) {
    unclass(VarCorr(object, ...)$cond)
}
.varcorr.glmmadmb <- function(object, ...) {
    suppressWarnings(VarCorr(object))
}

.numfixef <- function (object, ...)  UseMethod(".numfixef")
.numfixef.default <- function (object, ...) fixef(object, ...)
.numfixef.cpglm <- function (object, ...) coef(object, ...)
.numfixef.glmmTMB <- function (object, ...) fixef(object, ...)$cond

# RE model matrix
.remodmat <- function(object) UseMethod(".remodmat")

.remodmat.default <-
function(object) {
    env <- environment(formula(object))
    rval <- lapply(.findbars(formula(object)),
        function(f) model.matrix(as.formula(call("~", f[[2L]]), env = env),
                     data = model.frame(object)))
    rval <- do.call("cbind", rval)
    rval[, !duplicated(colnames(rval)), drop = FALSE]
}
    
.remodmat.merMod <-
function(object) {
    rval <- do.call("cbind", model.matrix(object, type = "randomListRaw"))
	rval[, !duplicated(colnames(rval)), drop = FALSE]
}

.remodmat.lme <-
function(object)
    model.matrix(object$modelStruct$reStruct, data = object$data[rownames(object$fitted), 
			, drop = FALSE])


.nullUpdateWarning <- 
function(message = 
"The null model is correct only if all variables used by the original model remain unchanged.",
Call = NULL) {
	if(!isTRUE(getOption("MuMIn.noUpdateWarning")))
		cry(Call, message, warn = TRUE)
}

# .nullFitRE: update `object` to intercept only model, keeping original RE terms.
# TODO: reOnlyModelCall or reOnlyFormula
.nullFitRE <- function(object, envir) UseMethod(".nullFitRE")
.nullFitRE.default <- 
function(object, envir = parent.frame()) {
    cl <- getCall(object)
	if(! "formula" %in% names(cl)) 
		stop("cannot create a null model for object without named \"formula\" argument in its call")
    cl$formula <- .nullREForm(formula(object))
	.nullUpdateWarning()
    eval(cl, envir)
}

.nullFitRE.lme <-
function(object, envir = parent.frame()) {
	cl <- getCall(object)
	cl$fixed <- update.formula(cl$fixed, . ~ 1)
    if(inherits(object, "glmmPQL")) cl$verbose <- FALSE
	.nullUpdateWarning()
	eval(cl, envir)
}

# sum up RE variance using VarCorr list
# For RE-intercept identical to a sum of diagonals of VC matrices.
# sum(sapply(lapply(vc, diag), sum))
.varRESum <- function(vc, X) {
	n <- nrow(X)
	sum(sapply(vc, function(sig) {
		mm1 <-  X[, rownames(sig), drop = FALSE]
		sum(matmultdiag(mm1 %*% sig, ty = mm1)) / n
	}))
}

# update model adding an observation level RE term
.OLREFit <- function(object) UseMethod(".OLREFit")
.OLREFit.default <- function(object) .NotYetImplemented()
.OLREFit.merMod <- function(object) {
    if (!any(sapply(object@flist, nlevels) == nobs(object))) {
		cl <- get_call(object)
        frm <- formula(object)
		nRowData <- eval(call("eval", as.expression(call("NROW", cl$formula[[2L]])),
            envir = cl$data), envir = environment(frm),
            enclos = parent.frame())
		fl <- length(frm)
		frx <- . ~ . + 1
		frx[[3L]][[3L]] <- call("(", call("|", 1, call("gl", nRowData, 1)))
		cl$formula <- update.formula(frm, frx)		
		object <- tryCatch(eval(cl, envir = environment(frm), enclos = parent.frame()),
			error = function(e) {
				cry(conditionCall(e), conditionMessage(e), warn = TRUE)
				cry(cl, "fitting model with the observation-level random effect term failed. Add the term manually")
		})
		.nullUpdateWarning("The result is correct only if all variables used by the model remain unchanged.")
    }
    object
}


# general function
r2glmm <-
function(family, vfe, vre, vol, link, pmean, lambda, omega, n) {
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
            vdelta <- (1 / lambda) + (1 / omega) # omega is theta
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
	
    #print(c(varFE = vfe, varRE = vre, varOL = vol))
    
	vtot <- sum(vfe, vre)
	matrix(c(vfe, vtot) / (vtot + rep(vol, each = 2L)),
        ncol = 2L, byrow = TRUE, dimnames = list(names(vol), c("R2m", "R2c")))
}


.binomial.sample.size <-
function(object) {
    tt <- terms(formula(object))
    y <- model.frame(object)[, rownames(attr(tt, "factors"))[attr(tt, "response")]]
    if(is.null(dim(y))) 1 else mean(rowSums(y))        
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
	
	for(i in seq.int(length(vc))) {
		a <- fixCoefNames(rownames(vc[[i]]))
		dimnames(vc[[i]]) <- list(a, a)
	}
	colnames(mmRE) <- fixCoefNames(colnames(mmRE))
    
	if(!all(unlist(sapply(vc, rownames), use.names = FALSE) %in% colnames(mmRE)))
		stop("RE term names do not match those in model matrix. \n",
			 "Have 'options(contrasts)' changed since the model was fitted?")
	
    varRE <- .varRESum(vc, mmRE) # == sum(as.numeric(VarCorr(fm)))
    
    familyName <- fam$family
    if(substr(familyName, 1L, 17L) == "Negative Binomial")
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
        lambda <- unname(exp(fixefnull + 0.5 * varRE))
        theta <- if(inherits(object, "glmerMod"))
            lme4::getME(object, "glmer.nb.theta") else
            sigma2(object)^-2 
        r2glmm(familyName, varFE, varRE, lambda = lambda, omega = theta, link = fam$link)
    }, Gamma = {
        nu <- sigma2(object)^-2
        omega <- 1
        r2glmm(fam, varFE, varRE, lambda = nu, omega = omega)
    }, quasipoisson = , nbinom1 = {
        vt <- .varRESum(.varcorr(null), mmRE)
        lambda <- unname(exp(fixefnull + 0.5 * vt))
		omega <- sigma2(object)^-2
        r2glmm(fam, varFE, varRE, lambda = lambda, omega = omega)
    }, poisson = {
        vt <- .varRESum(.varcorr(null), mmRE)
        lambda <- unname(exp(fixefnull + 0.5 * vt))
        omega <- 1
        rval <- r2glmm(fam, varFE, varRE, lambda = lambda, omega = omega)
        if(inherits(object, "merMod") &&
           familyName == "poisson" && pj2014) {
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
	fx <- fixef(object) # fixed effect estimates
	if(length(fx$zi) != 0L) # || length(fx$disp) != 0L)
		stop("r.squaredGLMM cannot (yet) handle 'glmmTMB' object with zero-inflation")
	r.squaredGLMM.merMod(object, null, envir, ...)
}

`r.squaredGLMM.glmmadmb` <-
function(object, null, envir = parent.frame(), ...) {
	if(object$zeroInflation)
		stop("r.squaredGLMM cannot (yet) handle 'glmmADMB' object with zero-inflation")
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
        r2glmm(familyName, varFE, 0, lambda = unname(exp(fixefnull)), omega = sigma2(object)^-2, link = fam$link)
    }, quasipoisson = , nbinom1 = {
        r2glmm(fam, varFE, 0, lambda = unname(exp(fixefnull)), omega = sigma2(object)^2)
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

