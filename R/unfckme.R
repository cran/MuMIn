
# vcov alternative that always returns a matrix (as it should be)
.vcov <- 
function(object, ...)
UseMethod(".vcov")

.vcov.default <- 
function(object, ...)
as.matrix(vcov(object, ...))

.vcov.glmmTMB <-
function(object, which = c("cond", "zi", "disp"), ...) {
    which <- match.arg(which)
    rval <- vcov(object, ...)[[which]]
    nm <- rownames(rval)
    nm[nm == "(Intercept)"] <- "(Int)"
    nm <- paste0(which, "(", nm, ")")
    dimnames(rval) <- list(nm, nm)
    rval
}



.numfixef <- 
function (object, ...) 
UseMethod(".numfixef")

.numfixef.default <- 
function (object, ...) fixef(object, ...)

.numfixef.cpglm <- 
function (object, ...) coef(object, ...)

.numfixef.cpglmm <- 
function (object, ...) cplm::fixef(object, ...)


.numfixef.glmmTMB <- 
function (object, ...) fixef(object, ...)$cond


# Consistent sigma (residual standard deviation)
sigma2 <-
function(object)
UseMethod("sigma2")

sigma2.default <-
function(object) {
    if(startsWith(family(object)$family, "Negative Binomial(")) {
        get(".Theta", environment(family(object)$aic))
    } else {
        sigma(object)
    }
}

sigma2.glmmPQL <-
function(object) {
    switch(family(object)$family,
        gaussian = , Gamma = object$sigma,
        object$sigma^2
    )
}
sigma2.glmmTMB <-
function(object) {
    if(family(object)$family == "nbinom1") sigma(object) + 1 else sigma(object)
}

sigma2.glmerMod <- 
function(object) {
    if(startsWith(family(object)$family, "Negative Binomial(")) {
        lme4::getME(object, "glmer.nb.theta")
    } else {
        NextMethod()
    }
}
 
    
# VarCorr wrapper returning consistent format (list of named RE variances) 
.varcorr <- 
function(object, ...) 
UseMethod(".varcorr")

.varcorr.default <- 
function(object, ...)
unclass(VarCorr(object, ...))



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


# Note: only for conditional model
.varcorr.glmmTMB <- function(object, ...) {
    unclass(VarCorr(object, ...)$cond)
}

.varcorr.glmmadmb <- function(object, ...) {
    suppressWarnings(VarCorr(object))
}
