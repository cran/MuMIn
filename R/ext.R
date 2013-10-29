# This is merely to get rid of the annoying behaviour in summary.glmML.
# it does not do anything except for printing the model output.
`summary.glmmML` <- function(object, ...) object

# family
`family.default` <- function (object, ...)  {
	cl <- .getCall(object)
	if(is.null(cl)) 
		return(NULL)
	fam <- cl$family
	if(is.null(fam)) 
		fam <- formals(match.fun(cl[[1L]]))$family
	if(is.null(fam))
		return(NA)
	switch(mode(fam), call = eval(fam), name =, character = match.fun(fam)())
}


`family.mer` <- function (object, ...)  {
	if(.getCall(object)[[1L]] == "lmer" && inherits(object, "mer"))
		gaussian() else family.default(object)
}

`family.gls` <-
`family.lme` <- function (object, ...) {
	if (inherits(object$family, "family")) object$family else gaussian()
}


# Classes 'coxme' and 'lmekin' from package 'coxme':


`formula.coxme` <-
function(x, ...)  {
	ret <- x$formulaList$fixed
	f <- ret[[3L]]
	for(f1 in x$formulaList$random) f <- call("+", f, f1)
	ret[[3L]] <- f
	ret
}

`formula.lmekin` <-
function(x, ...) eval(x$call$formula, parent.frame())


## Classes 'hurdle' and 'zeroinfl' from package 'pscl':

`family.zeroinfl` <-
function(object, ...) binomial(link = object$link)


#_______________________________________________________________________________

`formula.glimML` <- function(x, ...) x@formula

`family.glimML` <- function(object, ...) switch(object@method,
	"BB" = binomial(object@link),
	#"NB" = MASS::negative.binomial(theta = 1/object@param['phi.(Intercept)'],
	"NB" = get("negative.binomial", asNamespace("MASS"))(
		theta = 1 / object@param['phi.(Intercept)'], link = object@link))

#_______________________________________________________________________________

`terms.glimML` <- function (x, ...) terms.formula(x@formula, ...)

#_______________________________________________________________________________
# adds 'se.fit' argument
# https://stat.ethz.ch/pipermail/r-help/2004-April/050144.html
# http://web.archiveorange.com/archive/v/rOz2zbtjRgntPMuIDoIl

# based on the original 'predict.gls' in package 'nlme'
`predict.gls` <-
function (object, newdata, se.fit = FALSE, na.action = na.fail, ...) {
    if (missing(newdata) && missing(se.fit)) return(fitted(object))
	
    form <- formula(object)[-2L]
	if(length(form[[2L]]) == 3L && form[[2L]][[1L]] == "|" )
		form[[2L]] <- form[[2L]][[2L]] 
	
	dataMod <- model.frame(object, data = newdata, na.action = na.action,
										drop.unused.levels = TRUE)
	contr <- object$contrasts
	for(i in names(contr)) {
		levs <- levels(dataMod[, i])
		if (any(wch <- is.na(match(levs, rownames(contr[[i]]))))) {
                stop(sprintf(ngettext(sum(wch), "level %s not allowed for %s", 
                  "levels %s not allowed for %s"), paste(levs[wch], 
                  collapse = ",")), domain = NA)
            }
		attr(dataMod[[i]], "contrasts") <- contr[[i]][levs, , drop = FALSE]
	}
	X <- model.matrix(terms(form), data = dataMod)
	
    cf <- coef(object)
    val <- c(X[, names(cf), drop = FALSE] %*% cf)
	if(se.fit) {
		se <- sqrt(diag(X %*% vcov(object) %*% t(X)))
		val <- list(fit = val, se.fit = unname(se))
	}
	attr(val, "label") <- "Predicted values"
    if (!is.null(aux <- attr(object, "units")$y)) {
        attr(val, "label") <- paste(attr(val, "label"), aux)
    }
	val
}

`predict.lme` <-
function (object, newdata, level, asList = FALSE,
	na.action = na.fail, se.fit = FALSE, ...) {
	cl <- match.call()
	cl$se.fit <- NULL
	cl[[1L]] <- call("get", "predict.lme", asNamespace("nlme"))	
	res <- eval(cl, parent.frame())
	
	if(se.fit && (missing(level) || any(level > 0)))
		warning("cannot calculate standard errors for level > 0")
	if(se.fit && !missing(level) && length(level) == 1L && all(level == 0)) {
		if (missing(newdata) || is.null(newdata)) {
			X <- model.matrix(object, data = object$data)
		} else {
			tt <- delete.response(terms(formula(object)))
			xlev <- .getXlevels(tt, model.frame(object, data = object$data))
			X <- model.matrix(tt, data = newdata, contrasts.arg =
							  object$contrasts, xlev = xlev)
		}
		se <- sqrt(diag(X %*% vcov(object) %*% t(X)))
		names(se) <- names(res)
		list(fit = c(res), se.fit = se)
	} else res
}


`predict.mer` <- function (object, newdata, type = c("link", "response"), se.fit = FALSE, 
    ...)
.predict_glm(object, newdata, type, se.fit,
		trms =  delete.response(attr(object@frame, "terms")),
		coeff = object@fixef,
		offset = object@offset,
		...)

`predict.merMod` <- function (object, newdata, type = c("link", "response"), se.fit = FALSE, 
    ...)
.predict_glm(object, newdata, type, se.fit,
		trms = delete.response(terms(formula(object, fixed.only = TRUE))),
		coeff = fixef(object),
		offset = lme4::getME(object, "offset"),
		...)


.predict_glm <- function (object, newdata, type = c("link", "response"), se.fit = FALSE,
		trms, coeff, offset, ...) {
	
    type <- match.arg(type)
    if (!missing(newdata) && !is.null(newdata)) {
        xlev <- .getXlevels(trms, model.frame(trms, data = newdata))
        X <- model.matrix(trms, data = newdata, contrasts.arg = attr(model.matrix(object), 
            "contrasts"), xlev = xlev)
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(trms, "offset"))) 
            for (i in off.num) offset <- offset + eval(attr(trms, 
                "variables")[[i + 1L]], newdata)
		
		cl <- getCall(object)
        if (!is.null(cl$offset)) 
            offset <- offset + eval(cl$offset, newdata)
    } else {
        X <- model.matrix(object)
        if (!length(offset))  offset <- NULL
    }
    y <- (X %*% coeff)[, 1L]
    if (!is.null(offset)) 
        y <- y + offset
    fam <- family(object)
    if (se.fit) {
        covmat <- as.matrix(vcov(object))
        se <- sqrt(diag(X %*% covmat %*% t(X)))
        if (type == "response" && inherits(fam, "family")) 
            list(fit = fam$linkinv(y), se.fit = se * abs(fam$mu.eta(y)))
        else list(fit = y, se.fit = se)
    }
    else {
        if (type == "response" && inherits(fam, "family")) 
            fam$linkinv(y)
        else y
    }
}

# support for unmarked

#setMethod("logLik", "unmarkedFit", logLik.unmarkedFit)

`formula.unmarkedFit` <- function (x, ...) x@formula


##=============================================================================
## Classes: gee & geeglm
##=============================================================================

`coef.geese` <- 
function (object, ...) object$beta


## What if 'data' changed in the meantime?
# model.matrix.gee <-
# function (object, ...) {
	# cl <- .getCall(fgee)
	# cl[[1L]] <- as.name("model.matrix")
	# cl$object <- cl$formula
	# cl$id <- cl$corstr <- cl$formula <- NULL
	# eval(cl, parent.frame())
# }


##=============================================================================
## Class: yags
##=============================================================================

`coef.yagsResult` <-
function (object, ...)
structure(object@coefficients, names = object@varnames)


`getCall.yagsResult` <-
	function(x, ...) x@Call
	
	

`formula.yagsResult` <-
function (x, ...) 
eval(x@Call$formula, parent.frame())


##=============================================================================
## Class: MCMCglmm
##=============================================================================

`formula.MCMCglmm` <-
function (x, ...) x$Fixed$formula


`family.MCMCglmm` <-
function (object, ...) object$family


`formula.caic` <-
function(x, ...) formula(x$mod)

