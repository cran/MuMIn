# cbind list of data.frames omitting duplicated column (names)
`cbindDataFrameList` <-
function(x) {
	dfnames <- unlist(lapply(x, colnames))
	uq <- !duplicated(dfnames)
	res <- do.call("cbind", x)[,uq]
	colnames(res) <- dfnames[uq]
	return(res)
}

# same for rbind, check colnames and add NA's when any are missing
`rbindDataFrameList` <-
function(x) {
	all.colnames <- unique(unlist(lapply(x, colnames)))
	x <- lapply(x, function(y) {
		y[all.colnames[!(all.colnames %in% colnames(y))]] <- NA
		return(y[all.colnames])
	})
	return(do.call("rbind", x))
}

# test for marginality constraints
`formulaAllowed` <-
function(frm, except=NULL) {
	if(isTRUE(except)) return(TRUE)
	factors <- attr(terms(frm), "factors")
	if(length(factors) == 0) return(TRUE)
	if(is.character(except))
		factors <- factors[!(rownames(factors) %in% except), ]
	return(all(factors < 2))
}

# Calculate Akaike weights
`Weights` <-
function(aic, ...) {
	delta <- aic - min(aic)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))
	return (weight)
}

if (!existsFunction("nobs")) {

`nobs` <- function(object, ...) UseMethod("nobs")
`nobs.default` <- function(object, ...) NROW(resid(object, ...))

}

`coefDf` <- function(x) UseMethod("coefDf")
`coefDf.lme` <- function(x) x$fixDF$X
`coefDf.mer` <- function(x) rep(NA, x@dims[["p"]])
`coefDf.gls` <- function(x) rep(x$dims$N - x$dims$p, x$dims$p)
`coefDf.default` <- function(x) rep(tryCatch(df.residual(x), error=function(e) NA), length(coef(x)))

# Hidden functions

`.getLogLik` <- function()
	if ("stats4" %in% loadedNamespaces())
        stats4:::logLik else
		logLik

`.getCall` <- function(x) {
	if(mode(x) == "S4") {
		if ("call" %in% slotNames(x)) slot(x, "call") else
			NULL
	} else {
		if(!is.null(x$call)) {
			x$call
		} else if(!is.null(attr(x, "call"))) {
			attr(x, "call")
		} else
			NULL
	}
}

`.isREMLFit` <- function(x) {
	if (inherits(x, "mer")) return (x@dims[["REML"]] != 0)
	if (inherits(x, c("lme", "gls", "gam")) && !is.null(x$method))
		return(x$method %in% c("lme.REML", "REML"))
	if (any(inherits(x, c("lmer", "glmer"))))
		return(x@status["REML"] != 0)
	return(NA)
}


`.getRank` <- function(rank = NULL, rank.args = NULL, object = NULL, ...) {
	rank.args <- c(rank.args, list(...))

	if(is.null(rank)) {
		IC <- AICc
		attr(IC, "call") <- call("AICc", as.name("x"))
		return(IC)
	}
	srank <- substitute(rank, parent.frame())
	if(srank == "rank") srank <- substitute(rank)

	rank <- match.fun(rank)
	ICName <- switch(mode(srank), call=as.name("IC"), character=as.name(srank), name=, srank)
	ICarg <- c(list(as.name("x")), rank.args)
	ICCall <- as.call(c(ICName, ICarg)) 
	if(is.null(rank.args) || length(rank.args) == 0L) {
		IC <- rank
	} else {
		IC <- as.function(c(alist(x=), list(substitute(do.call("rank", ICarg), list(ICarg=ICarg)))))   
	}

	if(!is.null(object)) {
		test <- IC(object)
		if (!is.numeric(test) || length(test) != 1L)
			stop("'rank' should return numeric vector of length 1")
	}
	
	attr(IC, "call") <- ICCall
	IC
}

`matchCoef` <- function(m1, m2, all.terms = getAllTerms(m2)) {
	int <- attr(all.terms, "intercept")
	if(!is.null(int) && int != 0L) all.terms <- c("(Intercept)", all.terms)
	terms1 <- getAllTerms(m1)
	if(attr(terms1, "intercept")) terms1 <- c("(Intercept)", terms1)
	if(any((terms1 %in% all.terms) == FALSE)) stop("'m1' is not nested within 'm2")
	
	row <- structure(rep(NA, length(all.terms)), names=all.terms)
	coef1 <- coeffs(m1)
	row[terms1] <- NaN
	cf <- coef1[match(terms1, names(coef1), nomatch=0)]
	row[names(cf)]  <- cf
	row
}



#sorts alphabetically interaction components in model term names
`fixCoefNames` <-
function(x) {
	if(!is.character(x)) return(x)
	return(sapply(lapply(strsplit(x, ":"), sort), paste, collapse=":"))
}
