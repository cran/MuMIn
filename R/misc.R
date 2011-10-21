# compatibility with older versions of R

if(!("intercept" %in% names(formals(stats::reformulate)))) {
	`reformulate` <- function (termlabels, response = NULL, intercept = TRUE) {
		ret <- stats::reformulate(termlabels, response = response)
		if (!intercept) ret <- update.formula(ret, .~. -1)
		attr(ret, ".Environment") <- parent.frame()
		ret
	}
}

`DebugPrint` <- function(x) { cat(deparse(substitute(x)), "= \n") ; print(x) }


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
function(x)  UseMethod("Weights")

`Weights.model.selection` <-
function(x) x[, "weight"] / sum(x[, "weight"])

`Weights.averaging` <-
function(x) x$summary$Weight

`Weights.data.frame` <-
function(x) {
	if(ncol(x) == 2L && colnames(x)[2L] %in% c("AIC", "AICc", "BIC", "QAIC", "QAICc")
	&& is.numeric(x[, 2L]))
		Weights.default(x[, 2L])
	else NA
}

`Weights.default` <-
function(x) {
	delta <- x - min(x)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))
	return (weight)
}

if (!exists("nobs", mode = "function", where = "package:stats", inherits = FALSE)) {

`nobs` <- function(object, ...) UseMethod("nobs")
`nobs.default` <- function(object, ...) NROW(resid(object, ...))
`nobs.glm` <- function (object, ...) sum(!is.na(object$residuals))

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
		stats::logLik

`.getCall` <- function(x) {
	if(isS4(x)) {
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
		IC <- as.function(c(alist(x=, do.call("AICc", list(x)))))
		x <- NULL # just not to annoy R check
		as.function(c(alist(x=, do.call("AICc", list(x)))))
		attr(IC, "call") <- call("AICc", as.name("x"))
		class(IC) <- c("function", "ICWithCall")
		return(IC)
	} else if(inherits(rank, "ICWithCall") && length(rank.args) == 0L) {
		return(rank)
	}

	srank <- substitute(rank, parent.frame())
	if(srank == "rank") srank <- substitute(rank)

	rank <- match.fun(rank)
	ICName <- switch(mode(srank), call=as.name("IC"), character=as.name(srank), name=, srank)
	ICarg <- c(list(as.name("x")), rank.args)
	ICCall <- as.call(c(ICName, ICarg))
	IC <- as.function(c(alist(x=), list(substitute(do.call("rank", ICarg), list(ICarg=ICarg)))))

	if(!is.null(object)) {
		test <- IC(object)
		if (!is.numeric(test) || length(test) != 1L)
			stop("'rank' should return numeric vector of length 1")
	}

	attr(IC, "call") <- ICCall
	class(IC) <- c("function", "ICWithCall")
	IC
}

`matchCoef` <- function(m1, m2, all.terms = getAllTerms(m2, intercept = TRUE), beta=FALSE) {
	terms1 <- getAllTerms(m1, intercept = TRUE)
	if(any((terms1 %in% all.terms) == FALSE)) stop("'m1' is not nested within 'm2")

	row <- structure(rep(NA, length(all.terms)), names=all.terms)
	#coef1 <- coeffs(m1)
	coef1 <- if (beta) beta.weights(m1)[, 3L] else coeffs(m1)
	names(coef1) <- fixCoefNames(names(coef1))


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

#Tries to find out whether the models are fitted to the same data
.checkModels <- function(models, error = TRUE) {
	#
	cl <- sys.call(sys.parent())
	err <-  if (error) 	function(x) stop(simpleError(x, cl))
		else function(x) warning(simpleWarning(x, cl))
	res <- TRUE

	responses <- lapply(models, function(x) {
	  f <- formula(x)
	  if((length(f) == 2L) || (is.call(f[[2L]]) && f[[2L]][[1L]] == "~")) 0 else f[[2L]]
	})


 	if(!all(vapply(responses[-1L], "==", logical(1), responses[[1L]]))) {
		err("Response differs between models")
		res <- FALSE
	}


	datas <- lapply(models, function(x) .getCall(x)$data)
	# when using only 'nobs' - seems to be evaluated first outside of MuMIn namespace
	# which e.g. gives an error in glmmML - the glmmML::nobs method is faulty.
	nresid <- vapply(models, function(x) nobs(x), numeric(1L)) # , nall=TRUE

	if(!all(datas[-1L] == datas[[1]]) || !all(nresid[-1L] == nresid[[1L]])) {
		err("Models are not all fitted to the same data")
		res <- FALSE
	}
	invisible(res)
}

`videntical` <-
function(x) all(vapply(x[-1L], identical, logical(1), x[[1L]]))

# Check class for each object in a list
`linherits` <- function(x, whats) {
	as.logical(vapply(x, inherits, integer(length(whats)), names(whats),
		which=TRUE)) == whats
}


`.makeModelNames` <- function(models, all.terms) {
	allterms1 <- lapply(models, getAllTerms)
	if(missing(all.terms))	all.terms <- unique(unlist(allterms1))

	ret <- vapply(allterms1,
		function(x) paste(match(x, all.terms), collapse="+"),
		character(1L))
	if(videntical(ret)) {
		fam <- sapply(models, function(x) {
			tryCatch(unlist(family(x)[c("family", "link")]),
				error=function(e) c("", ""))
		})
		fam <- fam[!apply(fam, 1, videntical), , drop=FALSE]
		if(nrow(fam) > 0L)
			ret <- paste(ret, apply(fam, 2,  paste, collapse="."), sep="|")
	}

	cl <- lapply(models, .getCall)


	x <- lapply(cl, function(x) sapply(x[-1L], function(argval) {
		if(is.numeric(argval)) signif(argval, 3L) else deparse(argval, nlines=1)
	}))
	x <- rbindDataFrameList(lapply(lapply(x, t), as.data.frame))
	x$formula <- x$fixed <- x$model <- x$data <- NULL
	x <- as.matrix(x)
	x[is.na(x) | x == "NULL"] <- ""
	#x[x == "NULL"] <- NA

	#x <- x[, sapply(x, nlevels) != 1L, drop=FALSE]
	x <- x[, sapply(apply(x, 2, unique), length) != 1L, drop=FALSE]

	x[is.na(x)] <- ""
	if(ncol(x)) {
		ret <- paste(ret,
		gsub("([\"'\\s]+|\\w+ *=)","", apply(x, 1L, paste, collapse="/"), perl=TRUE),
		sep=":")
	}

	if(any(duplicated(ret))) {
		ret <- sprintf("Model %2$s (%1$s)", ret, format(seq_along(models)))
	}

	#ifelse(vapply(ret, is.null, boolean(1)), seq_along(ret), ret)
	ret
}
