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

`matchCoef` <- function(m1, m2, all.terms = getAllTerms(m2, intercept = TRUE),
	beta = FALSE) {
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
		err("response differs between models")
		res <- FALSE
	}


	datas <- lapply(models, function(x) .getCall(x)$data)
	# when using only 'nobs' - seems to be evaluated first outside of MuMIn namespace
	# which e.g. gives an error in glmmML - the glmmML::nobs method is faulty.
	nresid <- vapply(models, function(x) nobs(x), numeric(1L)) # , nall=TRUE

	if(!all(datas[-1L] == datas[[1]]) || !all(nresid[-1L] == nresid[[1L]])) {
		err("models are not all fitted to the same data")
		res <- FALSE
	}
	invisible(res)
}

`abbreviateTerms` <- function(x, n = 1L, capwords = FALSE) {
	ret <- x
	allVars <- all.vars(reformulate(x))
	pat <- paste("\\b", allVars, "\\b", sep="")
	if(capwords) {
		abx <- abbreviate(paste(toupper(substring(allVars, 1L, 1L)), 
			tolower(substring(allVars, 2L)), sep=""), n)
	} else {
		abx <- abbreviate(allVars, n)
	}
	
	for(i in seq_along(allVars)) ret <- gsub(pat[i], abx[i], ret, perl = TRUE)
	ret <- gsub("I\\((\\w+)\\)", "\\1", ret, perl = TRUE)
	attr(ret, "variables") <- structure(allVars, names = abx)
	ret
}

## do.call("substitute", list(f, sapply(abbreviate(all.names(f,  unique=T), 1), as.name)))

#models <- list(model1, model2)

`modelNames` <- function(models = NULL, allTerms, uqTerms, ...) {
	if(missing(allTerms)) 	allTerms <- lapply(models, getAllTerms)
	if(missing(uqTerms)) 	uqTerms <- unique(unlist(allTerms))
	
	ret <- sapply(allTerms, function(x) paste(sort(match(x, uqTerms)), collapse=""))

	dup <- table(ret)
	dup <- dup[dup > 1]
	
	if(length(dup) > 0) {
		idup <- which(ret %in% names(dup))

		ret[idup] <- sapply(idup, function(i) paste(ret[i], letters[sum(ret[seq.int(i)] == ret[i])], sep=""))
	}
	attr(ret, "variables") <- structure(seq_along(uqTerms), names = uqTerms)
	ret
}



`modelNames0` <- function(models, strict = FALSE, asNumeric = FALSE,
	withRandomTerms = TRUE, withFamily = TRUE, withArguments = TRUE,
	fmt = "Model %s %s"
	) {
	
	# sapply(tt, function(x) paste(sort(match(allt, x)), collapse=""))

	if(withRandomTerms) {
		allTermsList <- lapply(models, function(x) {
			tt <- getAllTerms(x)
			rtt <- attr(tt, "random.terms")
			c(tt, if(!is.null(rtt)) paste("(", rtt, ")", sep="") else NULL)
		})
	} else {
		allTermsList <- lapply(models, getAllTerms)
	}

	allTerms <- unique(unlist(allTermsList))

	if(asNumeric) {
		abvtt <- seq_along(allTerms)
	} else {
		abvtt <- abbreviateTerms(allTerms)
		variables <- attr(abvtt, "variables")
		if(withRandomTerms)
			abvtt <- gsub("\\(1 \\| (\\S+)(?: %in%.*)?\\)", "(\\1)", abvtt, perl = TRUE)
	}

	ret <- sapply(allTermsList, function(x) paste(abvtt[match(x, allTerms)], collapse="+"))

	if(withFamily && (strict || any(duplicated(ret)))) {
		fam <- sapply(models, function(x) {
					tryCatch(unlist(family(x)[c("family", "link")]),
						error=function(e) c("", ""))
				})
		fam <- paste(abbreviate(fam[1,], 4, strict=T), "(", abbreviate(fam[2,], 1, strict=FALSE), ")", sep="")
		ret <- paste(ret, fam, sep=" ")
	}

	if(withArguments && (strict || any(duplicated(ret)))) {
		cl <- lapply(models, .getCall)
		x <- lapply(cl, function(x) sapply(x[-1L], function(argval) {
			if(is.numeric(argval)) signif(argval, 3L) else deparse(argval, nlines=1)
		}))
		x <- rbindDataFrameList(lapply(lapply(x, t), as.data.frame))
		x$formula <- x$fixed <- x$model <- x$data <- x$family <- x$cluster <- NULL
		x <- as.matrix(x)
		x[is.na(x) | x == "NULL"] <- ""
		x <- as.data.frame(x)
		#x <- x[, sapply(apply(x, 2, unique), length) != 1L, drop=FALSE]
		x <- x[, sapply(lapply(x, unique), length) != 1L, drop=FALSE]

		if(ncol(x)) {
			ret <- paste(ret,
			gsub("([\"'\\s]+|\\w+ *=)","", apply(x, 1L, paste, collapse="/"), perl=TRUE),
			sep="'")
		}
	}
	if(strict || any(duplicated(ret))) {
		ret <- sprintf(fmt, format(seq_along(models)), ret)
	}

	if(!asNumeric) {
		attr(ret, "variables") <- variables
	}
	ret
}