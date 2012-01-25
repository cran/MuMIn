`dredge` <-
function(global.model, beta = FALSE, evaluate = TRUE, rank = "AICc",
		 fixed = NULL, m.max = NA, m.min = 0, subset, marg.ex = NULL,
		 trace = FALSE, varying, extra, ...) {

	# *** Rank ***
	rank.custom <- !missing(rank)
	rankArgs <- list(...)
	IC <- .getRank(rank, rankArgs)
	ICName <- as.character(attr(IC, "call")[[1L]])

	allTerms <- allTerms0 <- getAllTerms(global.model, intercept = TRUE)

	# Intercept(s)
	interceptLabel <- attr(allTerms, "interceptLabel")
	if(is.null(interceptLabel)) interceptLabel <- "(Intercept)"
	nInts <- sum(attr(allTerms, "intercept"))

	#XXX: use.ranef <- FALSE
	#if(use.ranef && inherits(global.model, "mer")) {
		#allTerms <- c(allTerms, paste("(", attr(allTerms0, "random.terms"), ")",
			#sep = ""))
	#}

	if(length(grep(":", all.vars(reformulate(allTerms))) > 0L))
		stop("variable names in the formula cannot contain \":\"")

	gmEnv <- parent.frame()
	gmCall <- .getCall(global.model)

	if (is.null(gmCall)) {
		gmCall <- substitute(global.model)
		if(!is.call(gmCall)) {
			if(inherits(global.model, c("gamm", "gamm4")))
				message("for 'gamm' models use 'MuMIn::gamm' wrapper")
			stop("could not retrieve the call to 'global.model'")
		}
		#"For objects without a 'call' component the call to the fitting function \n",
		#" must be used directly as an argument to 'dredge'.")
		# NB: this is unlikely to happen:
		if(!exists(as.character(gmCall[[1L]]), parent.frame(), mode="function"))
			 stop("could not find function '", gmCall[[1L]], "'")

	} else {
		# if 'update' method does not expand dots, we have a problem
		# with expressions like ..1, ..2 in the call.
		# So, try to replace them with respective arguments in the original call
		is.dotted <- grep("^\\.\\.", sapply(as.list(gmCall), deparse))
		if(length(is.dotted) > 0L) {
			substGmCall <- substitute(global.model)
			if(is.name(substGmCall)) {
				stop("call to 'global.model' contains '...' arguments and ",
					"cannot be updated: ", deparse(gmCall, control = NULL))
			} else gmCall[is.dotted] <-
				substitute(global.model)[names(gmCall[is.dotted])]
		}
	}
	logLik <- .getLogLik()

	# Check for na.omit
	if (!is.null(gmCall$na.action) &&
		as.character(gmCall$na.action) %in% c("na.omit", "na.exclude")) {
		stop("'global.model' should not use 'na.action' = ", gmCall$na.action)
	}

	if(names(gmCall)[2L] == "") names(gmCall)[2L] <-
		names(formals(deparse(gmCall[[1L]]))[1L])

	# TODO: other classes: model, fixed, etc...
		#gmCoefNames0 <- names(coeffs(global.model))
		#gmCoefNames <- fixCoefNames(gmCoefNames0)
	gmCoefNames <- fixCoefNames(names(coeffs(global.model)))
	#sglobCoefNames <- fixCoefNames(names(coeffs(global.model)))

	n.vars <- length(allTerms)

	if(isTRUE(rankArgs$REML) || (isTRUE(.isREMLFit(global.model)) && is.null(rankArgs$REML)))
		warning("comparing models fitted by REML")

	if (beta && is.null(tryCatch(beta.weights(global.model), error=function(e) NULL,
		warning = function(e) NULL))) {
		warning("do not know how to calculate beta weights for ",
				class(global.model)[1L], ", argument 'beta' ignored")
		beta <- FALSE
	}

	m.max <- if (missing(m.max)) (n.vars - nInts) else min(n.vars - nInts, m.max)

	# fixed variables:
	if (!is.null(fixed)) {
		if (inherits(fixed, "formula")) {
			if (fixed[[1L]] != "~" || length(fixed) != 2L)
				warning("'fixed' should be a one-sided formula")
			fixed <- c(getAllTerms(fixed))
		} else if (!is.character(fixed)) {
			stop ("'fixed' should be either a character vector with"
				  + " names of variables or a one-sided formula")
		}
		if (!all(fixed %in% allTerms)) {
			warning("not all terms in 'fixed' exist in 'global.model'")
			fixed <- fixed[fixed %in% allTerms]
		}
	}
	fixed <- c(fixed, allTerms[allTerms %in% interceptLabel])
	n.fixed <- length(fixed)
	termsOrder <- order(allTerms %in% fixed)
	allTerms <- allTerms[termsOrder]
	gmFormulaEnv <- environment(as.formula(formula(global.model), env = gmEnv))
	# TODO: gmEnv <- gmFormulaEnv ???

	### BEGIN:
	## varying BEGIN
	if(!missing(varying) && !is.null(varying)) {
		nvarying <- length(varying)
		varying.names <- names(varying)
		fvarying <- unlist(varying, recursive = FALSE)
		vlen <- vapply(varying, length, 1L)
		nvariants <- prod(vlen)
		variants <- as.matrix(expand.grid(split(seq_len(sum(vlen)),
			rep(seq_along(varying), vlen))))
	} else {
		variants <- NULL
		nvariants <- 1L
		nvarying <- 0L
		varying.names <- character(0L)
	}
	## varying END

	## extra BEGIN
	if(!missing(extra) && length(extra) != 0L) {
		extraNames <- sapply(extra, function(x) switch(mode(x),
			call = deparse(x[[1]]), name = deparse(x), character = , x))
		if(!is.null(names(extra)))
			extraNames <- ifelse(names(extra) != "", names(extra), extraNames)

		extra <- structure(as.list(unique(extra)), names = extraNames)

		if(any(c("adjR^2", "R^2") %in% extra)) {
			null.fit <- null.fit(global.model, TRUE, gmFormulaEnv)
			extra[extra == "R^2"][[1L]] <- function(x) r.squaredLR(x, null.fit)
			extra[extra == "adjR^2"][[1L]] <-
				function(x) attr(r.squaredLR(x, null.fit), "adj.r.squared")
		}

		extra <- sapply(extra, match.fun, simplify = FALSE)
		applyExtras <- function(x) unlist(lapply(extra, function(f) f(x)))
		extraResult <- applyExtras(global.model)

		if(!is.numeric(extraResult))
			stop("function in 'extra' returned non-numeric result")

		nextra <- length(extraResult)
		extraNames <- names(extraResult)
	} else {
		nextra <- 0L
		extraNames <- character(0L)
	}
	## extra END

	nov <- as.integer(n.vars - n.fixed)
	ncomb <- (2L ^ nov) * nvariants

	if(nov > 31L) stop(gettextf("maximum number of predictors is 31, but %d is given", nov))
	#if(nov > 10L) warning(gettextf("%d predictors will generate up to %.0f combinations", nov, ncomb))
	nmax <- ncomb * nvariants
	if(evaluate) {
		ret.nchunk <- 25L
		ret.ncol <- n.vars + nvarying + 3L + nextra
		ret <- matrix(NA_real_, ncol = ret.ncol, nrow = ret.nchunk)
		retCoefTable <- vector(ret.nchunk, mode = "list")
	} else {
		ret.nchunk <- nmax
	}

	calls <- vector(mode = "list", length = ret.nchunk)

	if(hasSubset <- !missing(subset))  {
		if(!tryCatch(is.language(subset), error = function(e) FALSE))
			subset <- substitute(subset)
		if(inherits(subset, "formula")) {
			if (subset[[1L]] != "~" || length(subset) != 2L)
				stop("'subset' should be a one-sided formula")
			subset <- subset[[2L]]
		}
		if(!all(all.vars(subset) %in% allTerms))
			warning("not all terms in 'subset' exist in 'global.model'")
	}

	comb.sfx <- rep(TRUE, n.fixed)
	comb.seq <- if(nov != 0L) seq_len(nov) else 0L
	k <- 0L
	extraResult1 <- integer(0L)
	ord <- integer(ret.nchunk)

	argsOptions <- list(
		response = attr(allTerms0, "response"),
		intercept = nInts,
		interceptLabel = interceptLabel,
		random = attr(allTerms0, "random"),
		gmCall = gmCall,
		gmEnv = gmEnv,
		allTerms = allTerms0,
		gmCoefNames = gmCoefNames,
		gmDataHead = if(!is.null(gmCall$data)) {
			if(eval(call("is.data.frame", gmCall$data), gmEnv))
				eval(call("head", gmCall$data, 1L), gmEnv) else gmCall$data
			} else NULL,
		gmFormulaEnv = gmFormulaEnv
		)

	retColIdx <- if(nvarying) -n.vars - seq_len(nvarying) else TRUE

	prevJComb <- 0L
	for(iComb in seq.int(ncomb)) {
		jComb <- ceiling(iComb / nvariants)
		if(jComb != prevJComb) {
			isok <- TRUE
			prevJComb <- jComb
			comb <- c(as.logical(intToBits(jComb - 1L)[comb.seq]), comb.sfx)
			nvar <- sum(comb) - nInts

			if(nvar > m.max || nvar < m.min || (hasSubset &&
				!eval(subset, structure(as.list(comb), names = allTerms)))) {
				isok <- FALSE
				next;
			}
			newArgs <- makeArgs(global.model, allTerms[comb], comb, argsOptions)
			formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
				attr(newArgs, "formulaList")
			if(!all(vapply(formulaList, formulaAllowed, logical(1L), marg.ex)))  {
				isok <- FALSE; next;
			}
			if(!is.null(attr(newArgs, "problems"))) {
				print.warnings(structure(vector(mode = "list",
					length = length(attr(newArgs, "problems"))),
						names = attr(newArgs, "problems")))
			} # end if <problems>

			cl <- gmCall
			cl[names(newArgs)] <- newArgs
		} #  end if(jComb != prevJComb)

		if(!isok) next;
		## --- Variants ---------------------------
		clVariant <- cl
		if (nvarying) clVariant[varying.names] <- 
			fvarying[variants[(iComb - 1L) %% nvariants + 1L, ]]

		if(trace) {
			cat(iComb, ": "); print(clVariant)
			utils::flush.console()
		}

		if(evaluate) {
			# begin row1: (clVariant, gmEnv, modelId, IC(), applyExtras(),
			#              nextra, allTerms, beta,
			#              if(nvarying) variantsIdx[v] else NULL
			fit1 <- tryCatch(eval(clVariant, gmEnv), error = function(err) {
				err$message <- paste(conditionMessage(err), "(model",
					iComb, "skipped)", collapse = "")
				class(err) <- c("simpleError", "warning", "condition")
				warning(err)
				return(NULL)
			})

			if (is.null(fit1)) next;

			if(nextra != 0L) {
				extraResult1 <- applyExtras(fit1)
				if(length(extraResult1) < nextra) {
					tmp <- rep(NA_real_, nextra)
					tmp[match(names(extraResult1), names(extraResult))] <- extraResult1
					extraResult1 <- tmp
				}
				#row1 <- c(row1, extraResult1)
			}

			mcoef1 <- matchCoef(fit1, all.terms = allTerms, beta = beta,
				allCoef = TRUE)

			ll <- logLik(fit1)
			row1 <- c(mcoef1[allTerms], extraResult1,
				df = attr(ll, "df"), ll = ll, ic = IC(fit1)
			)
			## end -> row1

			k <- k + 1L # all OK, add model to table

			ret.nrow <- nrow(ret)
			if(k > ret.nrow) { # append if necesarry
				nadd <- min(ret.nchunk, nmax - ret.nrow)
				retCoefTable <- c(retCoefTable, vector(nadd, mode = "list"))
				ret <- rbind(ret, matrix(NA, ncol = ret.ncol, nrow = nadd),
					deparse.level = 0L)
					calls <- c(calls, vector("list", nadd))
				ord <- c(ord, integer(nadd))
			}

			ord[k] <- iComb
			ret[k, retColIdx] <- row1
			retCoefTable[[k]] <- attr(mcoef1, "coefTable")
		} else { # if evaluate
			k <- k + 1L # all OK, add model to table
		}
		calls[[k]] <- clVariant
	} ### for (iComb ...)

	if(k == 0L) stop("the result is empty")
	names(calls) <- ord
	if(!evaluate) return(calls[seq_len(k)])

	if(k < nrow(ret)) {
		i <- seq_len(k)
		ret <- ret[i, , drop = FALSE]
		ord <- ord[i]
		calls <- calls[i]
		retCoefTable <- retCoefTable[i]
	}

	if(nvarying) {
		varlev <- ord %% nvariants; varlev[varlev == 0L] <- nvariants
		ret[, n.vars + seq_len(nvarying)] <- variants[varlev, ]
	}

	ret <- as.data.frame(ret)
	row.names(ret) <- ord

	# Convert columns with presence/absence of terms to factors
	tfac <- which(!(allTerms %in% gmCoefNames))

	ret[tfac] <- lapply(ret[tfac], factor, levels = NaN, labels="+")

	i <- seq_along(allTerms)
	v <- order(termsOrder)
	ret[, i] <- ret[, v]
	allTerms <- allTerms[v]
	colnames(ret) <- c(allTerms, varying.names, extraNames, "df", "logLik", ICName)

	if(nvarying) {
		variant.names <- lapply(varying, function(x)
			make.unique(if(is.null(names(x))) as.character(x) else names(x)))
		for (i in varying.names) ret[, i] <-
			factor(ret[, i], labels = variant.names[[i]])
	}

	o <- order(ret[, ICName], decreasing = FALSE)
	ret <- ret[o, ]
	retCoefTable <- retCoefTable[o]

	ret$delta <- ret[, ICName] - min(ret[, ICName])
	ret$weight <- exp(-ret$delta / 2) / sum(exp(-ret$delta / 2))

	structure(ret,
		class = c("model.selection", "data.frame"),
		calls = calls[o],
		global = global.model,
		global.call = gmCall,
		terms = structure(allTerms, interceptLabel = interceptLabel),
		rank = IC,
		rank.call = attr(IC, "call"),
		beta = beta,
		call = match.call(expand.dots = TRUE),
		coefTables = retCoefTable,
		nobs = nobs(global.model),
		vCols = varying.names
	)
} ######

`subset.model.selection` <-
function(x, subset, select, recalc.weights = TRUE, ...) {
	if (missing(select)) {
		if(missing(subset)) return(x)
		e <- .substHas(substitute(subset))
		i <- eval(e, x, parent.frame())
		return(`[.model.selection`(x, i, recalc.weights = recalc.weights, ...))
	} else {
		cl <- match.call(expand.dots = FALSE)
	    cl <- cl[c(1L, match(names(formals("subset.data.frame")), names(cl), 0L))]
	    cl[[1L]] <- as.name("subset.data.frame")
		ret <- eval(cl, parent.frame())
		if(recalc.weights && ("weight" %in% colnames(ret)))
			ret[, 'weight'] <- ret[, 'weight'] / sum(ret[, 'weight'])
	    return(ret)
	}
}

`[.model.selection` <-
function (x, i, j, recalc.weights = TRUE, ...) {
	ret <- `[.data.frame`(x, i, j, ...)
	if (missing(j)) {
		s <- c("row.names", "calls", "coefTables", "random.terms", "order")
		k <- match(dimnames(ret)[[1L]], dimnames(x)[[1L]])
		attrib <- attributes(x)
		attrib[s] <- lapply(attrib[s], `[`, k)
		attributes(ret) <- attrib
		if(recalc.weights)
			ret$weight <- ret$weight / sum(ret$weight)
			#ret[, 'weight'] <- ret[, 'weight'] / sum(ret[, 'weight'])

		if(!is.null(warningList <- attr(ret, "warnings")))
			attr(ret, "warnings") <- warningList[sapply(warningList, attr, "id") %in% rownames(ret)]
	} else {
		cls <- class(ret)
		class(ret) <- cls[cls != "model.selection"] # numeric or data.frame
	}
	return(ret)
}

`print.model.selection` <-
function(x, abbrev.names = TRUE, warnings = getOption("warn") != -1L, ...) {
	orig.x <- x
	if(!is.null(x$weight)) x$weight <- round(x$weight, 3L)
	xterms <- attr(x, "terms")
	if(is.null(xterms)) {
		print.data.frame(x, ...)
	} else {
		xterms <- gsub(" ", "", xterms, fixed = TRUE)
		if(abbrev.names) xterms <- abbreviateTerms(xterms, 3L)

		colnames(x)[seq_along(xterms)] <-  xterms

		globcl <- attr(x, "global.call")
		if(!is.null(globcl)) {
			cat("Global model call: ")
			print(globcl)
			cat("---\n")
			random.terms <- attr(getAllTerms(attr(x, "global")), "random.terms")
			if(!is.null(random.terms)) random.terms <- list(random.terms)
		} else random.terms <- attr(x, "random.terms")

		cat("Model selection table \n")
		dig <- c("R^2" = 4L, df = 0L, logLik = 3L, AICc = 1L, AICc = 1L,
			AIC = 1L, BIC = 1L, QAIC = 1L, QAICc = 1L, ICOMP = 1L, Cp = 1L,
			delta = 2L,	weight = 3L)

		j <- match(colnames(x), names(dig), nomatch = 0L)
		i <- sapply(x, is.numeric) & (j == 0L)
		x[, i] <- signif(x[, i], 4L)
		for(i in names(dig)[j]) x[, i] <- round(x[, i], digits = dig[i])

		vLegend <- character(0L)
		if(abbrev.names) {
			vCols <- attr(x, "vCols")
			vlen <- nchar(vCols)
			i <- vCols[1L]

			if(!is.null(vCols)) {
				for(i in vCols) {
					z <- x[, i]
					lev <- levels(z)
					lev <- lev[!(lev %in% c("", "NULL"))]
					z <- factor(z, levels = lev)
					shlev <- gsub("((?<=F)ALSE|(?<=T)RUE|~(1\\|)?| +)", "", lev,
						perl = TRUE, ignore.case = TRUE)
					shlev <- abbreviate(shlev, nchar(i))
					#shlev <- paste(substr(i, 1, 1), seq(nlevels(z)), sep="")
					x[, i] <- factor(z, labels = shlev)
					if(any(j <- shlev != lev)) vLegend <- c(vLegend, paste(i,
						": ", paste(shlev[j], "=", sQuote(lev[j]),
						collapse = ", "), sep = ""))
				}
			}
		}

		uqran <- unique(unlist(random.terms, use.names = FALSE))
		abbran <- abbreviate(gsub("1 | ", "", uqran, fixed = TRUE), 1L)
		colran <- vapply(random.terms, function(s) paste(abbran[match(s, uqran)],
			collapse = "+"), "")

		if(addrandcol <- length(unique(colran)) > 1L) {
			k <- which(colnames(x) == "df")[1L]
			x <- cbind(x[, 1L:(k - 1L)], random = colran, x[, k:ncol(x)])
		}

		print.default(as.matrix(x)[, !sapply(x, function(.x) all(is.na(.x))),
			drop = FALSE], na.print = "", quote = FALSE)

		if(abbrev.names && length(vLegend))
			cat("Abbreviations:", vLegend, sep = "\n")

		if(!is.null(random.terms)) {
			if(addrandcol) {
				cat("Random terms: \n")
				cat(paste(abbran, "=", sQuote(uqran)), sep = "\n")
			} else {
				cat("Random terms (all models): \n")
				cat(paste(sQuote(uqran)), sep = ", ")
				cat("\n")
			}

		}
		if (warnings && !is.null(attr(x, "warnings"))) {
			cat("\n"); print.warnings(attr(x, "warnings"))
		}
	}
	invisible(orig.x)
}

`update.model.selection` <- function (object, global.model, ..., evaluate = TRUE) {
    cl <- attr(object, "call")
    if (is.null(cl)) stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...

	if(!missing(global.model))
		extras <- c(list(global.model = substitute(global.model)), extras)

    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(cl)))
        for (a in names(extras)[existing]) cl[a] <- extras[a]
        if (any(!existing)) {
            cl <- c(as.list(cl), extras[!existing])
            cl <- as.call(cl)
        }
    }
    return(if (evaluate) eval(cl, parent.frame()) else cl)
}

`coef.model.selection` <- function (object, ...)
	object[, attr(object, "terms")]

`logLik.model.selection` <- function (object, ...) {
	nobs <- attr(object, "nobs")
	n <- nrow(object)
	ret <- vector(n, mode = "list")
	for(i in 1:n) ret[[i]] <-
		structure(object[i, "logLik"], df = object[i, "df"], nobs = nobs,
			class = "logLik")
	ret
}
