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

	if(length(grep(":", all.vars(reformulate(allTerms))) > 0L))
		stop("variable names in the model formula cannot contain \":\"")

	gmEnv <- parent.frame()
	gmCall <- .getCall(global.model)

	if (is.null(gmCall)) {
		gmCall <- substitute(global.model)
		if(!is.call(gmCall)) {
			if(inherits(global.model, c("gamm", "gamm4")))
				message("to use gamm models with 'dredge', use 'MuMIn::gamm' wrapper")
			stop("could not retrieve the call to 'global.model'")
		}
		#"For objects without a 'call' component the call to the fitting function \n",
		#" must be used directly as an argument to 'dredge'.")

		# this is unlikely to happen:
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
				stop("call to 'global.model' contains '...' arguments and cannot be updated: ",
					deparse(gmCall, control = NULL))
			} else {
				gmCall[is.dotted] <- substitute(global.model)[names(gmCall[is.dotted])]
			}
		}
	}

	# TODO: other classes: model, fixed, etc...
	gmFormula <- as.formula(formula(global.model))
	gmCoefNames0 <- names(coeffs(global.model))

	# Check for na.omit
	if (!is.null(gmCall$na.action) &&
		as.character(gmCall$na.action) %in% c("na.omit", "na.exclude")) {
		stop("'global.model' should not use 'na.action' = ", gmCall$na.action)
	}

	#if(names(gmCall)[2] == "") names(gmCall)[2] <- "formula"
	if(names(gmCall)[2L] == "") names(gmCall)[2L] <- names(formals(deparse(gmCall[[1]]))[1])

	gmCoefNames <- fixCoefNames(gmCoefNames0)
	#sglobCoefNames <- fixCoefNames(names(coeffs(global.model)))

	n.vars <- length(allTerms)

	if(isTRUE(rankArgs$REML) || (isTRUE(.isREMLFit(global.model)) && is.null(rankArgs$REML)))
		warning("comparing models fitted by REML")

	if (beta && is.null(tryCatch(beta.weights(global.model), error=function(e) NULL,
		warning=function(e) NULL))) {
		warning("do not know how to calculate B-weights for ",
				sQuote(class(global.model)[1L]), ", argument 'beta' ignored")
		beta <- FALSE
	}

	m.max <- if (missing(m.max)) (n.vars - nInts) else min(n.vars - nInts, m.max)

	# fixed variables:
	if (!is.null(fixed)) {
		if (inherits(fixed, "formula")) {
			if (fixed[[1]] != "~" || length(fixed) != 2L)
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
	ordAllTerms <- allTerms[termsOrder]
	#mostattributes(ordAllTerms) <- 	attributes(allTerms)
	allTerms <- ordAllTerms

	logLik <- .getLogLik()

	isMER <- any(inherits(global.model, c("mer", "lmer", "glmer")))
	gmFormulaEnv <- attr(gmFormula, ".Environment")

	### BEGIN:
	## varying BEGIN
	if(!missing(varying) && !is.null(varying)) {
		#variants <- apply(expand.grid(lapply(varying, seq_along)), 1L,
			#function(x) sapply(names(varying), function(y) varying[[y]][[x[y]]],
			#simplify=FALSE))
		variantsIdx <- expand.grid(lapply(varying, seq_along))
		#nvariants <- nrow(variantsIdx)
		seq.variants <- seq.int(nrow(variantsIdx))
		nvarying <- length(varying)
		varying.names <- names(varying)
	} else {
		variantsIdx <- NULL
		seq.variants <- 1L
		nvarying <- 0L
		varying.names <- character(0L)
	}
	nvariants <- length(seq.variants)
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
	ncomb <- 2L ^ nov

	if(nov > 31L) stop(gettextf("maximum number of predictors is 31, but %d is given", nov))
	if(nov > 10L) warning(gettextf("%d predictors will generate up to %.0f possible combinations", nov, ncomb))

	if(evaluate) {
		ret.nchunk <- 25L
		ret.ncol <- n.vars + nvarying + 3L + nextra
		ret <- matrix(NA_real_, ncol = ret.ncol, nrow = ret.nchunk)
	} else {
		ret.nchunk <- ncomb * nvariants
	}

	calls <- vector(mode = "list", length = ret.nchunk)

	if(hasSubset <- !missing(subset))  {
		if(!tryCatch(is.language(subset), error = function(e) FALSE))
			subset <- substitute(subset)
		if(inherits(subset, "formula")) {
			if (subset[[1]] != "~" || length(subset) != 2L)
				stop("'subset' should be a one-sided formula")
			subset <- subset[[2L]]
		}
		if(!all(all.vars(subset) %in% allTerms))
			warning("not all terms in 'subset' exist in 'global.model'")
	}

	comb.sfx <- rep(TRUE, n.fixed)
	comb.seq <- if(nov != 0L) seq.int(nov) else 0L
	k <- 0L
	ord <- extraResult1 <- integer(0L)

	if(!is.null(gmCall$data)) {
		if(eval(call("is.data.frame", gmCall$data), gmEnv))
			gmDataHead <- eval(call("head", gmCall$data, 1), gmEnv) else
			gmDataHead <- gmCall$data
	} else gmDataHead <- NULL

	argsOptions <- list(
		response = attr(allTerms0, "response"),
		intercept = nInts,
		interceptLabel = interceptLabel,
		random = attr(allTerms0, "random"),
		gmCall = gmCall,
		gmEnv = gmEnv,
		allTerms = allTerms0,
		gmCoefNames = gmCoefNames,
		gmDataHead = gmDataHead,
		gmFormulaEnv = gmFormulaEnv
		)

	for(j in seq.int(ncomb)) {
		comb <- c(as.logical(intToBits(j - 1L)[comb.seq]), comb.sfx)

		nvar <- sum(comb) - nInts
		if(nvar > m.max || nvar < m.min) next;
		if(hasSubset && !eval(subset, structure(as.list(comb), names=allTerms)))
			next;

		#terms1 <- allTerms[comb]
		newArgs <- makeArgs(global.model, allTerms[comb], comb, argsOptions)

		formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
			attr(newArgs, "formulaList")
		if(!all(vapply(formulaList, formulaAllowed, logical(1), marg.ex))) next;

		if(!is.null(attr(newArgs, "problems"))) {
			print.warnings(structure(vector(mode = "list",
				length = length(attr(newArgs, "problems"))),
					names = attr(newArgs, "problems")))
		}

		cl <- gmCall
		cl[names(newArgs)] <- newArgs

		for (ivar in seq.variants) { ## --- Variants ---------------------------
			clVariant <- cl
			if(nvarying) {
				newVaryingArgs <- sapply(varying.names, function(x)
					varying[[x]][[variantsIdx[ivar, x]]], simplify = FALSE)
				#for(i in varying.names) clVariant[i] <- newVaryingArgs[i]
				clVariant[varying.names] <- newVaryingArgs
			}

			modelId <- ((j - 1L) * nvariants) + ivar
			if(trace) {
				cat(modelId, ": ")
				print(clVariant)
				utils::flush.console()
			}

			if(evaluate) {
				# begin row1: (clVariant, gmEnv, modelId, IC(), applyExtras(),
				#              nextra, allTerms, beta,
				#              if(nvarying) variantsIdx[ivar] else NULL
				fit1 <- tryCatch(eval(clVariant, gmEnv), error = function(err) {
					err$message <- paste(conditionMessage(err), "(model",
						modelId, "skipped)", collapse = "")
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
				ll <- logLik(fit1)
				row1 <- c(
					matchCoef(fit1, all.terms = allTerms, beta = beta)[allTerms],
					if(nvarying) unlist(variantsIdx[ivar, ]),
					extraResult1, df = attr(ll, "df"), ll = ll, ic = IC(fit1)
				)
				## end -> row1

				k <- k + 1L # all OK, add model to table
				ord[k] <- modelId

				ret.nrow <- nrow(ret)
				if(k > ret.nrow) {
					nadd <- min(ret.nchunk, (ncomb * nvariants) - ret.nrow)
					ret <- rbind(ret, matrix(NA, ncol=ret.ncol, nrow=nadd))
					calls <- c(calls, vector("list", nadd))
				}

				ret[k, ] <- row1
			} else { # if evaluate
				k <- k + 1L # all OK, add model to table
			}
			calls[[k]] <- clVariant

		} # for (ivar ...)
	} ### for (j ...)

	names(calls) <- ord
	if(!evaluate) return(calls[seq.int(k)])

	if(k < nrow(ret)) ret <- ret[seq.int(k), , drop=FALSE]

	ret <- as.data.frame(ret)
	row.names(ret) <- ord

	# Convert columns with presence/absence of terms to factors
	tfac <- which(!(allTerms %in% gmCoefNames))

	ret[tfac] <- lapply(ret[tfac], factor, levels=NaN, labels="+")

	i <- seq_along(allTerms)
	v <- order(termsOrder)
	ret[, i] <- ret[, v]
	#ret[, seq_along(allTerms)] <- ret[, order(termsOrder)]
	allTerms <- allTerms[v]
	#colnames(ret) <- c(allTerms0, varying.names, "df", "logLik", ICName)
	colnames(ret) <- c(allTerms, varying.names, extraNames, "df", "logLik", ICName)

	if(nvarying) {
		variant.names <- lapply(varying, function(x) make.unique(if(is.null(names(x))) as.character(x) else names(x)))
		for (i in varying.names) ret[, i] <-
			factor(ret[, i], levels = seq_along(variant.names[[i]]),
				labels = variant.names[[i]])
	}

	o <- order(ret[, ICName], decreasing = FALSE)
	ret <- ret[o, ]
	ret$delta <- ret[, ICName] - min(ret[, ICName])
	ret$weight <- exp(-ret$delta / 2) / sum(exp(-ret$delta / 2))

	ret <- structure(ret,
		class = c("model.selection", "data.frame"),
		calls = calls[o],
		global = global.model,
		global.call = gmCall,
		terms = allTerms,
		rank = IC,
		rank.call = attr(IC,"call"),
		call = match.call(expand.dots = TRUE)
	)

	if (!is.null(attr(allTerms0, "random.terms")))
		attr(ret, "random.terms") <- attr(allTerms0, "random.terms")

	return(ret)
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
		s <- c("row.names", "calls")
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
	if(!is.null(x$weight))
		x$weight <- round(x$weight, 3L)
	xterms <- attr(x, "terms")
	if(is.null(xterms)) {
		print.data.frame(x, ...)
	} else {
		xterms <- gsub(" ", "", xterms)
		if(abbrev.names) xterms <- abbreviateTerms(xterms, 3L)

		colnames(x)[seq_along(xterms)] <-  xterms
		cl <- attr(x, "global.call")
		if(!is.null(cl)) {
			cat("Global model call: ")
			print(cl)
			cat("---\n")
		}

		cat("Model selection table \n")
		dig <- c("R^2" = 4L, df = 0, logLik = 3, AICc = 1, AICc = 1, AIC = 1,
			BIC = 1, QAIC = 1, QAICc = 1, ICOMP = 1, Cp = 1, delta = 2L, weight = 3L)

		j <- match(colnames(x), names(dig), nomatch = 0)
		i <- sapply(x, is.numeric) & (j == 0L)
		x[, i] <- signif(x[, i], 4L)
		for(i in names(dig)[j]) x[, i] <- round(x[, i], digits = dig[i])

		print.default(as.matrix(x)[, !sapply(x, function(.x) all(is.na(.x))),
			drop = FALSE], na.print = "", quote = FALSE)
		if (!is.null(attr(x, "random.terms"))) {
			cat("Random terms:", paste(attr(x, "random.terms"), collapse=", "),
				"\n")
		}
		if (warnings && !is.null(attr(x, "warnings"))) {
			cat("\n"); print.warnings(attr(x, "warnings"))
		}
	}
}

`update.model.selection` <- function (object, global.model, ..., evaluate = TRUE) {
    cl <- attr(object, "call")
    if (is.null(cl))
        stop("need an object with call component")
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
