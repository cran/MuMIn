`pdredge` <-
function(global.model, cluster = NA, beta = FALSE, evaluate = TRUE,
	rank = "AICc", fixed = NULL, m.max = NA, m.min = 0, subset, marg.ex = NULL,
	trace = FALSE, varying, extra, check = FALSE,  ...) {
#FIXME: m.max cannot be 0 - e.g. for intercept only model

	qlen <- 25L
	# Imports: clusterCall, clusterApply
	doParallel <- inherits(cluster, "cluster")
	if(doParallel) {
		.parallelPkgCheck() # XXX: workaround to avoid importing from 'parallel'
		clusterCall <- get("clusterCall")
		clusterApply <- get("clusterApply")
		clusterCall(cluster, "require", "MuMIn", character.only = TRUE)

		clusterCall(cluster, assign, "assignFromNs", function(name, asName = name,
			ns = "MuMIn") {
			assign(asName, get(name, loadNamespace(ns)), envir = .GlobalEnv)
			invisible(NULL)
		}, envir = .GlobalEnv)

		clusterCall(cluster, "assignFromNs", ".getLogLik")
		clusterCall(cluster, "assignFromNs", "tryCatchWE")
		clusterCall(cluster, "assignFromNs", "matchCoef")
		clusterCall(cluster, "assignFromNs", "parGetMsRow")

		.getRow <- function(X) clusterApply(cluster, X, fun = "parGetMsRow")

	} else {
		.getRow <- function(X) lapply(X, parGetMsRow, parCommonProps)
		clusterCall <- function(...) NULL
		message("Not using cluster.")
	}



	gmEnv <- parent.frame()
	gmCall <- .getCall(global.model)
	gmNobs <- nobs(global.model)

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

	# *** Rank ***
	rank.custom <- !missing(rank)
	rankArgs <- list(...)
	IC <- .getRank(rank, rankArgs)
	ICName <- as.character(attr(IC, "call")[[1L]])

	allTerms <- allTerms0 <- getAllTerms(global.model, intercept = TRUE,
		data = eval(gmCall$data, envir = gmEnv))

	# Intercept(s)
	interceptLabel <- attr(allTerms, "interceptLabel")
	if(is.null(interceptLabel)) interceptLabel <- "(Intercept)"
	nInts <- sum(attr(allTerms, "intercept"))

	if(length(grep(":", all.vars(reformulate(allTerms))) > 0L))
		stop("variable names in the formula cannot contain \":\"")
	logLik <- .getLogLik()

	# parallel: check whether the models would be identical:
	if(doParallel) testUpdatedObj(cluster, global.model, gmCall, do.eval = check)

	# Check for na.omit
	if (!is.null(gmCall$na.action) &&
		as.character(gmCall$na.action) %in% c("na.omit", "na.exclude")) {
		stop("'global.model' should not use 'na.action' = ", gmCall$na.action)
	}

	if(names(gmCall)[2L] == "") names(gmCall)[2L] <-
		names(formals(deparse(gmCall[[1L]]))[1L])

	gmCoefNames <- fixCoefNames(names(coeffs(global.model)))

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
			fixed <- as.vector(getAllTerms(fixed))
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
		variants <- varying.names <- NULL
		nvariants <- 1L
		nvarying <- 0L
	}

	## varying END

	## extra BEGIN
	if(!missing(extra) && length(extra) != 0L) {
		extraNames <- sapply(extra, function(x) switch(mode(x),
			call = deparse(x[[1L]]), name = deparse(x), character = , x))
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

	if(nov > 31L) stop(gettextf("number of predictors (%d) exceeds allowed maximum (31)"), nov, domain = "MuMIn")
	#if(nov > 10L) warning(gettextf("%d predictors will generate up to %.0f combinations", nov, ncomb))
	nmax <- ncomb * nvariants
	if(evaluate) {
		ret.nchunk <- 25L
		ret.ncol <- n.vars + nvarying + 3L + nextra
		ret <- matrix(NA_real_, ncol = ret.ncol, nrow = ret.nchunk)
		coefTables <- vector(ret.nchunk, mode = "list")
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

	# TODO: allow for 'marg.ex' per formula in multi-formula models
	if(missing(marg.ex) || (!is.null(marg.ex) && is.na(marg.ex))) {
		newArgs <- makeArgs(global.model, allTerms, rep(TRUE, length(allTerms)),
							argsOptions)
		formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs
			else attr(newArgs, "formulaList")

		marg.ex <- unique(unlist(lapply(sapply(formulaList, formulaAllowed,
			simplify = FALSE), attr, "marg.ex")))
		if(!length(marg.ex)) marg.ex <- NULL
		#cat("Marginality exceptions:", marg.ex, "\n")
	}
	###

	# BEGIN parallel
	qi <- 0L
	queued <- vector(qlen, mode = "list")
	parCommonProps <- list(gmEnv = gmEnv, IC = IC, beta = beta,
		allTerms = allTerms, nextra = nextra)
	if(nextra) {
		parCommonProps$applyExtras <- applyExtras
		parCommonProps$extraResult <- extraResult
	}
	if(doParallel) clusterVExport(cluster, clustDredgeProps = parCommonProps,
		tryCatchWE)
	# END parallel

	retColIdx <- if(nvarying) -n.vars - seq_len(nvarying) else TRUE

	warningList <- list()
	# qlen <- 4 ## DEBUG: !!!!

	prevJComb <- 0L
	for(iComb in seq.int(ncomb)) {
		jComb <- ceiling(iComb / nvariants)
		if(jComb != prevJComb) {
			isok <- TRUE
			prevJComb <- jComb
			comb <- c(as.logical(intToBits(jComb - 1L)[comb.seq]), comb.sfx)
			nvar <- sum(comb) - nInts
			if(!(nvar > m.max || nvar < m.min) && (!hasSubset || eval(subset,
				structure(as.list(comb), names = allTerms)))) {
				newArgs <- makeArgs(global.model, allTerms[comb], comb, argsOptions)
				formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
					attr(newArgs, "formulaList")
				if(all(vapply(formulaList, formulaAllowed, logical(1L), marg.ex))) {
					if(!is.null(attr(newArgs, "problems"))) {
						print.warnings(structure(vector(mode = "list",
							length = length(attr(newArgs, "problems"))),
								names = attr(newArgs, "problems")))
					} # end if <problems>
					cl <- gmCall
					cl[names(newArgs)] <- newArgs
				} else isok <- FALSE # end if <formulaAllowed>
			} else isok <- FALSE # end if <subset, m.max >= nvar >= m.min>
		} #  end if(jComb != prevJComb)

		if(isok) {
			## --- Variants ---------------------------
			clVariant <- cl
			if(nvarying) clVariant[varying.names] <-
				fvarying[variants[(iComb - 1L) %% nvariants + 1L, ]]

			if(trace) {
				cat(iComb, ": "); print(clVariant)
				utils::flush.console()
			}
			if(evaluate) {
				qi <- qi + 1L
				queued[[(qi)]] <- list(call = clVariant, id = iComb)
			} else { # if !evaluate
				k <- k + 1L # all OK, add model to table
				calls[[k]] <- clVariant
			}
		} # if isok

		#if(evaluate && qi && (qi + nvariants > qlen || iComb == ncomb)) {
		if(evaluate && qi && (qi > qlen || iComb == ncomb)) {
			#DebugPrint(paste(qi, nvariants, qlen, iComb, ncomb))
			qseq <- seq_len(qi)
			qresult <- .getRow(queued[qseq])
			#cat(sprintf("queue done: %d\n", qi)) # DEBUG

			if(any(vapply(qresult, is.null, TRUE)))
				stop("some results returned from cluster node(s) are NULL. \n",
					"This should not happen and indicates problems with ",
					"the cluster node", domain = "R-MuMIn")
			haveProblems <- logical(qi)

			nadd <- sum(sapply(qresult, function(x) inherits(x$value, "condition")
				+ length(x$warnings)))
			wi <- length(warningList)
			if(nadd) warningList <- c(warningList, vector(nadd, mode = "list"))

			# DEBUG: print(sprintf("Added %d warnings, now is %d", nadd, length(warningList)))

			for (i in qseq)
				for(cond in c(qresult[[i]]$warnings,
					if(inherits(qresult[[i]]$value, "condition"))
						list(qresult[[i]]$value))) {
						wi <- wi + 1L
						warningList[[wi]] <- if(is.null(conditionCall(cond)))
							queued[[i]]$call else conditionCall(cond)
						if(inherits(cond, "error")) {
							haveProblems[i] <- TRUE
							msgsfx <- "(model %d skipped)"
						} else
							msgsfx <- "(in model %d)"
						names(warningList)[wi] <- paste(conditionMessage(cond),
							 gettextf(msgsfx, queued[[i]]$id))
						attr(warningList[[wi]], "id") <- queued[[i]]$id
				}

			withoutProblems <- which(!haveProblems)
			qrows <- lapply(qresult[withoutProblems], "[[", "value")
			qresultLen <- length(qrows)
			retNrow <- nrow(ret)
			if(k + qresultLen > retNrow) {
				nadd <- min(ret.nchunk, nmax - retNrow)
				coefTables <- c(coefTables, vector(nadd, mode = "list"))
				ret <- rbind(ret, matrix(NA, ncol = ret.ncol, nrow = nadd),
					deparse.level = 0L)
				calls <- c(calls, vector("list", nadd))
				ord <- c(ord, integer(nadd))
			}
			qseqOK <- seq_len(qresultLen)
			for(m in qseqOK) ret[k + m, retColIdx] <- qrows[[m]]
			ord[k + qseqOK] <- vapply(queued[withoutProblems], "[[", 1L, "id")
			calls[k + qseqOK] <- lapply(queued[withoutProblems], "[[", "call")
			coefTables[k + qseqOK] <- lapply(qresult[withoutProblems], "[[", "coefTable")
			k <- k + qresultLen
			qi <- 0L
		}
	} ### for (iComb ...)

	if(k == 0L) {
		if(length(warningList)) print.warnings(warningList)
		stop("the result is empty")
	}

	names(calls) <- ord
	if(!evaluate) return(calls[seq_len(k)])

	if(k < nrow(ret)) {
		i <- seq_len(k)
		ret <- ret[i, , drop = FALSE]
		ord <- ord[i]
		calls <- calls[i]
		coefTables <- coefTables[i]
	}

	if(nvarying) {
		varlev <- ord %% nvariants; varlev[varlev == 0L] <- nvariants
		ret[, n.vars + seq_len(nvarying)] <- variants[varlev, ]
	}

	ret <- as.data.frame(ret)
	row.names(ret) <- ord

	# Convert columns with presence/absence of terms to factors
	tfac <- which(!(allTerms %in% gmCoefNames))
	ret[tfac] <- lapply(ret[tfac], factor, levels = NaN, labels = "+")

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
	coefTables <- coefTables[o]

	ret$delta <- ret[, ICName] - min(ret[, ICName])
	ret$weight <- exp(-ret$delta / 2) / sum(exp(-ret$delta / 2))

	ret <- structure(ret,
		class = c("model.selection", "data.frame"),
		calls = calls[o],
		global = global.model,
		global.call = gmCall,
		terms = structure(allTerms, interceptLabel = interceptLabel),
		rank = IC,
		rank.call = attr(IC, "call"),
		beta = beta,
		call = match.call(expand.dots = TRUE),
		coefTables = coefTables,
		nobs = gmNobs,
		vCols = varying.names
	)

	if(length(warningList)) {
		class(warningList) <- c("warnings", "list")
		attr(ret, "warnings") <- warningList
	}

	if (!is.null(attr(allTerms0, "random.terms")))
		attr(ret, "random.terms") <- attr(allTerms0, "random.terms")

	if(doParallel) clusterCall(cluster, "rm",
		list = c("parGetMsRow", "clustDredgeProps", ".getLogLik", "tryCatchWE",
			"matchCoef", "parGetMsRow"), envir = .GlobalEnv)
	return(ret)
} ######


`parGetMsRow` <- function(modv, Z = get("clustDredgeProps", .GlobalEnv)) {
	### modv == list(call = clVariant, id = modelId)
	result <- tryCatchWE(eval(modv$call, Z$gmEnv))
	if (inherits(result$value, "condition")) return(result)

	fit1 <- result$value
	if(Z$nextra != 0L) {
		extraResult1 <- Z$applyExtras(fit1)
		if(length(extraResult1) < Z$nextra) {
			tmp <- rep(NA_real_, Z$nextra)
			tmp[match(names(extraResult1), names(Z$extraResult))] <-
				extraResult1
			extraResult1 <- tmp
		}
	} else extraResult1 <- NULL
	ll <- .getLogLik()(fit1)
	mcoef <- matchCoef(fit1, all.terms = Z$allTerms, beta = Z$beta, allCoef = TRUE)

	list(value = c(mcoef, extraResult1, df = attr(ll, "df"), ll = ll,
		ic = Z$IC(fit1)),
		nobs = nobs(fit1),
		coefTable = attr(mcoef, "coefTable"),
		warnings = result$warnings)
}

.test_pdredge <- function(dd) {
	cl <- attr(dd, "call")
	cl$cluster <- cl$check <- NULL
	cl[[1]] <- as.name("dredge")
	if(!identical(c(dd), c(eval(cl)))) stop("Whoops...")
	dd
}
