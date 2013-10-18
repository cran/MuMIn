`pdredge` <-
function(global.model, cluster = NA, beta = FALSE, evaluate = TRUE,
	rank = "AICc", fixed = NULL, m.max = NA, m.min = 0, subset, marg.ex = NULL,
	trace = FALSE, varying, extra, ct.args = NULL, check = FALSE, ...) {
#FIXME: m.max cannot be 0 - e.g. for intercept only model

	qlen <- 25L
	# Imports: clusterCall, clusterApply
	doParallel <- inherits(cluster, "cluster")
	if(doParallel) {
		.parallelPkgCheck() # XXX: workaround to avoid importing from 'parallel'
		clusterCall <- get("clusterCall")
		clusterApply <- get("clusterApply")
		clusterCall(cluster, "require", "MuMIn", character.only = TRUE)
		.getRow <- function(X) clusterApply(cluster, X, fun = ".pdredge_process_model")
	} else {
		.getRow <- function(X) lapply(X, pdredge_process_model, envir = props)
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
		if(!exists(as.character(gmCall[[1L]]), parent.frame(), mode = "function"))
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
		## object from 'run.mark.model' has $call of 'make.mark.model' - fixing it here:
		if(inherits(global.model, "mark") && gmCall[[1L]] == "make.mark.model") {
			gmCall <- call("run.mark.model", model = gmCall, invisible = TRUE)
		}
	}

    LL <- .getLik(global.model)
	logLik <- LL$logLik
	lLName <- LL$name
	# *** Rank ***
	rank.custom <- !missing(rank)
	if(!rank.custom && lLName == "qLik") {
		rank <- "QIC"
		warning("using 'QIC' instead of 'AICc'")
	}
	rankArgs <- list(...)
	IC <- .getRank(rank, rankArgs)
	ICName <- as.character(attr(IC, "call")[[1L]])

	allTerms <- allTerms0 <- getAllTerms(global.model, intercept = TRUE,
		data = eval(gmCall$data, envir = gmEnv))

	# Intercept(s)
	interceptLabel <- attr(allTerms, "interceptLabel")
	if(is.null(interceptLabel)) interceptLabel <- "(Intercept)"
	nInts <- sum(attr(allTerms, "intercept"))


	


	# parallel: check whether the models would be identical:
	if(doParallel && check) testUpdatedObj(cluster, global.model, gmCall, level = check)

	# Check for na.omit
	if (!is.null(gmCall$na.action) &&
		as.character(gmCall$na.action) %in% c("na.omit", "na.exclude")) {
		stop("'global.model' should not use 'na.action' = ", gmCall$na.action)
	}

	if(names(gmCall)[2L] == "") gmCall <-
		match.call(gmCall, definition = eval(gmCall[[1L]], envir = parent.frame()), expand.dots = TRUE)


	gmCoefNames <- fixCoefNames(names(coeffs(global.model)))

	n.vars <- length(allTerms)

	if(isTRUE(rankArgs$REML) || (isTRUE(.isREMLFit(global.model)) && is.null(rankArgs$REML)))
		warning("comparing models fitted by REML")

	if (beta && is.null(tryCatch(beta.weights(global.model), error = function(e) NULL,
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
		} else if (identical(fixed, TRUE)) {
			fixed <- as.vector(allTerms[!(allTerms %in% interceptLabel)])
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

	### BEGIN Manage 'varying'
	## @param:	varying
	## @value:	varying, varying.names, variants, nvariants, nvarying
	if(!missing(varying) && !is.null(varying)) {
		nvarying <- length(varying)
		varying.names <- names(varying)
		fvarying <- unlist(varying, recursive = FALSE, use.names = FALSE)
		vlen <- vapply(varying, length, 1L)
		nvariants <- prod(vlen)
		variants <- as.matrix(expand.grid(split(seq_len(sum(vlen)),
			rep(seq_len(nvarying), vlen))))
		
		flat.variant.Vvals <- unlist(lapply(varying, .makeListNames),
			recursive = FALSE, use.names = FALSE)
		
	} else {
		variants <- varying.names <- NULL
		nvariants <- 1L
		nvarying <- 0L
	}
	## END: varying

	## BEGIN Manage 'extra'
	## @param:	extra, global.model, gmFormulaEnv, 
	## @value:	extra, nextra, extraNames, null.fit
	if(!missing(extra) && length(extra) != 0L) {
		extraNames <- sapply(extra, function(x) switch(mode(x),
			call = deparse(x[[1L]]), name = deparse(x), character = , x))
		if(!is.null(names(extra)))
			extraNames <- ifelse(names(extra) != "", names(extra), extraNames)
		extra <- structure(as.list(unique(extra)), names = extraNames)

		if(any(c("adjR^2", "R^2") %in% extra)) {
			null.fit <- null.fit(global.model, evaluate = TRUE, envir = gmFormulaEnv)
			extra[extra == "R^2"][[1L]] <- function(x) r.squaredLR(x, null = null.fit)
			extra[extra == "adjR^2"][[1L]] <-
				function(x) attr(r.squaredLR(x, null = null.fit), "adj.r.squared")
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
	## END: manage 'extra'

	nov <- as.integer(n.vars - n.fixed)
	ncomb <- (2L ^ nov) * nvariants

	if(nov > 31L) stop(gettextf("number of predictors (%d) exceeds allowed maximum (31)",
								nov, domain = "R-MuMIn"))
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

	## BEGIN: Manage 'subset'
	## @param:	hasSubset, subset, allTerms, [interceptLabel], 
	## @value:	hasSubset, subset
	if(missing(subset))  {
		hasSubset <- 1L
	} else {
		if(!tryCatch(is.language(subset) || is.matrix(subset), error = function(e) FALSE))
			subset <- substitute(subset)

		if(is.matrix(subset)) {
			dn <- dimnames(subset)
			#at <- allTerms[!(allTerms %in% interceptLabel)]
			n <- length(allTerms)
			if(is.null(dn) || any(sapply(dn, is.null))) {
				di <- dim(subset)
				if(any(di != n)) stop("unnamed 'subset' matrix does not have ",
					"both dimensions equal to number of terms in 'global.model': %d", n)
				dimnames(subset) <- list(allTerms, allTerms)
			} else {
				if(!all(unique(unlist(dn)) %in% allTerms))
					warning("at least some dimnames of 'subset' matrix do not ",
					"match term names in 'global.model'")
				subset <- matrix(subset[match(allTerms, rownames(subset)),
					match(allTerms, colnames(subset))],
					dimnames = list(allTerms, allTerms),
					nrow = n, ncol = n)
			}
			if(any(!is.na(subset[!lower.tri(subset)]))) {
				warning("non-missing values exist outside the lower triangle of 'subset'")
				subset[!lower.tri(subset)] <- NA
			}
			mode(subset) <- "logical"
			hasSubset <- 2L # subset as matrix
		} else {
			if(inherits(subset, "formula")) {
				if (subset[[1L]] != "~" || length(subset) != 2L)
					stop("'subset' formula should be one-sided")
				subset <- subset[[2L]]
			}
			subset <- as.expression(subset)
			ssValidNames <- c("comb", "*nvar*")
			subsetExpr <- .subst4Vec(subset[[1L]], allTerms, as.name("comb"))

			if(nvarying) {
			ssValidNames <- c("cVar", "comb", "*nvar*")
			subsetExpr <- .substFun4Fun(subsetExpr, "V", function(x, cVar, fn) {
				if(length(x) > 2L)
					warning("discarding extra arguments for 'V' in 'subset' expression")
				i <- which(fn == x[[2L]])[1L]
					if(is.na(i)) stop(sQuote(x[[2L]]), " is not a valid name of 'varying' element")
				call("[[", cVar, i)
			}, as.name("cVar"), varying.names)
			if(!all(all.vars(subsetExpr) %in% ssValidNames))
				subsetExpr <- .subst4Vec(subsetExpr, varying.names, as.name("cVar"), fun = "[[")
			}
			ssVars <- all.vars(subsetExpr)
			okVars <- ssVars %in% ssValidNames
			if(!all(okVars)) stop("unrecognized names in 'subset' expression: ",
				prettyEnumStr(ssVars[!okVars]))
			
			ssEnv <- new.env(parent = .GlobalEnv)
			ssFunc <- setdiff(all.vars(subsetExpr, functions = TRUE), ssVars)
			if("dc" %in% ssFunc) assign("dc", .subset_dc, ssEnv)
			
			hasSubset <- if(any(ssVars == "cVar")) 4L else # subset as expression
				3L # subset as expression using 'varying' variables

		}
	} # END: manage 'subset'

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
		formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
			attr(newArgs, "formulaList")

		marg.ex <- unique(unlist(lapply(sapply(formulaList, formulaMargChk,
			simplify = FALSE), attr, "marg.ex")))
		if(!length(marg.ex)) marg.ex <- NULL else
			cat("Marginality exceptions:", sQuote(marg.ex), "\n")
	}
	###

	# BEGIN parallel
	qi <- 0L
	queued <- vector(qlen, mode = "list")
	props <- list(
				gmEnv = gmEnv,
				IC = IC, 
				# beta = beta,
				# allTerms = allTerms, 
				nextra = nextra,
				matchCoefCall = as.call(c(list(
					as.name("matchCoef"), as.name("fit1"), 
					all.terms = allTerms, beta = beta, 
					allCoef = TRUE), ct.args))
				# matchCoefCall = as.call(c(alist(matchCoef, fit1, all.terms = Z$allTerms, 
				#   beta = Z$beta, allCoef = TRUE), ct.args))
		)
	if(nextra) {
		props$applyExtras <- applyExtras
		props$extraResultNames <- names(extraResult)
	}
	props <- as.environment(props)
	
	if(doParallel) clusterVExport(cluster,
								  pdredge_props = props,
								  .pdredge_process_model = pdredge_process_model
								  )
	# END parallel

	retColIdx <- if(nvarying) -n.vars - seq_len(nvarying) else TRUE

	warningList <- list()

	####
	prevJComb <- 0L
	for(iComb in seq.int(ncomb)) {
		jComb <- ceiling(iComb / nvariants)
		if(jComb != prevJComb) {
			isok <- TRUE
			prevJComb <- jComb
			comb <- c(as.logical(intToBits(jComb - 1L)[comb.seq]), comb.sfx)
			nvar <- sum(comb) - nInts
			
			
			# !!! POSITIVE condition for 'pdredge', NEGATIVE for 'dredge':
			if((nvar >= m.min && nvar <= m.max) && switch(hasSubset,
					# 1 - no subset, 2 - matrix, 3 - expression
					TRUE,                                    # 1 
					all(subset[comb, comb], na.rm = TRUE),   # 2
					.evalExprIn(subsetExpr, env = ssEnv, enclos = parent.frame(),
						comb = comb, `*nvar*` = nvar),		 # 3
					TRUE
					)
				) {

				newArgs <- makeArgs(global.model, allTerms[comb], comb, argsOptions)
				formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
					attr(newArgs, "formulaList")
				if(all(vapply(formulaList, formulaMargChk, logical(1L), marg.ex))) {
					if(!is.null(attr(newArgs, "problems"))) {
						print.warnings(structure(vector(mode = "list",
							length = length(attr(newArgs, "problems"))),
								names = attr(newArgs, "problems")))
					} # end if <problems>
					cl <- gmCall
					cl[names(newArgs)] <- newArgs
				} else isok <- FALSE # end if <formulaMargChk>
			} else isok <- FALSE # end if <subset, m.max >= nvar >= m.min>
		} #  end if(jComb != prevJComb)

		if(isok) {
			## --- Variants ---------------------------
			clVariant <- cl
			isok2 <- TRUE
			if(nvarying) {
				cvi <- variants[(iComb - 1L) %% nvariants + 1L, ]				
				isok2 <- (hasSubset != 4L) || .evalExprIn(subsetExpr, env = ssEnv,
					enclos = parent.frame(), comb = comb, `*nvar*` = nvar,
					cVar = flat.variant.Vvals[cvi])
				clVariant[varying.names] <- fvarying[cvi]
			}
			
			if(isok2) {
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
					ord[k] <- iComb
				}
			}
		} # if isok

		#if(evaluate && qi && (qi + nvariants > qlen || iComb == ncomb)) {
		if(evaluate && qi && (qi > qlen || iComb == ncomb)) {
			#DebugPrint(paste(qi, nvariants, qlen, iComb, ncomb))
			qseq <- seq_len(qi)
			qresult <- .getRow(queued[qseq])
			utils::flush.console()
		
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
				nadd <- min(max(ret.nchunk, qresultLen), nmax - retNrow)
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
	colnames(ret) <- c(allTerms, varying.names, extraNames, "df", lLName, ICName)

	if(nvarying) {
		variant.names <- vapply(flat.variant.Vvals, function(x) if(is.character(x)) x else
			deparse(x, control = NULL, width.cutoff = 20L)[1L], character(1L))
		vnum <- split(seq_len(sum(vlen)), rep(seq_len(nvarying), vlen))
		names(vnum) <- varying.names
		for (i in varying.names) ret[, i] <-
			factor(ret[, i], levels = vnum[[i]], labels = variant.names[vnum[[i]]])

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
		list = c(".pdredge_process_model", "pdredge_props"), envir = .GlobalEnv)
		

	return(ret)
} ######



`pdredge_process_model` <- function(modv, envir = get("pdredge_props", .GlobalEnv)) {
	### modv == list(call = clVariant, id = modelId)
	result <- tryCatchWE(eval(modv$call, get("gmEnv", envir)))
	if (inherits(result$value, "condition")) return(result)

	fit1 <- result$value
	if(get("nextra", envir) != 0L) {
		extraResult1 <- get("applyExtras", envir)(fit1)
		nextra <- get("nextra", envir)
		if(length(extraResult1) < nextra) {
			tmp <- rep(NA_real_, nextra)
			tmp[match(names(extraResult1), get("extraResultNames", envir))] <-
				extraResult1
			extraResult1 <- tmp
		}
	} else extraResult1 <- NULL
	ll <- .getLik(fit1)$logLik(fit1)

	#mcoef <- matchCoef(fit1, all.terms = get("allTerms", envir),
	# beta = get("beta", envir), allCoef = TRUE)
	mcoef <- eval(get("matchCoefCall", envir))

	list(value = c(mcoef, extraResult1, df = attr(ll, "df"), ll = ll,
		ic = get("IC", envir)(fit1)),
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
