`dredge` <-
function(global.model, beta = FALSE, evaluate = TRUE, rank = "AICc",
		 fixed = NULL, m.max = NA, m.min = 0, subset, marg.ex = NULL,
		 trace = FALSE, varying, extra, ct.args = NULL, ...) {

	gmEnv <- parent.frame()
	gmCall <- .getCall(global.model)
	gmNobs <- nobs(global.model)

	if (is.null(gmCall)) {
		gmCall <- substitute(global.model)
		if(!is.call(gmCall)) {
			stop("need a 'global.model' with call component. Consider using ", 
				if(inherits(global.model, c("gamm", "gamm4"))) 
					"'uGamm'" else "'updateable'")
		}
		#"For objects without a 'call' component the call to the fitting function \n",
		#" must be used directly as an argument to 'dredge'.")
		# NB: this is unlikely to happen:
		if(!is.function(eval(gmCall[[1L]], parent.frame())))
			gettext('could not find function "%s"', deparse(gmCall[[1L]], 
				control = NULL), domain = "R")
	} else {
		# if 'update' method does not expand dots, we have a problem
		# with expressions like ..1, ..2 in the call.
		# So, try to replace them with respective arguments in the original call
		is.dotted <- grep("^\\.\\.", sapply(as.list(gmCall), deparse))
		if(length(is.dotted) > 0L) {
			substGmCall <- substitute(global.model)
			if(is.name(substGmCall)) {
				.cry(NA, "call to 'global.model' contains '...' arguments and cannot be updated: %s", deparse(gmCall, control = NULL))
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
		.cry(NA, "using 'QIC' instead of 'AICc'", warn = TRUE)
	}
	
	rankArgs <- list(...)
	IC <- .getRank(rank, rankArgs)
	ICName <- as.character(attr(IC, "call")[[1L]])
	
	tryCatch(IC(global.model), error = function(e) {
		e$call <- do.call(substitute, list(attr(IC, "call"), list(x = as.name("global.model"))))
		stop(e)
	})

	allTerms <- allTerms0 <- getAllTerms(global.model, intercept = TRUE,
		data = eval(gmCall$data, envir = gmEnv))

	# Intercept(s)
	interceptLabel <- attr(allTerms, "interceptLabel")
	if(is.null(interceptLabel)) interceptLabel <- "(Intercept)"
	nInts <- sum(attr(allTerms, "intercept"))

	#XXX: use.ranef <- FALSE
	#if(use.ranef && inherits(global.model, "mer")) {
		#allTerms <- c(allTerms, paste("(", attr(allTerms0, "random.terms"), ")",
			#sep = ""))
	#}


	# Check for na.omit
	if (!is.null(gmCall$na.action) &&
		as.character(gmCall$na.action) %in% c("na.omit", "na.exclude")) {
		.cry(NA, "'global.model' should not use 'na.action' = \"%s\"", gmCall$na.action)
	}
	
	if(names(gmCall)[2L] == "") gmCall <-
		match.call(gmCall, definition = eval(gmCall[[1]], envir = parent.frame()), expand.dots = TRUE)
		
		
	# TODO: other classes: model, fixed, etc...
	gmCoefNames <- fixCoefNames(names(coeffs(global.model)))

	n.vars <- length(allTerms)

	if(isTRUE(rankArgs$REML) || (isTRUE(.isREMLFit(global.model)) && is.null(rankArgs$REML)))
		.cry(NA, "comparing models fitted by REML", warn = TRUE)

	if (beta && is.null(tryCatch(beta.weights(global.model), error = function(e) NULL,
		warning = function(e) NULL))) {
		.cry(NA, "do not know how to calculate beta weights for '%s', argument 'beta' ignored",
			 class(global.model)[1L], warn = TRUE)
		beta <- FALSE
	}

	m.max <- if (missing(m.max)) (n.vars - nInts) else min(n.vars - nInts, m.max)

	# fixed variables:
	if (!is.null(fixed)) {
		if (inherits(fixed, "formula")) {
			if (fixed[[1L]] != "~" || length(fixed) != 2L)
				.cry(NA, "'fixed' should be a one-sided formula", warn = TRUE)
			fixed <- as.vector(getAllTerms(fixed))
		} else if (identical(fixed, TRUE)) {
			fixed <- as.vector(allTerms[!(allTerms %in% interceptLabel)])
		} else if (!is.character(fixed)) {
			.cry(NA, "'fixed' should be either a character vector with"
				  + " names of variables or a one-sided formula")
		}
		if (!all(fixed %in% allTerms)) {
			.cry(NA, "not all terms in 'fixed' exist in 'global.model'", warn = TRUE)
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
			.cry(NA, "function in 'extra' returned non-numeric result")

		nextra <- length(extraResult)
		extraNames <- names(extraResult)
	} else {
		nextra <- 0L
		extraNames <- character(0L)
	}
	## END: manage 'extra'

	nov <- as.integer(n.vars - n.fixed)
	ncomb <- (2L ^ nov) * nvariants

	if(nov > 31L) .cry(NA, "number of predictors (%d) exceeds allowed maximum (31)", nov)
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
				if(any(di != n)) stop("unnamed 'subset' matrix does not have both dimensions",
					" equal to number of terms in 'global.model': %d", n)
				
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

			subsetExpr <- subset[[1L]]

			## subset X
			#gloFactorTable <- t(attr(terms(global.model), "factors")[-1L, ] != 0)
			gloFactorTable <- t(attr(terms(reformulate(allTerms0[!(allTerms0
				%in% interceptLabel)])), "factors") != 0)
			
			rownames(gloFactorTable) <- allTerms0[!(allTerms0 %in% interceptLabel)]
	
			
			subsetExpr <- .substFun4Fun(subsetExpr, ".", function(x, fac, at, vName) {
				if(length(x) != 2L) .cry(x, "exactly one argument needed, %d given.", length(x) - 1L)
				if(length(x[[2L]]) == 2L && x[[2L]][[1L]] == "+") {
					fun <- "all"
					sx <- as.character(x[[2L]][[2L]])
				} else {
					fun <- "any"
					sx <- as.character(x[[2L]])
				}
				#print(sx)
				dn <- dimnames(fac)
				#print(dn)
				#browser()
				if(!(sx %in% dn[[2L]])) .cry(x, "unknown variable name '%s'", sx)
				as.call(c(as.name(fun), call("[", vName, as.call(c(as.name("c"), 
					match(dn[[1L]][fac[, sx]], at))))))
			}, gloFactorTable, allTerms, as.name("comb"))
			
			subsetExpr <- .subst4Vec(subsetExpr, allTerms, as.name("comb"))
			
			
			if(nvarying) {
				ssValidNames <- c("cVar", "comb", "*nvar*")
				subsetExpr <- .substFun4Fun(subsetExpr, "V", function(x, cVar, fn) {
					if(length(x) > 2L) .cry(x, "discarding extra arguments", warn = TRUE)
					i <- which(fn == x[[2L]])[1L]
					if(is.na(i)) .cry(x, "'%s' is not a valid name of 'varying' element",
									  as.character(x[[2L]]), warn = TRUE)
					call("[[", cVar, i)
				}, as.name("cVar"), varying.names)
				if(!all(all.vars(subsetExpr) %in% ssValidNames))
					subsetExpr <- .subst4Vec(subsetExpr, varying.names,
											 as.name("cVar"), fun = "[[")
			}
			
			
			# DebugPrint(subsetExpr)
			
			ssVars <- all.vars(subsetExpr)
			okVars <- ssVars %in% ssValidNames
			if(!all(okVars)) stop("unrecognized names in 'subset' expression: ",
				prettyEnumStr(ssVars[!okVars]))

			ssEnv <- new.env(parent = .GlobalEnv)
			ssFunc <- setdiff(all.vars(subsetExpr, functions = TRUE), ssVars)
			if("dc" %in% ssFunc) assign("dc", .subset_dc, ssEnv)
			
			hasSubset <- if(any(ssVars == "cVar")) 4L else # subset as expression
				3L # subset as expression using 'varying' variables
				
			#DebugPrint(subsetExpr)
			

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
	
	matchCoefCall <- as.call(c(alist(matchCoef, fit1, all.terms = allTerms,
		  beta = beta, allCoef = TRUE), ct.args))
	
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

	retColIdx <- if(nvarying) -n.vars - seq_len(nvarying) else TRUE

	prevJComb <- 0L
	#DebugPrint(ncomb)
	#DebugPrint(nvariants)
	for(iComb in seq.int(ncomb)) {
		jComb <- ceiling(iComb / nvariants)
		#DebugPrint(iComb)
		#DebugPrint(jComb)
		if(jComb != prevJComb) {
			isok <- TRUE
			prevJComb <- jComb
			
			#DebugPrint(hasSubset)
			comb <- c(as.logical(intToBits(jComb - 1L)[comb.seq]), comb.sfx)
			nvar <- sum(comb) - nInts
			#DebugPrint(nInts)
				
			if(nvar > m.max || nvar < m.min ||
			   switch(hasSubset,
					FALSE,
					!all(subset[comb, comb], na.rm = TRUE),
					!.evalExprIn(subsetExpr, env = ssEnv, enclos = parent.frame(),
						comb = comb, `*nvar*` = nvar),
					FALSE
			   )) {
				isok <- FALSE
				next;
			}
			newArgs <- makeArgs(global.model, allTerms[comb], comb, argsOptions)
			formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
				attr(newArgs, "formulaList")

			
			if(!all(vapply(formulaList, formulaMargChk, logical(1L), marg.ex)))  {
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
		if (nvarying) {
			cvi <- variants[(iComb - 1L) %% nvariants + 1L, ]

			if(hasSubset == 4L &&
				!.evalExprIn(subsetExpr, env = ssEnv, enclos = parent.frame(),
					comb = comb, `*nvar*` = nvar, cVar = flat.variant.Vvals[cvi]))
						next;
			clVariant[varying.names] <- fvarying[cvi]
		}

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

			#mcoef1 <- matchCoef(fit1, all.terms = allTerms, beta = beta,
			#	allCoef = TRUE)
			mcoef1 <- eval(matchCoefCall)

			ll <- logLik(fit1)
			nobs1 <- nobs(fit1)
			if(nobs1 != gmNobs) warning(gettextf(
				"number of observations in model #%d (%d) differs from that in the global model (%d)",
				iComb, nobs1, gmNobs))

			row1 <- c(mcoef1[allTerms], extraResult1,
				df = attr(ll, "df"), ll = ll, ic = IC(fit1)
			)
			## end -> row1

			k <- k + 1L # all OK, add model to table

			ret.nrow <- nrow(ret)
			if(k > ret.nrow) { # append if necesarry
				nadd <- min(ret.nchunk, nmax - ret.nrow)
				coefTables <- c(coefTables, vector(nadd, mode = "list"))
				ret <- rbind(ret, matrix(NA, ncol = ret.ncol, nrow = nadd),
					deparse.level = 0L)
					calls <- c(calls, vector("list", nadd))
				ord <- c(ord, integer(nadd))
			}
			ret[k, retColIdx] <- row1
			coefTables[[k]] <- attr(mcoef1, "coefTable")
		} else { # if !evaluate
			k <- k + 1L # all OK, add model to table
		}
		ord[k] <- iComb
		calls[[k]] <- clVariant
	} ### for (iComb ...)

	if(k == 0L) stop("result is empty")
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
		coefTables = coefTables,
		nobs = gmNobs,
		vCols = varying.names
	)
} ######
