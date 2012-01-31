`dredge` <-
function(global.model, beta = FALSE, evaluate = TRUE, rank = "AICc",
		 fixed = NULL, m.max = NA, m.min = 0, subset, marg.ex = NULL,
		 trace = FALSE, varying, extra, ...) {

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

	#XXX: use.ranef <- FALSE
	#if(use.ranef && inherits(global.model, "mer")) {
		#allTerms <- c(allTerms, paste("(", attr(allTerms0, "random.terms"), ")",
			#sep = ""))
	#}

	if(length(grep(":", all.vars(reformulate(allTerms))) > 0L))
		stop("variable names in the formula cannot contain \":\"")

	logLik <- .getLogLik()

	# Check for na.omit
	if (!is.null(gmCall$na.action) &&
		as.character(gmCall$na.action) %in% c("na.omit", "na.exclude")) {
		stop("'global.model' should not use 'na.action' = ", gmCall$na.action)
	}

	if(names(gmCall)[2L] == "") names(gmCall)[2L] <-
		names(formals(deparse(gmCall[[1L]]))[1L])

	# TODO: other classes: model, fixed, etc...
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
		variants <- varying.names <- NULL
		nvariants <- 1L
		nvarying <- 0L
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

# XXX XXX: adjusting 'marg.ex'
		# newArgs <- makeArgs(global.model, allTerms, rep(TRUE, length(allTerms), argsOptions)
		# formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
			# attr(newArgs, "formulaList")
		# marg.ex <- unlist(lapply(lapply(formulaList, formulaAllowed), "attr", "marg.ex"))

		
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

			ord[k] <- iComb
			ret[k, retColIdx] <- row1
			coefTables[[k]] <- attr(mcoef1, "coefTable")
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

