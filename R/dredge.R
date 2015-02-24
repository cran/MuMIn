`dredge` <-
function(global.model, beta = FALSE, evaluate = TRUE, rank = "AICc",
		 fixed = NULL, m.max = NA, m.min = 0, subset,
		 trace = FALSE, varying, extra, ct.args = NULL,
		 ...) {
	
	trace <- min(as.integer(trace), 2L)

	gmEnv <- parent.frame()
	gmCall <- get_call(global.model)
	gmNobs <- nobs(global.model)

	if (is.null(gmCall)) {
		gmCall <- substitute(global.model)
		if(!is.call(gmCall)) {
			stop("need a 'global.model' with a call component. Consider using ", 
				if(inherits(global.model, c("gamm", "gamm4"))) 
					"'uGamm'" else "'updateable'")
		}
		#"For objects without a 'call' component the call to the fitting function \n",
		#" must be used directly as an argument to 'dredge'.")
		# NB: this is unlikely to happen
		if(!is.function(eval.parent(gmCall[[1L]])))
			cry(NA, "could not find function '%s'", asChar(gmCall[[1L]]))
	} else {
		# if 'update' method does not expand dots, we have a problem with
		# expressions like ..1, ..2 in the call. So try to replace them with
		# respective arguments in the original call
		isDotted <- grep("^\\.\\.", sapply(as.list(gmCall), deparse))
		if(length(isDotted) != 0L) {
			if(is.name(substitute(global.model))) {
				cry(NA, "call stored in 'global.model' contains unexpanded dots and cannot be updated: \n    %s",
					 asChar(gmCall))
			} else gmCall[isDotted] <-
				substitute(global.model)[names(gmCall[isDotted])]
		}
		
		# object from 'run.mark.model' has $call of 'make.mark.model' - fixing
		# it here:
		if(inherits(global.model, "mark") && gmCall[[1L]] == "make.mark.model") {
			gmCall <- call("run.mark.model", model = gmCall, invisible = TRUE)
		}
		
	}

	lik <- .getLik(global.model)
	logLik <- lik$logLik
	logLikName <- lik$name
	
	# *** Rank ***
	rank.custom <- !missing(rank)
	
	if(!rank.custom && logLikName == "qLik") {
		rank <- "QIC"
		cry(NA, "using 'QIC' instead of 'AICc'", warn = TRUE)
	}
	
	rankArgs <- list(...)

	if(any(badargs <- names(rankArgs) == "marg.ex")) {
		cry(NA, "argument \"marg.ex\" is defunct and has been ignored",
			 warn = TRUE)
		rankArgs <- rankArgs[!badargs]
	}
	if(any(names(rankArgs) == "na.action"))
		cry("RTFM", "argument \"na.action\" is inappropriate here",
			 warn = FALSE)

	IC <- .getRank(rank, rankArgs)
	
	if(any(badargs <- is.na(match(names(rankArgs),
		c(names(formals(get("rank", environment(IC))))[-1L], names(formals()))))))
		cry("RTFM", ngettext(sum(badargs),
			"argument %s is not a name of formal argument of %s",
			"arguments %s are not names of formal arguments of %s"),
			prettyEnumStr(names(rankArgs[badargs])), "'dredge' or 'rank'", 
			warn = TRUE)
	
	ICName <- as.character(attr(IC, "call")[[1L]])
	
	if(length(tryCatch(IC(global.model), error = function(e) {
		stop(simpleError(conditionMessage(e), subst(attr(IC, "call"),
			x = as.name("global.model"))))
	})) != 1L) {
		cry(NA, "result of '%s' is not of length 1", asChar(attr(IC, "call")))
	}

	allTerms <- allTerms0 <- getAllTerms(global.model, intercept = TRUE,
		data = eval(gmCall$data, envir = gmEnv))

	# Intercept(s)
	interceptLabel <- attr(allTerms, "interceptLabel")
	if(is.null(interceptLabel)) interceptLabel <- "(Intercept)"
	nIntercepts <- sum(attr(allTerms, "intercept"))

	#XXX: use.ranef <- FALSE
	#if(use.ranef && inherits(global.model, "mer")) {
		#allTerms <- c(allTerms, paste0("(", attr(allTerms0, "random.terms"), ")"))
	#}


	# Check for na.omit
	if(!(gmNaAction <- .checkNaAction(cl = gmCall, what = "'global.model'")))
		cry(NA, attr(gmNaAction, "message"))
	
	if(names(gmCall)[2L] == "") gmCall <-
		match.call(gmCall, definition = eval.parent(gmCall[[1L]]),
				   expand.dots = TRUE)

		
	# TODO: other classes: model, fixed, etc...
	gmCoefNames <- fixCoefNames(names(coeffs(global.model)))

	nVars <- length(allTerms)

	if(isTRUE(rankArgs$REML) || (isTRUE(.isREMLFit(global.model)) && is.null(rankArgs$REML)))
		cry(NA, "comparing models fitted by REML", warn = TRUE)

	if (beta && is.null(tryCatch(beta.weights(global.model), error = function(e) NULL,
		warning = function(e) NULL))) {
		cry(NA, "do not know how to calculate beta weights for '%s', argument 'beta' ignored",
			 class(global.model)[1L], warn = TRUE)
		beta <- FALSE
	}

	m.max <- if (missing(m.max)) (nVars - nIntercepts) else min(nVars - nIntercepts, m.max)

	# fixed variables:
	if (!is.null(fixed)) {
		if (inherits(fixed, "formula")) {
			if (fixed[[1L]] != "~" || length(fixed) != 2L)
				cry(NA, "'fixed' should be a one-sided formula", warn = TRUE)
			fixed <- as.vector(getAllTerms(fixed))
		} else if (identical(fixed, TRUE)) {
			fixed <- as.vector(allTerms[!(allTerms %in% interceptLabel)])
		} else if (!is.character(fixed)) {
			cry(NA, paste("'fixed' should be either a character vector with",
						   " names of variables or a one-sided formula"))
		}
		if (!all(i <- (fixed %in% allTerms))) {
			cry(NA, "some terms in 'fixed' do not exist in 'global.model': %s",
				 prettyEnumStr(fixed[!i]), warn = TRUE)
			fixed <- fixed[i]
		}
	}
	
	deps <- attr(allTerms0, "deps")
	fixed <- union(fixed, rownames(deps)[rowSums(deps, na.rm = TRUE) == ncol(deps)])
	fixed <- c(fixed, allTerms[allTerms %in% interceptLabel])
	
	nFixed <- length(fixed)
	if(nFixed != 0L) message(sprintf(ngettext(nFixed, "Fixed term is %s", "Fixed terms are %s"),
		prettyEnumStr(fixed)))

	termsOrder <- order(allTerms %in% fixed)
	allTerms <- allTerms[termsOrder]
	
	di <- match(allTerms, rownames(deps))
	deps <- deps[di, di]	
	
	gmFormulaEnv <- environment(as.formula(formula(global.model), env = gmEnv))
	# TODO: gmEnv <- gmFormulaEnv ???

	### BEGIN Manage 'varying'
	## @param:	varying
	## @value:	varying, varyingNames, variants, nVariants, nVarying
	if(!missing(varying) && !is.null(varying)) {
		nVarying <- length(varying)
		varyingNames <- names(varying)
		fvarying <- unlist(varying, recursive = FALSE, use.names = FALSE)
		vlen <- vapply(varying, length, 1L)
		nVariants <- prod(vlen)
		variants <- as.matrix(expand.grid(split(seq_len(sum(vlen)),
			rep(seq_len(nVarying), vlen))))
		
		variantsFlat <- unlist(lapply(varying, .makeListNames),
			recursive = FALSE, use.names = FALSE)
		
	} else {
		variants <- varyingNames <- NULL
		nVariants <- 1L
		nVarying <- 0L
	}
	## END: varying
	
	## BEGIN Manage 'extra'
	## @param:	extra, global.model, gmFormulaEnv, 
	## @value:	extra, nextra, extraNames, nullfit_
	if(!missing(extra) && length(extra) != 0L) {
		# a cumbersome way of evaluating a non-exported function in a parent frame:
		extra <- eval(as.call(list(call("get", ".get.extras",
			envir = call("asNamespace", .packageName), inherits = FALSE),
				substitute(extra), r2nullfit = TRUE)), parent.frame())
		
		#extra <- eval(call(".get.extras", substitute(extra), r2nullfit = TRUE), parent.frame())
		if(any(c("adjR^2", "R^2") %in% names(extra))) {
			nullfit_ <- null.fit(global.model, evaluate = TRUE, envir = gmFormulaEnv)
		}
		applyExtras <- function(x) unlist(lapply(extra, function(f) f(x)))
		extraResult <- applyExtras(global.model)
		if(!is.numeric(extraResult))
			cry(NA, "function in 'extra' returned non-numeric result")

		nextra <- length(extraResult)
		extraNames <- names(extraResult)
	} else {
		nextra <- 0L
		extraNames <- character(0L)
	}
	## END: manage 'extra'

	nov <- as.integer(nVars - nFixed)
	ncomb <- (2L ^ nov) * nVariants

	if(nov > 31L) cry(NA, "number of predictors (%d) exceeds allowed maximum of 31", nov)
	#if(nov > 10L) warning(gettextf("%d predictors will generate up to %.0f combinations", nov, ncomb))
	nmax <- ncomb * nVariants
	rvChunk <- 25L
	if(evaluate) {
		rvNcol <- nVars + nVarying + 3L + nextra
		ret <- matrix(NA_real_, ncol = rvNcol, nrow = rvChunk)
		coefTables <- vector(rvChunk, mode = "list")
	}

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
				
				subset0 <- subset
				subset <- matrix(subset[
					match(allTerms, rownames(subset)),
					match(allTerms, colnames(subset))],
					dimnames = list(allTerms, allTerms),
					nrow = n, ncol = n)
				tsubset <- t(subset)
				nas <- is.na(subset)
				i <- lower.tri(subset) & is.na(subset) & !t(nas)
				ti <- t(i)
				subset[i] <- subset[ti]
				subset[ti] <- NA
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

			#gloFactorTable <- t(attr(terms(global.model), "factors")[-1L, ] != 0)
			gloFactorTable <- t(attr(terms(reformulate(allTerms0[!(allTerms0
				%in% interceptLabel)])), "factors") != 0)
			rownames(gloFactorTable) <- allTerms0[!(allTerms0 %in% interceptLabel)]
			
			subsetExpr <- subset[[1L]]
			subsetExpr <- .exprapply(subsetExpr, ".", .sub_dot, gloFactorTable, 
				allTerms, as.name("comb"))
			subsetExpr <- .exprapply(subsetExpr, c("{", "Term"), .sub_Term)

			#@@@ TODO has subsetExpr <- .exprapply(subsetExpr, "has", .sub_Term)
			
			tmp <- updateDeps(subsetExpr, deps)
			subsetExpr <- tmp$expr
			deps <- tmp$deps

			subsetExpr <- .exprapply(subsetExpr, "dc", .sub_args_as_vars)
			subsetExpr <- .subst4Vec(subsetExpr, allTerms, as.name("comb"))
			
			if(nVarying) {
				ssValidNames <- c("cVar", "comb", "*nvar*")
				subsetExpr <- .exprapply(subsetExpr, "V", .sub_V, 
					as.name("cVar"), varyingNames)
				if(!all(all.vars(subsetExpr) %in% ssValidNames))
					subsetExpr <- .subst4Vec(subsetExpr, varyingNames,
											 as.name("cVar"), fun = "[[")
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
	
	comb.sfx <- rep(TRUE, nFixed)
	comb.seq <- if(nov != 0L) seq_len(nov) else 0L
	k <- 0L
	extraResult1 <- integer(0L)
	calls <- vector(mode = "list", length = rvChunk)
	ord <- integer(rvChunk)

	argsOptions <- list(
		response = attr(allTerms0, "response"),
		intercept = nIntercepts,
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
	
	retColIdx <- if(nVarying) -nVars - seq_len(nVarying) else TRUE
	
	if(trace > 1L) {
		progressBar <- if(.Platform$GUI == "Rgui") {
			 winProgressBar(max = ncomb, title = "'dredge' in progress")
		#} else if(capabilities("tcltk") && ("package:tcltk" %in% search())) {
			 #tkProgressBar(max = ncomb, title = "'dredge' in progress")
		} else txtProgressBar(max = ncomb, style = 3)
		setProgressBar <- switch(class(progressBar),
			    txtProgressBar = setTxtProgressBar,
			   #tkProgressBar = setTkProgressBar,
			   winProgressBar = setWinProgressBar,
			   function(...) {})
		on.exit(close(progressBar))
	}
	iComb <- -1L
	while((iComb <- iComb + 1L) < ncomb) {
		varComb <- iComb %% nVariants
		jComb <- (iComb - varComb) / nVariants

		if(varComb == 0L) {
			isok <- TRUE
			
			## comb : logical term indexes
			comb <- c(as.logical(intToBits(jComb)[comb.seq]), comb.sfx)
			nvar <- sum(comb) - nIntercepts
							
			if(nvar > m.max || nvar < m.min ||
			   !formula_margin_check(comb, deps) ||
			   switch(hasSubset,
					FALSE,
					!all(subset[comb, comb], na.rm = TRUE),
					!evalExprInEnv(subsetExpr, env = ssEnv, enclos = parent.frame(),
						comb = comb, `*nvar*` = nvar),
					FALSE
			   )) {
				isok <- FALSE
				next
			}
			
			newArgs <- makeArgs(global.model, allTerms[comb], comb, argsOptions)
			formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
				attr(newArgs, "formulaList")

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
		if (nVarying) {
			#cvi <- variants[(iComb - 1L) %% nvariants + 1L, ]
			cvi <- variants[varComb + 1L, ]
			if(hasSubset == 4L &&
				!evalExprInEnv(subsetExpr, env = ssEnv, enclos = parent.frame(),
					comb = comb, `*nvar*` = nvar, cVar = variantsFlat[cvi]))
						next;
			clVariant[varyingNames] <- fvarying[cvi]
		}

		if(trace == 1L) {
			cat(iComb, ": "); print(clVariant)
			utils::flush.console()
		} else if(trace == 2L) {
			setProgressBar(progressBar, value = iComb,
				title = sprintf("dredge: %d of %.0f subsets", k, (k / iComb) * ncomb))
		}
	
		if(evaluate) {
			# begin row1: (clVariant, gmEnv, modelId, IC(), applyExtras(),
			#              nextra, allTerms, beta,
			#              if(nVarying) variantsIdx[v] else NULL
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
			}

			mcoef1 <- eval(matchCoefCall)
			ll1 <- logLik(fit1)
			nobs1 <- nobs(fit1)
			if(nobs1 != gmNobs) cry(NA, "number of observations in model #%d (%d) different from global model (%d)",
				iComb, nobs1, gmNobs, warn = TRUE)

			row1 <- c(mcoef1[allTerms], extraResult1,
				df = attr(ll1, "df"), ll = ll1, ic = IC(fit1)
			)
			## end -> row1

			k <- k + 1L # all OK, add model to table
			rvlen <- nrow(ret)
			if(retNeedsExtending <- k > rvlen) { # append if necesarry
				nadd <- min(rvChunk, nmax - rvlen)
				ret <- rbind(ret, matrix(NA_real_, ncol = rvNcol, nrow = nadd),
					deparse.level = 0L)
				addi <- seq.int(rvlen + 1L, length.out = nadd)
				coefTables[addi] <- vector("list", nadd)
			}
			ret[k, retColIdx] <- row1
			coefTables[[k]] <- attr(mcoef1, "coefTable")
		} else { # if !evaluate
			k <- k + 1L
			rvlen <- length(ord)	
			if(retNeedsExtending <- k > rvlen) {
				nadd <- min(rvChunk, nmax - rvlen)
				addi <- seq.int(rvlen + 1L, length.out = nadd)
			}
		}
		if(retNeedsExtending) {
			calls[addi] <- vector("list", nadd)
			ord[addi] <- integer(nadd)
		}
		ord[k] <- iComb
		calls[[k]] <- clVariant
	} ### for (iComb ...)
	
	if(k == 0L) stop("result is empty")
	ord <- ord + 1L
	names(calls) <- ord
	if(!evaluate) return(calls[seq_len(k)])

	if(k < nrow(ret)) {
		i <- seq_len(k)
		ret <- ret[i, , drop = FALSE]
		ord <- ord[i]
		calls <- calls[i]
		coefTables <- coefTables[i]
	}

	if(nVarying) {
		varlev <- ord %% nVariants
		varlev[varlev == 0L] <- nVariants
		ret[, nVars + seq_len(nVarying)] <- variants[varlev, ]
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
	colnames(ret) <- c(allTerms, varyingNames, extraNames, "df", logLikName, ICName)

	if(nVarying) {
		variant.names <- vapply(variantsFlat, asChar, "", width.cutoff = 20L)

		vnum <- split(seq_len(sum(vlen)), rep(seq_len(nVarying), vlen))
		names(vnum) <- varyingNames
		for (i in varyingNames) ret[, i] <-
			factor(ret[, i], levels = vnum[[i]], labels = variant.names[vnum[[i]]])	
	}

	o <- order(ret[, ICName], decreasing = FALSE)
	ret <- ret[o, ]
	coefTables <- coefTables[o]

	ret$delta <- ret[, ICName] - min(ret[, ICName])
	ret$weight <- exp(-ret$delta / 2) / sum(exp(-ret$delta / 2))

	structure(ret,
		class = c("model.selection", "data.frame"),
		model.calls = calls[o],
		global = global.model,
		global.call = gmCall,
		terms = structure(allTerms, interceptLabel = interceptLabel),
		rank = IC,
		rank.call = attr(IC, "call"),
		beta = beta,
		call = match.call(expand.dots = TRUE),
		coefTables = coefTables,
		nobs = gmNobs,
		vCols = varyingNames
	)
} ######

`dredgeAll` <-
function(global.model, beta = FALSE, ...) {
	cl <- match.call(definition = dredge)
	cl$evaluate <- FALSE
	cl[[1L]] <- as.name("dredge")
	models <- lapply(eval.parent(cl), eval, parent.frame())
	ret <- model.sel(models)
	attr(ret, "modelList") <- models
	attr(ret, "global") <- global.model
	attr(ret, "global.call") <- get_call(global.model)
	attr(ret, "call") <- cl
	ret
}