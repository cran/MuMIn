
# code snippet to handle argument 'beta'
.expr_beta_arg <- 
expression({
	if(is.logical(beta) && beta) {
		betaMode <- as.integer(beta)
		strbeta <- if(beta) "sd" else "none"
	} else if(is.character(beta)) {
		strbeta <- match.arg(beta)
		beta <- strbeta != "none"
		betaMode <- (strbeta != "none") + (strbeta == "partial.sd")
	} else {
        cry(, "invalid value for 'beta' : the argument is taken to be \"none\"",
            warn = TRUE)
		betaMode <- 0L
		strbeta <- "none"
	}
})


dredge <-
function(global.model, beta = c("none", "sd", "partial.sd"), 
         evaluate = TRUE, rank = "AICc",
		 fixed = NULL, m.lim = NULL, m.min, m.max, subset,
		 trace = FALSE, varying, extra, ct.args = NULL,
		 deps = attr(allTerms0, "deps"),
         cluster = NULL,
		 ...) {
         
         
    if(isTRUE(evaluate) && inherits(cluster, "cluster")) {
        cl <- match.call()
        cl[[1L]] <- quote(.dredge.par)
        return(eval(cl))
    }
    

	trace <- min(as.integer(trace), 2L)
	strbeta <- betaMode <- NULL
	eval(.expr_beta_arg)

	gmEnv <- parent.frame()
	gmNobs <- nobs(global.model)

	gmCall <- get_call(global.model)
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
		if(!is.function(eval(gmCall[[1L]], gmEnv)))
			cry(, "could not find function '%s'", asChar(gmCall[[1L]]))
	} else {
		# if 'update' method does not expand dots, we have a problem with
		# expressions like ..1, ..2 in the call. So try to replace them with
		# respective arguments in the original call
		isDotted <- grep("^\\.\\.", sapply(as.list(gmCall), asChar))
		if(length(isDotted) != 0L) {
			if(is.name(substitute(global.model))) {
				cry(, "the call stored in 'global.model' contains dotted names and cannot be updated. \n    Consider using 'updateable' on the modelling function")
			} else gmCall[isDotted] <-
				substitute(global.model)[names(gmCall[isDotted])]
		}
		# object from 'run.mark.model' has $call of 'make.mark.model' - fixing
		# it here:
		if(inherits(global.model, "mark") && gmCall[[1L]] == "make.mark.model") {
			gmCall <- call("run.mark.model", model = gmCall, invisible = TRUE)
		}
	}
	
	thisCall <- sys.call()
	exprApply(gmCall[["data"]], NA, function(expr) {
		if(!is.symbol(expr[[1L]]) || all(expr[[1L]] != c("@", "$", "::")))
			cry(thisCall, "'global.model' uses \"data\" that is a function value: use a variable instead")
	})

	lik <- .getLik(global.model)
	logLik <- lik$logLik

	# *** Rank ***
	rank.custom <- !missing(rank)

	if(!rank.custom && lik$name == "qLik") {
		rank <- "QIC"
		cry(, "using 'QIC' instead of 'AICc'", warn = TRUE)
	}

	rankArgs <- list(...)

	if(any(badargs <- names(rankArgs) == "marg.ex")) {
		cry(, "argument \"marg.ex\" is defunct and has been ignored",
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
		cry(, "result of '%s' is not of length 1", asChar(attr(IC, "call")))
	}

	allTerms <- allTerms0 <- getAllTerms(global.model, intercept = TRUE,
		data = eval(gmCall$data, envir = gmEnv))

	# Intercept(s)
	interceptLabel <- attr(allTerms, "interceptLabel")
	if(is.null(interceptLabel)) interceptLabel <- "(Intercept)"
	nIntercepts <- sum(attr(allTerms, "intercept"))

	# Check for na.omit
	if(!(gmNaAction <- .checkNaAction(cl = gmCall, what = "'global.model'", envir = gmEnv)))
		cry(, attr(gmNaAction, "message"))

	if(!nzchar(names(gmCall)[2L]))
		gmCall <-
			match.call(gmCall, definition = eval(gmCall[[1L]], envir = gmEnv),
				expand.dots = TRUE)


	# TODO: other classes: model, fixed, etc...
    gmCoefNames <- names(coeffs(global.model))
    if(any(dup <- duplicated(gmCoefNames)))
        cry(, "model cannot have duplicate coefficient names: ",
             prettyEnumStr(gmCoefNames[dup]))

	gmCoefNames <- fixCoefNames(gmCoefNames)

	nVars <- length(allTerms)

	if(isTRUE(rankArgs$REML) || (isTRUE(.isREMLFit(global.model)) && is.null(rankArgs$REML)))
		cry(, "comparing models fitted by REML", warn = TRUE)

	if ((betaMode != 0L) && is.null(tryCatch(std.coef(global.model, betaMode == 2L),
		error = return_null, warning = return_null))) {
		cry(, "do not know how to standardize coefficients of '%s': argument 'beta' ignored",
			 class(global.model)[1L], warn = TRUE)
		betaMode <- 0L
		strbeta <- "none"
	}

	if(nomlim <- is.null(m.lim)) m.lim <- c(0, NA)
	## XXX: backward compatibility:
	if(!missing(m.max) || !missing(m.min)) {
		warning("arguments 'm.min' and 'm.max' are deprecated, use 'm.lim' instead")
		if(!nomlim) stop("cannot use both 'm.lim' and 'm.min' or 'm.max'")
		if(!missing(m.min)) m.lim[1L] <- m.min[1L]
		if(!missing(m.max)) m.lim[2L] <- m.max[1L]
	}
	if(!is.numeric(m.lim) || length(m.lim) != 2L || any(m.lim < 0, na.rm = TRUE))
		stop("invalid 'm.lim' value")
	m.lim[2L] <- if (!is.finite(m.lim[2L])) (nVars - nIntercepts) else
		min(nVars - nIntercepts, m.lim[2L])
	if (!is.finite(m.lim[1L])) m.lim[1L] <- 0
	m.min <- m.lim[1L]
    m.max <- m.lim[2L]

	# fixed variables:
	if (!is.null(fixed)) {
		if (inherits(fixed, "formula")) {
			if (fixed[[1L]] != "~" || length(fixed) != 2L)
				cry(, "'fixed' should be a one-sided formula", warn = TRUE)
			fixed <- as.vector(getAllTerms(fixed))
		} else if (identical(fixed, TRUE)) {
			fixed <- as.vector(allTerms[!(allTerms %in% interceptLabel)])
		} else if (!is.character(fixed)) {
			cry(, paste("'fixed' should be either a character vector with",
						   " names of variables or a one-sided formula"))
		}
		if (!all(i <- (fixed %in% allTerms))) {
			cry(, "some terms in 'fixed' do not exist in 'global.model': %s",
				 prettyEnumStr(fixed[!i]), warn = TRUE)
			fixed <- fixed[i]
		}
	}
	fixed <- union(fixed, rownames(deps)[rowSums(deps, na.rm = TRUE) == ncol(deps)])
	fixed <- c(fixed, allTerms[allTerms %in% interceptLabel])
    fixed <- fixed[!duplicated(fixed)]
	nFixed <- length(fixed)
	if(nFixed != 0L) message(sprintf(ngettext(nFixed, "Fixed term is %s", "Fixed terms are %s"),
		prettyEnumStr(fixed)))

	termsOrder <- order(allTerms %in% fixed)
	allTerms <- allTerms[termsOrder]

	di <- match(allTerms, rownames(deps))
	deps <- deps[di, di, drop = FALSE]

	gmFormulaEnv <- environment(as.formula(formula(global.model), env = gmEnv))
	# TODO: gmEnv <- gmFormulaEnv ???

	### BEGIN 'varying'
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

	## BEGIN 'extra'
	## @param:	extra, global.model, gmFormulaEnv,
	## @value:	extra, nExtra, extraNames, nullfit_
	if(!missing(extra) && length(extra) != 0L) {
    
        r2inExtra <- is.vector(extra) && any(c("adjR^2", "R^2") %in% extra)
		
		if (r2inExtra) {
            if(nVariants > 1L)
                stop("\"R^2\" in 'extra' cannot be used when 'varying' is given")
            nullfit <- null.fit(global.model, evaluate = TRUE, envir = gmFormulaEnv)
        } else {
            nullfit <- NULL
        }
          
        extra <- eval.parent(call(".get.extras", substitute(extra), r2nullfit = nullfit))

		extraResult <- .applyExtras(global.model, extra)
		nExtra <- length(extraResult)
		extraNames <- names(extraResult)
	} else {
        extra <- NULL
		nExtra <- 0L
		extraNames <- character()
	}
	## END: 'extra'

	nov <- as.integer(nVars - nFixed)
	ncomb <- (2L ^ nov) * nVariants
    novMax <- log2(.Machine$integer.max %/% nVariants)
    if(nov > novMax)
		cry(, "number of non-fixed predictors [%d] exceeds the allowed maximum of %.0f (with %d variants)", 
            nov, novMax, nVariants)
    resultChunkSize <- 25L
	if(evaluate) {
		rvNcol <- nVars + nVarying + 3L + nExtra
		rval <- matrix(NA_real_, ncol = rvNcol, nrow = resultChunkSize)
		coefTables <- vector(resultChunkSize, mode = "list")
	}

	## BEGIN: 'subset'
	## @param:	hasSubset, subset, allTerms, [interceptLabel],
	## @value:	hasSubset, subset
	if(missing(subset))  {
		hasSubset <- 1L
	} else {
		if(!tryCatch(is.language(subset) || is.matrix(subset) || 
           is.null(subset) || isTRUE(subset), error = function(e) FALSE))
			subset <- substitute(subset)
            
        if(is.null(subset) || isTRUE(subset)) {
            hasSubset <- 1L
        } else if(is.matrix(subset)) {
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
				nas <- is.na(subset)
				lotri <- lower.tri(subset)
				i <- lotri & nas & !t(nas)
				subset[i] <- t(subset)[i]
				subset[!lotri] <- NA

			}
			if(any(!is.na(subset[!lower.tri(subset)]))) {
				warning("non-missing values exist outside the lower triangle of 'subset'")
				subset[!lower.tri(subset)] <- NA
			}
			mode(subset) <- "logical"
			hasSubset <- 2L # subset as matrix

		} else { # subset is not null, TRUE or a matrix
			if(inherits(subset, "formula")) {
				if (subset[[1L]] != "~" || length(subset) != 2L)
					stop("'subset' formula should be one-sided")
				subset <- subset[[2L]]
			}
			subset <- as.expression(subset)
			ssValidNames <- c("comb", "*nvar*")


			tmpTerms <- terms(reformulate(allTerms0[!(allTerms0 %in% interceptLabel)]))
			gloFactorTable <- t(attr(tmpTerms, "factors") != 0)

			offsetNames <- sapply(attr(tmpTerms, "variables")[attr(tmpTerms, "offset") + 1L], asChar)

			if(length(offsetNames) != 0L) {
				gloFactorTable <- rbind(gloFactorTable,
					matrix(FALSE, ncol = ncol(gloFactorTable), nrow = length(offsetNames),
						dimnames = list(offsetNames, NULL)))
				for(i in offsetNames) gloFactorTable[offsetNames, offsetNames] <- TRUE
				#Note `diag<-` does not work for x[1x1] matrix:
				# diag(gloFactorTable[offsetNames, offsetNames, drop = FALSE]) <- TRUE
			}
			
			.DebugPrint(gloFactorTable)

			# fix interaction names in rownames:
			rownames(gloFactorTable) <- allTerms0[!(allTerms0 %in% interceptLabel)]

			subsetExpr <- subset[[1L]]
			subsetExpr <- exprapply0(subsetExpr, c("with", "."), .subst.with, gloFactorTable,
				allTerms, as.name("comb"), gmEnv)

			subsetExpr <- exprapply0(subsetExpr, c("{", "Term"), .subst.term)

			#@@@ TODO has subsetExpr <- exprapply0(subsetExpr, "has", .subst.term)

			tmp <- updateDeps(subsetExpr, deps)
			subsetExpr <- tmp$expr
			deps <- tmp$deps
			
			subsetExpr <- exprapply0(subsetExpr, "dc", .subst.vars.for.args)
			subsetExpr <- .subst.names.for.items(subsetExpr, allTerms, "comb")

			if(nVarying) {
				ssValidNames <- c("cVar", "comb", "*nvar*")
				subsetExpr <- exprapply0(subsetExpr, "V", .subst.v,
					as.name("cVar"), varyingNames)
				if(!all(all.vars(subsetExpr) %in% ssValidNames))
					subsetExpr <- .subst.names.for.items(subsetExpr, varyingNames,
											 "cVar", fun = "[[")
			}
			ssVars <- all.vars(subsetExpr)
			okVars <- ssVars %in% ssValidNames
			if(!all(okVars)) stop("unrecognized names in 'subset' expression: ",
				prettyEnumStr(ssVars[!okVars]))

			ssEnv <- new.env(parent = parent.frame())
			ssFunc <- setdiff(all.vars(subsetExpr, functions = TRUE), ssVars)
			if("dc" %in% ssFunc) assign("dc", .subset_dc, ssEnv)

			hasSubset <- if(any(ssVars == "cVar")) 4L else # subset as expression
				3L # subset as expression using 'varying' variables

		}
	} # END: 'subset'

	comb.sfx <- rep(TRUE, nFixed)
	comb.seq <- if(nov != 0L) seq_len(nov) else 0L
	k <- 0L
	extraResult1 <- integer(0L)
	calls <- vector(mode = "list", length = resultChunkSize)
	ord <- integer(resultChunkSize)

	argsOptions <- list(
		response = attr(allTerms0, "response"),
		intercept = nIntercepts,
		interceptLabel = interceptLabel,
		random = attr(allTerms0, "random"),
		gmCall = gmCall,
		gmEnv = gmEnv,
		allTerms = allTerms0,
		gmCoefNames = gmCoefNames,
		## TODO: is 'gmDataHead' needed anymore?
		gmDataHead = if(!is.null(gmCall$data)) {
			if(eval(call("is.data.frame", gmCall$data), gmEnv))
				eval(call("head", gmCall$data, 1L), gmEnv) else gmCall$data
			} else NULL,
		gmFormulaEnv = gmFormulaEnv
		)

## [[end of common code]]

	matchCoefCall <- as.call(c(alist(matchCoef, fit1, all.terms = allTerms,
		  beta = betaMode, allCoef = TRUE), ct.args))

	retColIdx <- if(nVarying) -nVars - seq_len(nVarying) else TRUE

	dotrace <- if(trace == 1L) {
		dotrace <- function()  {
			cat(iComb, ": "); print(clVariant)
			utils::flush.console()
		}
	} else if(trace > 1L) {
		progressBar <- .progbar(max = ncomb, title = "\"dredge\" working...")
		on.exit(.closeprogbar(progressBar))
		function() progressBar(value = iComb,
			title = sprintf("dredge: %d of ca. %.0f subsets (%d total)", k, (k / iComb) * ncomb, iComb))
	} else function() {}
	

	iComb <- -1L
	while((iComb <- iComb + 1L) < ncomb) {
		varComb <- iComb %% nVariants
		jComb <- (iComb - varComb) %/% nVariants

		#if(iComb %% 100L == 0L) setProgressBar(progressBar, value = iComb, title = sprintf("dredge: %d/%d total", k, iComb))

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
			
			newArgs <- makeArgs(global.model, allTerms[comb], argsOptions) # comb
				
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

		dotrace()

		if(evaluate) {
			# begin row1: (clVariant, gmEnv, modelId, IC(), applyExtras(),
			#              nExtra, allTerms, beta,
			#              if(nVarying) variantsIdx[v] else NULL
			fit1 <- tryCatch(eval(clVariant, gmEnv), error = function(err) {
				err$message <- paste(conditionMessage(err), "(model",
					iComb, "skipped)", collapse = "")
				class(err) <- c("simpleError", "warning", "condition")
				warning(err)
				return(NULL)
			})

			if (is.null(fit1)) next;

			if(nExtra != 0L) {
				extraResult1 <- .applyExtras(fit1, extra)
				if(length(extraResult1) < nExtra) {
					tmp <- rep(NA_real_, nExtra)
					tmp[match(names(extraResult1), names(extraResult))] <- extraResult1
					extraResult1 <- tmp
				}
			}

			mcoef1 <- eval(matchCoefCall)

			ll1 <- logLik(fit1)
			nobs1 <- nobs(fit1)
            
             
			if(nobs1 != gmNobs) cry(, "model #%d [%d] is fitted to a different number of observations than the global model [%d]",
				iComb, nobs1, gmNobs, warn = TRUE)

			row1 <- c(mcoef1[allTerms], extraResult1,
				df = attr(ll1, "df"), ll = ll1, ic = IC(fit1)
			)
			## end -> row1

			k <- k + 1L # all OK, add model to table
			rvlen <- nrow(rval)
			if(retNeedsExtending <- k > rvlen) { # append if necesarry
				nadd <- min(resultChunkSize, ncomb - rvlen)
				rval <- rbind(rval, matrix(NA_real_, ncol = rvNcol, nrow = nadd),
					deparse.level = 0L)
				addi <- seq.int(rvlen + 1L, length.out = nadd)
				coefTables[addi] <- vector("list", nadd)
			}
			rval[k, retColIdx] <- row1
			coefTables[[k]] <- attr(mcoef1, "coefTable")
		} else { # if !evaluate
			k <- k + 1L
			rvlen <- length(ord)
			if(retNeedsExtending <- k > rvlen) {
				nadd <- min(resultChunkSize, ncomb - rvlen)
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

	if(k < nrow(rval)) {
		i <- seq_len(k)
		rval <- rval[i, , drop = FALSE]
		ord <- ord[i]
		calls <- calls[i]
		coefTables <- coefTables[i]
	}

	if(nVarying) {
		varlev <- ord %% nVariants
		varlev[varlev == 0L] <- nVariants
		rval[, nVars + seq_len(nVarying)] <- variants[varlev, ]
	}

	rval <- as.data.frame(rval, stringsAsFactors = TRUE)
	row.names(rval) <- ord

	# Convert columns with presence/absence of terms to factors
	tfac <- which(!(allTerms %in% gmCoefNames))
	rval[tfac] <- lapply(rval[tfac], factor, levels = NaN, labels = "+")
	rval[, seq_along(allTerms)] <- rval[, v <- order(termsOrder)]
	allTerms <- allTerms[v]

    colnames(rval) <- c(allTerms, varyingNames, extraNames, "df", lik$name, ICName)
	if(nVarying) {
		variant.names <- vapply(variantsFlat, asChar, "", width.cutoff = 20L)

		vnum <- split(seq_len(sum(vlen)), rep(seq_len(nVarying), vlen))
		names(vnum) <- varyingNames
		for (i in varyingNames) rval[, i] <-
			factor(rval[, i], levels = vnum[[i]], labels = variant.names[vnum[[i]]])
	}

	rval <- rval[o <- order(rval[, ICName], decreasing = FALSE), ]
	coefTables <- coefTables[o]

	rval$delta <- rval[, ICName] - rval[1L, ICName]
	rval$weight <- Weights(rval$delta)
    mode(rval$df) <- "integer"
	
	structure(rval,
		model.calls = calls[o],
		global = global.model,
		global.call = gmCall,
		terms = structure(allTerms, interceptLabel = interceptLabel),
		rank = IC,
		beta = strbeta, #eval(formals(sys.function())[["beta"]])[betaMode + 1L],
		call = match.call(expand.dots = TRUE),
		coefTables = coefTables,
		nobs = gmNobs,
		vCols = varyingNames,
		column.types = {
			colTypes <- c(terms = length(allTerms), varying = length(varyingNames),
				extra = length(extraNames), df = 1L, loglik = 1L, ic = 1L, delta = 1L,
				weight = 1L)
			column.types <- rep(1L:length(colTypes), colTypes)
			names(column.types) <- colnames(rval)
			lv <- 1L:length(colTypes)
			factor(column.types, levels = lv, labels = names(colTypes)[lv])
		},
        extra = extra,
        class = c("model.selection", "data.frame")
	)
} #



`dredgeAll` <-
function(global.model, beta = FALSE, ...) {
	cl <- match.call(definition = dredge)
	cl$evaluate <- FALSE
	cl[[1L]] <- as.name("dredge")
	models <- lapply(eval.parent(cl), eval, parent.frame())
	rval <- model.sel(models)
	attr(rval, "modelList") <- models
	attr(rval, "global") <- global.model
	attr(rval, "global.call") <- get_call(global.model)
	attr(rval, "call") <- cl
	rval
}
