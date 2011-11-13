`dredge.old` <-
function(global.model, beta = FALSE, evaluate = TRUE, rank = "AICc",
		 fixed = NULL, m.max = NA, m.min = 0, subset, marg.ex = NULL, trace = FALSE,
		 varying = NULL, ...) {

	rank.custom <- !missing(rank)

	# *** Rank ***
	# switch to QAICc if quasi* family and no rank
	#if(inherits(global.model, "glm") && family(global.model)$family %in%
	#	c("quasi", "quasibinomial", "quasipoisson") && !rank.custom) {
	#	rankArgs <- list(chat=summary(global.model)$dispersion)
	#	rank <- "QAICc"
	#	warning("QAICc used for '", family(global.model)$family,
	#			"' family with c-hat = ", signif(rankArgs$chat))
	#	rank.custom <- TRUE
	#} else {
		rankArgs <- list(...)
	#}
	IC <- .getRank(rank, rankArgs)

	ICName <- as.character(attr(IC, "call")[[1L]])
	# *** Rank ***

	all.terms <- getAllTerms(global.model, intercept = TRUE)
	interceptLabel <- attr(all.terms, "interceptLabel")
	if(is.null(interceptLabel)) interceptLabel <- "(Intercept)"

	nInts <- length(interceptLabel)

		# Just in case:
	gterms <- tryCatch(terms(formula(global.model)),
		error=function(...) terms(global.model))

	response <- if(attr(gterms, "response") == 0L) NULL else "."


	if(length(grep(":", all.vars(delete.response(gterms) > 0L))))
		stop("Variable names in the model formula cannot contain \":\"")

	global.call <- .getCall(global.model)

	if (is.null(global.call)) {
		global.call <- substitute(global.model)
		if(!is.call(global.call)) {
			if(inherits(global.model, c("gamm", "gamm4")))
				message("To use gamm models with 'dredge', use 'MuMIn::gamm' wrapper")
			stop("Could not retrieve the call to 'global.model'")
		}
		#"For objects without a 'call' component the call to the fitting function \n",
		#" must be used directly as an argument to 'dredge'.")

		# this is unlikely to happen:
		if(!exists(as.character(global.call[[1L]]), parent.frame(), mode="function"))
			 stop("Could not find function '", global.call[[1L]], "'")


		formula.arg <- 2L # assume is a first argument
	} else {


		# if 'update' method does not expand dots, we have a problem
		# with expressions like ..1, ..2 in the call.
		# So, try to replace them with respective arguments in the original call
		is.dotted <- grep("^\\.\\.", sapply(as.list(global.call), deparse))
		if(length(is.dotted) > 0L)
			global.call[is.dotted] <-
				substitute(global.model)[names(global.call[is.dotted])]

		# if we have a call, try to figure out the 'formula' argument name
		formula.arg <- if(inherits(global.model, "lme"))	"fixed" else
			if(inherits(global.model, "gls")) "model" else {
				tryCatch({
					arg <- names(formals(match.fun(global.call[[1L]])))
					if ("formula" %in% arg) "formula" else 2L
				}, error = function(e) {
					2L
				})
			}
	}

	global.formula <- tryCatch(formula(global.model), error = function(...) global.call[[formula.arg]])
	# if (!missing(formula)) {
		# global.formula <- update.formula(global.formula, formula)
		# rm(formula)
	# }

	globCoefNames0 <- names(coeffs(global.model))
	if(has.start <- !is.null(global.call$start)) {
		if(is.null(names(global.call$start))) {
			names(global.call$start)[-1] <- names(coeffs(global.model))
		}
	}

	# Check for na.omit
	if (!is.null(global.call$na.action) &&
		as.character(global.call$na.action) %in% c("na.omit", "na.exclude")) {
		stop("'global.model' should not use 'na.action' = ", global.call$na.action)
	}

	has.int <- attr(all.terms, "intercept")
	globCoefNames <- fixCoefNames(globCoefNames0)
	#sglobCoefNames <- fixCoefNames(names(coeffs(global.model)))

	n.vars <- length(all.terms)
	ret <- numeric(0L)
	formulas <- character(0L)

	is.lm <- !inherits(global.model, "glm") & inherits(global.model, "lm")

	if(isTRUE(rankArgs$REML) || (isTRUE(.isREMLFit(global.model)) && is.null(rankArgs$REML)))
		warning("Comparing models with different fixed effects fitted by REML")

	if (beta && is.null(tryCatch(beta.weights(global.model), error=function(e) NULL,
		warning=function(e) NULL))) {
		warning("Do not know how to calculate B-weights for ",
				class(global.model)[1L], ", argument 'beta' ignored")
		beta <- FALSE
	}

	summary.globmod <- summary(global.model)

	has.rsq <- is.list(summary.globmod) && is.numeric(summary.globmod$r.squared)
	has.dev <- !is.null(deviance(global.model))

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
		if (!all(fixed %in% all.terms)) {
			warning("Not all terms in 'fixed' exist in 'global.model'")
			fixed <- fixed[fixed %in% all.terms]
		}
	}

	fixed <- c(fixed, all.terms[all.terms %in% interceptLabel])
	#fixed <- c(fixed, interceptLabel)

	n.fixed <- length(fixed)
	termsOrder <- order(all.terms %in% fixed)
	ordAllTerms <- all.terms[termsOrder]
	mostattributes(ordAllTerms) <- 	attributes(all.terms)
	all.terms <- ordAllTerms
	rm(ordAllTerms)

	llik <- .getLogLik()
	getK <- function(x) as.vector(attr(llik(x), "df"))

	isMER <- any(inherits(global.model, c("mer", "lmer", "glmer")))
	env <- attr(gterms, ".Environment")

	# DEBUG <- function(x) cat("*", deparse(substitute(x)), "=", x, "*\n")

	### BEGIN:

	## varying BEGIN
	if(!missing(varying)) {
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

	#DebugPrint(n.vars)
	#DebugPrint(n.fixed)

	nov <- as.integer(n.vars - n.fixed)
	ncomb <- 2L ^ nov

	if(nov > 31L) stop(gettextf("Maximum number of predictors is 31, but %d is given", nov))
	if(nov > 10L) warning(gettextf("%d predictors will generate up to %.0f possible combinations", nov, ncomb))

	#binPos <- 2L ^ seq.int(0L, nov - 1L)

	if(evaluate) {
		ret.nchunk <- 25L
		ret.ncol <- n.vars + (2L * has.rsq) + has.dev + 2L + nvarying
		ret <- matrix(NA_real_, ncol=ret.ncol, nrow=ret.nchunk)
	} else {
		ret.nchunk <- ncomb * nvariants
	}

	calls <- vector(mode="list", length = ret.nchunk)

	if(hasSubset <- !missing(subset))  {
		if(!tryCatch(is.language(subset), error = function(e) FALSE))
			subset <- substitute(subset)
		if(inherits(subset, "formula")) {
			if (subset[[1]] != "~" || length(subset) != 2L)
				stop("'subset' should be a one-sided formula")
			subset <- subset[[2L]]
		}
		if(!all(all.vars(subset) %in% all.terms))
			warning("Not all terms in 'subset' exist in 'global.model'")
	}

	comb.sfx <- rep(TRUE, n.fixed)
	comb.seq <- if(nov != 0L) seq.int(nov) else 0L
	k <- 0L
	ord <- integer(0L)

	#whichNotInt <- !(all.terms %in% interceptLabel)
	#nullComb <- rep(FALSE, length(all.terms))

	# TODO: has.int is not logical length(has.int) > 1

	for(j in seq.int(ncomb)) {
		#comb <- c(bitAnd(binPos, j - 1L) != 0L, rep(TRUE, n.fixed))
		comb <- c(as.logical(intToBits(j - 1L)[comb.seq]), comb.sfx)

		nvar <- sum(comb) - nInts
		if(nvar > m.max || nvar < m.min) next;
		#comb <- { v <- nullComb; v[whichNotInt] <- comb; v }

		if(hasSubset && !eval(subset, structure(as.list(comb), names=all.terms)))
			next;

		terms1 <- all.terms[comb]
		terms1[terms1 %in% interceptLabel] <- "1"
		#DebugPrint(terms1)

		frm <- reformulate(c("1", terms1), response=response, intercept=has.int)

		if(!formulaAllowed(frm, marg.ex)) next;
		attr(frm, ".Environment") <- env

		if (isMER) frm <- update.formula(frm, attr(all.terms, "random"))

		cl <- global.call
		cl[[formula.arg]] <- update.formula(global.formula, frm)
		###

		if(has.start) {
			clMmat <- cl
			clMmat[[1L]] <- as.name("model.matrix")
			clMmat[[formula.arg]] <- update.formula(global.formula, frm)
			names(clMmat)[2L] <- "object" # XXX
			cl$start <- cl$start[c(TRUE, globCoefNames0 %in% colnames(eval(clMmat)))]
		}

		for (ivar in seq.variants) { ## --- Variants ---------------------------
			clVariant <- cl
			if(nvarying) {
				updateArgs <- sapply(varying.names, function(x)
					varying[[x]][[variantsIdx[ivar, x]]], simplify=FALSE)
				for(i in varying.names) {

					clVariant[i] <- updateArgs[i]
				}
			}

			modelId <- ((j - 1L) * nvariants) + ivar
			if(trace) { cat(modelId, ": "); print(clVariant)	}

			if(evaluate) {
				fit1 <- tryCatch(eval(clVariant, parent.frame()), error=function(err) {
					err$message <- paste(conditionMessage(err), "(model",
						modelId, "skipped)", collapse="")
					class(err) <- c("simpleError", "warning", "condition")
					warning(err)
					return(NULL)
				})

				if (is.null(fit1)) next;

				k <- k + 1L # all OK, add model to table
				ord[k] <- modelId

				row1 <- c(
					matchCoef(fit1, all.terms=all.terms, beta=beta)[all.terms],
					if(nvarying) unlist(variantsIdx[ivar, ]),
					k=getK(fit1)
				)

				if (has.rsq) {
					fit1.summary <- summary(fit1)
					row1 <- c(row1, r.squared=fit1.summary$r.squared,
						adj.r.squared = fit1.summary$adj.r.squared)
				}
				if (has.dev)	row1 <- c(row1, dev=deviance(fit1))

				ic <- IC(fit1)
				row1 <- c(row1, IC = ic)

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

	if(!evaluate) return(calls[seq.int(k)])

	if(k < nrow(ret)) ret <- ret[seq.int(k), , drop=FALSE]

	ret <- as.data.frame(ret)
	#row.names(ret) <- seq_len(NROW(ret))
	row.names(ret) <- ord

	# Convert columns with presence/absence of terms to factors
	tfac <- which(!(all.terms %in% globCoefNames))

	ret[tfac] <- lapply(ret[tfac], factor, levels=NaN, labels="+")

	i <- seq_along(all.terms)
	v <- order(termsOrder)
	ret[, i] <- ret[, v]
	all.terms <- all.terms[v]

	colnames(ret) <- c(all.terms, varying.names, "k",
		if (has.rsq) c("R.sq", "Adj.R.sq"),
		if (has.dev) if (is.lm) "RSS" else "Dev.",
		#if (rank.custom) rank else c("AIC", "AICc")
		ICName
	)


	if(nvarying) {
		variant.names <- lapply(varying, function(x) make.unique(if(is.null(names(x))) as.character(x) else names(x)))
		for (i in varying.names)
			ret[, i] <-
			factor(ret[, i], levels = seq_along(variant.names[[i]]),
				labels = variant.names[[i]])
	}


	o <- order(ret[, ICName], decreasing = FALSE)
	ret <- ret[o, ]
	ret$delta <- ret[, ICName] - min(ret[, ICName])
	ret$weight <- exp(-ret$delta / 2) / sum(exp(-ret$delta / 2))

	class(ret) = c("model.selection", "data.frame")

	attr(ret, "calls") <- calls[o]
	attr(ret, "global") <- global.model
	attr(ret, "global.call") <- global.call
	attr(ret, "terms") <- all.terms

	attr(ret, "rank") <- IC
	attr(ret, "rank.call") <- attr(IC,"call")
	attr(ret, "call") <- match.call(expand.dots = TRUE)

	if (!is.null(attr(all.terms, "random.terms")))
		attr(ret, "random.terms") <- attr(all.terms, "random.terms")

	return(ret)
} ######



`subset.model.selection` <-
function(x, subset, select, recalc.weights = TRUE, ...) {
	if (missing(select)) {
		e <- substitute(subset)
		i <- eval(e, x, parent.frame(2L))
		return(`[.model.selection`(x, i, recalc.weights = recalc.weights, ...))
	} else {
		cl <- match.call(expand.dots = FALSE)
	    cl <- cl[c(1L, match(names(formals("subset.data.frame")), names(cl), 0L))]
	    cl[[1L]] <- as.name("subset.data.frame")
	    return(eval(cl, parent.frame(2L)))
	}
}


`[.model.selection` <-
function (x, i, j, recalc.weights = TRUE, ...) {
	res <- `[.data.frame`(x, i, j, ...)
	if (missing(j)) {
		att <- attributes(x)
		s <- c("row.names", "calls")
		att[s] <- lapply(att[s], `[`, i)
		attributes(res) <- att
		if(recalc.weights)
			res$weight <- res$weight / sum(res$weight)
	} else {
		cls <- class(res)
		class(res) <- cls[cls != "model.selection"] # numeric or data.frame
	}
	return(res)
}

`print.model.selection` <-
function(x, abbrev.names = TRUE, ...) {
	if(!is.null(x$weight))
		x$weight <- round(x$weight, 3)

	xterms <- attr(x, "terms")

	if(is.null(xterms)) {
		print.data.frame(x, ...)
	} else {

		xterms <- gsub(" ", "", xterms)

		if(abbrev.names) {
			xterms <- sapply(xterms, function(z) {
				spl <- strsplit(rep(z, 2L), c("[^\\w\\.]+",  "[\\w\\.$]+"),
					perl=TRUE)
				spl[[1]] <- abbreviate(spl[[1]], minlength=3)
				o <- order(sapply(lapply(spl, `==`, ""), `[[`, 1), decreasing=TRUE)
				x2 <- unsplit(spl, rep(order(sapply(lapply(spl, `==`, ""),
					`[[`, 1), decreasing=T), length.out=length(unlist(spl))))
				return(paste(x2, collapse=""))
			})
		}

		colnames(x)[seq_along(xterms)] <-  xterms


		cl <- attr(x, "global.call")
		if(!is.null(cl)) {
			cat("Global model: ")
			print(cl)
      cat("---\n")
		}

		cat ("Model selection table \n")
		i <- sapply(x, is.numeric)
		x[,i] <- signif(x[,i], 4L)
		print.default(as.matrix(x[, !sapply(x, function(.x) all(is.na(.x)))]),
					  na.print="", quote=FALSE)
		if (!is.null(attr(x, "random.terms"))) {
			cat("Random terms:", paste(attr(x, "random.terms"), collapse=", "),
				"\n")
		}
	}
}

`update.model.selection` <- function (object, ..., evaluate = TRUE) {
    call <- attr(object, "call")
    if (is.null(call))
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    return(if (evaluate) eval(call, parent.frame()) else call)
}
