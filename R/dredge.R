`dredge` <-
function(global.model, beta = FALSE, eval = TRUE, rank = "AICc",
		 fixed = NULL, m.max = NA, subset, marg.ex = NULL, trace = FALSE, ...) {

	# XXX
	if(!missing(eval) && !eval) .NotYetUsed("eval")
	rank.custom <- !missing(rank)
	#dots <- list(...)

	# *** Rank ***
	# switch to QAICc if quasi* family and no rank
	if(inherits(global.model, "glm") && family(global.model)$family %in%
		c("quasi", "quasibinomial", "quasipoisson") && !rank.custom) {
		rankArgs <- list(chat=summary(global.model)$dispersion)
		rank <- "QAICc"
		warning("QAICc used for '", family(global.model)$family,
				"' family with c-hat = ", signif(rankArgs$chat))
		rank.custom <- TRUE
	} else {
		rankArgs <- list(...)
	}
	IC <- .getRank(rank, rankArgs)
	ICName <- as.character(attr(IC, "call")[[1]])
	# *** Rank ***

	intercept <- "(Intercept)"
	all.terms <- getAllTerms(global.model)

	# Just in case:
	gterms <- tryCatch(terms(formula(global.model)),
		error=function(...) terms(global.model))

	response <- if(attr(gterms, "response") == 0) NULL else "."

	if(length(grep(":", all.vars(delete.response(gterms) > 0))))
		stop("Variable names in the model can not contain \":\"")

	global.call <- if(mode(global.model) == "S4") {
		if ("call" %in% slotNames(global.model)) slot(global.model, "call") else
			NULL
	} else {
		if(!is.null(global.model$call)) {
			global.model$call
		} else if(!is.null(attr(global.model, "call"))) {
			attr(global.model, "call")
		} else
			NULL
	}

	if (is.null(global.call)) {
		global.call <- substitute(global.model)
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

	is.glm <- inherits(global.model, "glm")
	is.lm <- !is.glm & inherits(global.model, "lm")


	if(isTRUE(rankArgs$REML) || (isTRUE(.isREMLFit(global.model)) && is.null(rankArgs$REML))) {
		warning("Comparing models with different fixed effects fitted by REML")
	}


	if (beta && is.null(tryCatch(beta.weights(global.model), error=function(e) NULL,
		warning=function(e) NULL))) {
		warning("Do not know how to calculate B-weights for ",
				class(global.model)[1L], ", argument 'beta' ignored")
         beta <- FALSE
	}

	summary.globmod <- summary(global.model)

	has.rsq <- is.list(summary.globmod) && is.numeric(summary.globmod$r.squared)
	has.dev <- !is.null(deviance(global.model))

	m.max <- if (missing(m.max)) n.vars else min(n.vars, m.max)

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


	#int.term <- if (has.int) "1" else "0"

	n.fixed <- length(fixed)
	all.terms <- do.call("structure", c(list(all.terms[order(all.terms %in%
		fixed)]), attributes(all.terms)))


	######

	#if (!eval) return(formulas)

	llik <- .getLogLik()
	getK <- function(x) as.vector(attr(llik(x), "df"))

	isMER <- any(inherits(global.model, c("mer", "lmer", "glmer")))
	env <- attr(gterms, ".Environment")


	# DEBUG <- function(x) cat("*", deparse(substitute(x)), "=", x, "*\n")

	### BEGIN:
	nov <- as.integer(n.vars - n.fixed)
	ncomb <- 2L ^ nov

	if(nov > 31L) stop(gettextf("Maximum number of varying predictors is 31, but %d is given", nov))
	if(nov > 10L) warning(gettextf("%d varying predictors will generate up to %d possible combinations", nov, ncomb))

	#binPos <- 2L ^ seq.int(0L, nov - 1L)

	ret.nchunk <- 25L
	ret.ncol <- length(all.terms) + (2L * has.rsq) + has.dev + 3L
	ret <- matrix(NA, ncol=ret.ncol, nrow=ret.nchunk)
	calls <- vector(mode="list", length=ret.nchunk)

	hasSubset <- !missing(subset)
	if(hasSubset) subset <- substitute(subset)

	comb.sfx <- rep(TRUE, n.fixed)
	comb.seq <- seq.int(nov)
	k <- 0L
	ord <- integer(0L)

	for(j in seq.int(ncomb)) {
		#comb <- c(bitAnd(binPos, j - 1L) != 0L, rep(TRUE, n.fixed))
		comb <- c(as.logical(intToBits(j - 1L)[comb.seq]), comb.sfx)
		if(sum(comb) > m.max) next;
		if(hasSubset && !eval(subset, structure(as.list(comb), names=all.terms))) next;

		terms1 <- all.terms[comb]
		frm <- reformulate(c("1", terms1), response=response, intercept=has.int)

		if(!formulaAllowed(frm, marg.ex)) next;
		attr(frm, ".Environment") <- env
		###

		row1 <- rep(NA, n.vars)
		row1[match(terms1, all.terms)] <- rep(1L, length(terms1))

		cl <- global.call
		if(has.start) {
			cl2 <- cl
			cl2[[1]] <- as.name("model.matrix")
			cl2[[formula.arg]] <- update.formula(global.formula, frm)
			names(cl2)[2] <- "object" # XXX
			cl$start <- cl$start[c(TRUE, globCoefNames0 %in% colnames(eval(cl2)))]
		}
		if (isMER) frm <- update.formula(frm, attr(all.terms, "random"))
		cl[[formula.arg]] <- update.formula(global.formula, frm)

		if(trace) { cat(j, ": "); print(cl)	}

		fit1 <- tryCatch(eval(cl, parent.frame()), error=function(err) {
			err$message <- paste(conditionMessage(err), "(model",
				j, "skipped by dredge)", collapse="")
			class(err) <- c("simpleError", "warning", "condition")

			warning(err)
			return(NULL)
		})

		if (is.null(fit1)) next;
		k <- k + 1L # all OK, add model to table
		ord[k] <- j

		coef1 <- if (beta) beta.weights(fit1)[, 3L] else coeffs(fit1)
		names(coef1) <- fixCoefNames(names(coef1))

		icept <- if (has.int) coef1[intercept] else NA

		coef1 <- coef1[c(na.omit(match(all.terms, names(coef1))))]
		row1[match(names(coef1), all.terms)] <- coef1

		row1 <- c(icept, row1, k=getK(fit1))
		if (has.rsq) {
			fit1.summary <- summary(fit1)
			row1 <- c(row1, r.squared=fit1.summary$r.squared,
				adj.r.squared=fit1.summary$adj.r.squared)
		}
		if (has.dev)	row1 <- c(row1, deviance(fit1))

		ic <- IC(fit1)
		row1 <- c(row1, IC=ic)

		ret.nrow <- nrow(ret)
		if(k > ret.nrow) {
			nadd <- min(ret.nchunk, ncomb - ret.nrow)
			ret <- rbind(ret, matrix(NA, ncol=ret.ncol, nrow=nadd))
			calls <- c(calls, vector("list", nadd))
		}
		calls[[k]] <- cl

		ret[k, ] <- row1
	} ### END

	if(k < nrow(ret)) ret <- ret[seq.int(k), , drop=FALSE]

	ret <- as.data.frame(ret)
	#row.names(ret) <- seq_len(NROW(ret))
	row.names(ret) <- ord

	# Convert columns with presence/absence of terms to factors
	tfac <- which(c(FALSE, !(all.terms %in% globCoefNames)))

	ret[tfac] <- lapply(ret[tfac], factor, levels=1, labels="+")

	colnames(ret) <- c("(int.)", all.terms, "k",
		if (has.rsq) c("R.sq", "Adj.R.sq"),
		if (has.dev) if (is.lm) "RSS" else "Dev.",
		#if (rank.custom) rank else c("AIC", "AICc")
		ICName
	)

	o <- order(ret[, ICName], decreasing = FALSE)
	ret <- ret[o,]
	ret$delta <- ret[, ICName] - min(ret[, ICName])
	ret$weight <- exp(-ret$delta / 2) / sum(exp(-ret$delta / 2))

	class(ret) = c("model.selection", "data.frame")

	attr(ret, "calls") <- calls[o]
	attr(ret, "global") <- global.model
	attr(ret, "global.call") <- global.call
	attr(ret, "terms") <- c(intercept, all.terms)

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
