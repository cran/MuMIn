`model.avg` <-
function(object, ..., beta = FALSE,
	rank = NULL, rank.args = NULL, revised.var = TRUE,
	dispersion = NULL) {

	if (isTRUE("method" %in% names(match.call())))
		stop("the argument 'method' is no longer accepted")

	if(inherits(object, "model.selection")) {
		if(!("subset" %in% names(match.call(get.models))))
			warning("'subset' argument is missing. Using all models by default")
		object <- get.models(object, ...)
	}

	if (inherits(object, "list")) {
		models <- object
		object <- object[[1L]]
        if (!is.null(rank) || is.null(rank <- attr(models, "rank"))) {
            rank <- .getRank(rank, rank.args = rank.args, object = object)
      	}
	} else {
		models <- list(object, ...)
        rank <- .getRank(rank, rank.args = rank.args, object = object)
	}

	nmodels <- length(models)
	if(nmodels == 1L) stop("only one model supplied. Nothing to do")
	.checkModels(models)

    alpha <- 0.05
	.fnull <- function(...) return(NULL)
	ICname <- deparse(attr(rank,"call")[[1L]])

	allterms1 <- lapply(models, getAllTerms)
	all.terms <- unique(unlist(allterms1))
	# sort by level (main effects first)
	all.terms <- all.terms[order(vapply(gregexpr(":", all.terms),
		function(x) if(x[1L] == -1L) 0L else length(x), numeric(1L)), all.terms)]

	# all.model.names <- modelNames(models, asNumeric = FALSE,
		# withRandomTerms = FALSE, withFamily = FALSE)
	all.model.names <- .modelNames(allTerms = allterms1, uqTerms = all.terms)

	#if(is.null(names(models))) names(models) <- all.model.names

	# check if models are unique:
	mcoeffs <- lapply(models, coeffs)
	dup <- unique(sapply(mcoeffs, function(i) which(sapply(mcoeffs, identical, i))))
	dup <- dup[sapply(dup, length) > 1L]
	ndups <- length(dup)

	if (ndups > 0L) stop("models are not unique. Duplicates: ",
		prettyEnumStr(sapply(dup, paste, sep = "", collapse = " = "),
			quote = "'"))

	# workaround for different behavior of model.matrix with lme: data argument is required
	if(any(linherits(models, c(lme = TRUE)))) {
		model.matrix.lme <- function(object, data = object$data, ...)
			model.matrix.default(object, data = data, ...)
	}

	ic <- vapply(models, rank, numeric(1L))
	#dev <- if (!is.null(tryCatch(deviance(models[[1L]]), error = .fnull)))
		#vapply(models, deviance, numeric(1L)) else NA
	ll <- vapply(models, logLik, numeric(1L))
	delta <- ic - min(ic)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))
	model.order <- order(weight, decreasing = TRUE)

	#if(!is.null(dispersion)) dispersion <- rep(dispersion, length.out = nmodels)

	# DEBUG:
	# sapply(sapply(sapply(models[model.order], coef), names), paste, collapse="+")
	# sapply(models, function(x) paste(match(getAllTerms(x), all.terms), collapse="+"))

	# ----!!! From now on, everything MUST BE ORDERED by 'weight' !!!-----------
	mstab <- cbind(logLik = ll, IC = ic, Delta = delta, Weight = weight,
		deparse.level = 0)
	rownames(mstab) <- all.model.names
	mstab <- mstab[model.order, ]
	if(!is.null(dispersion)) {
		dispersion <- rep(dispersion, length.out = nmodels)[model.order]
		mstab <- cbind(mstab, Dispersion = dispersion)
	}

	weight <- mstab[, "Weight"] # has been sorted in table
	models <- models[model.order]

	allCoefNames <- unique(unlist(lapply(mcoeffs, names)))
	allCoefNames <- allCoefNames[order(vapply(gregexpr(":", allCoefNames),
		function(x) if(x[1L] == -1L) 0L else length(x), numeric(1L)), allCoefNames)]
	npar <- length(allCoefNames)

	if (beta)	response.sd <- sd(model.response(model.frame(object)))

	all.coef <- all.se <- all.df <- matrix(NA_real_, nrow = nmodels, ncol = npar)
	for (i in seq_len(nmodels)) {
		m <- models[[i]]
		coefmat <- coefTable(m, dispersion = dispersion[i])
		ncoef <- NROW(coefmat)
		if(ncoef > 0L) {
			#if(is.null(mdf) || !length(mdf)) mdf <- rep(NA, ncoef)
			if (beta) coefmat[, 1L:2L] <- coefmat[, 1L:2L] * sd(model.matrix(m)) / response.sd
			j <- match(rownames(coefmat), allCoefNames)
			all.coef[i, j] <- coefmat[, 1L]
			all.se[i, j] <- coefmat[, 2L]
			mdf <- coefDf(m)[!is.na(coeffs(m))]
			if(!is.null(mdf) && length(mdf)) all.df[i, j] <- mdf
		}
	}
	rownames(all.se) <- rownames(all.coef) <- rownames(mstab)

	# Benchmark: 3.7x faster
	#system.time(for(i in 1:10000) t(array(unlist(p), dim=c(length(all.terms),length(models)))))
	#system.time(for(i in 1:10000) do.call("rbind", p))

	vpresent <- do.call("rbind", lapply(models, function(x)
		all.terms %in% getAllTerms(x)))

	importance <- apply(weight * vpresent, 2L, sum)
	names(importance) <- all.terms
	importance <- sort(importance, decreasing = TRUE)

	missing.par <- is.na(all.coef)

	avg.model <- t(vapply(seq_along(allCoefNames),
		function(i) par.avg(
			all.coef[, i],
			se = all.se[, i],
			weight = weight,
			df = all.df[, i],
			alpha = alpha,
			revised.var = revised.var),
		double(5L)))

	avg.model[is.nan(avg.model)] <- NA


	dimnames(avg.model) <- list(allCoefNames, c("Estimate", "Std. Err.",
		"Adjusted SE", "Lower CI", "Upper CI"))

	colnames(all.coef) <- colnames(all.se) <- colnames(all.df) <- allCoefNames
    names(all.terms) <- seq_along(all.terms)

	colnames(mstab)[2L] <- ICname

	coef.shrink <- avg.model[, 1] *
		#colSums(weight * !is.na(all.coef))
		colSums(array(weight * as.double(!missing.par), dim = dim(all.coef)))

	#global.mm <- model.matrix(fit)
	#cnmmxx2 <- unique(unlist(lapply(gm, function(x) names(coef(x)))))
	#mmx <- gmm[, cnmmxx[match(colnames(gmm), cnmmxx, nomatch = 0)]]

	mmxs <- tryCatch(cbindDataFrameList(lapply(models, model.matrix)),
					 error = .fnull)

	# Far less efficient:
	#mmxs <- lapply(models, model.matrix)
	#mx <- mmxs[[1]];
	#for (i in mmxs[-1])
	#	mx <- cbind(mx, i[,!(colnames(i) %in% colnames(mx)), drop=FALSE])

	# residuals averaged (with brute force)
	rsd <- tryCatch(apply(vapply(models, residuals, residuals(object)), 1L,
		weighted.mean, w = weight), error = .fnull)
	trm <- tryCatch(terms(models[[1L]]),
			error = function(e) terms(formula(models[[1L]])))
	frm <- reformulate(all.terms,
				response = attr(trm, "variables")[-1L][[attr(trm, "response")]])

	#print(colSums(weight) * array(as.double(!is.na(all.coef)), dim=dim(all.coef)))

	ret <- list(
		summary = as.data.frame(mstab),
		coefficients = all.coef,
		se = all.se,
		#variable.codes2 = all.terms,
		term.codes = attr(all.model.names, "variables"),
		avg.model = avg.model,
		coef.shrinkage = coef.shrink,
		importance = importance,
		beta = beta,
		term.names = allCoefNames,
		x = mmxs,
		residuals = rsd,
		formula = frm,
		call = match.call(),
		dfs = all.df
	)

	attr(ret, "modelList") <- models
	#attr(ret, "method") <- method # c("Coefficient" = method, "Variance" = method.var)
	attr(ret, "beta") <- beta
	attr(ret, "revised.var") <- revised.var
	class(ret) <- "averaging"
	return(ret)
}

`coef.averaging` <-
function(object, full = FALSE, ...) if(full) object$coef.shrinkage else
	object$avg.model[, 1L]


# TODO: predict type="response" + average on response scale
`predict.averaging` <-
function(object, newdata = NULL, se.fit = FALSE, interval = NULL,
	type = c("link", "response"), full = TRUE, ...) {

	type <- match.arg(type)
	if (!missing(interval)) .NotYetUsed("interval", error = FALSE)

	models <- attr(object, "modelList")

	# Benchmark: vapply is ~4x faster
	#system.time(for(i in 1:1000) sapply(models, inherits, what="gam")) /
	#system.time(for(i in 1:1000) vapply(models, inherits, logical(1L), what="gam"))

	# If all models inherit from lm:
	if ((missing(se.fit) || !se.fit)
		&& all(linherits(models, c(gam = FALSE, lm = TRUE)))
		&& !any(is.na(object$coef.shrinkage))
		) {
		coeff <- coef(object, full = full)
		frm <- formula(object)

		tt <- delete.response(terms(frm))
		X <- object$x

		if (missing(newdata) || is.null(newdata)) {
			Xnew <- X
		} else {
			xlev <- unlist(unname(lapply(models, "[[", "xlevels")),
						   recursive = FALSE, use.names = TRUE)
			Xnew <- model.matrix(tt, data = newdata, xlev = xlev)
		}

		Xnew <- Xnew[, match(names(coeff), colnames(Xnew), nomatch = 0L)]
		ret <- (Xnew %*% coeff)[, 1L]

		#if (se.fit) {
		#	scale <- 1
		#	covmx <- solve(t(X) %*% X)
		#	se <- sqrt(diag(Xnew %*% covmx %*% t(Xnew))) * sqrt(scale)
		#	return(list(fit = y, se.fit = se))
		#}
	} else {
		# otherwise, use brute force:

		if(full == FALSE) warning("argument 'full' ignored")

		#pred <- if(!missing(newdata))
		#	lapply(models, predict, newdata = newdata, se.fit = se.fit,...) else
		#	lapply(models, predict, se.fit = se.fit, ...)

		cl <- as.list(match.call())
		cl[[1L]] <- as.name("predict")
		if("type" %in% names(cl)) cl$type <- "link"
		if(!missing(newdata)) cl$newdata <- as.name("newdata")
		#pred <- do.call("lapply", c(as.name("models"), cl[-2]))

		cl <- as.call(cl)
		pred <- lapply(models, function(x) {
			cl[[2L]] <- x
			tryCatch(eval(cl), error=function(e) e)
		})

		err <- sapply(pred, inherits, "condition")
		if (any(err)) {
			lapply(pred[err], warning)
			stop(sprintf(ngettext(sum(err), "'predict' for model %s caused error",
				"'predict' for models %s caused errors"),
				paste(sQuote(names(models[err])), collapse = ", ")))
		}


		if(all(sapply(pred, function(x) c("fit", "se.fit") %in% names(x)))) {
			fit <- do.call("cbind", lapply(pred, "[[", "fit"))
			se.fit <- do.call("cbind", lapply(pred, "[[", "se.fit"))

			#npar <- sapply(models, function(x) as.vector(attr(logLik(x), "df")))

			revised.var <- attr(object, "revised.var")

			apred <- sapply(seq(nrow(fit)), function(i)
				par.avg(fit[i, ], se.fit[i, ], weight = object$summary$Weight,
					df = NA, revised.var = revised.var))

			# TODO: ase!
			#no.ase <- all(is.na(object$avg.model[,3]))
			# if(no.ase) 2 else 3
			ret <- list(fit = apred[1L, ], se.fit = apred[2L, ])

			if (type == "response") {
				fam <- tryCatch(vapply(models, function(z)
					unlist(family(z)[c("family", "link")]), character(2L)),
								error = function(e) NULL)

				if(!is.null(fam)) {
					if(any(fam[,1] != fam[, -1L]))
					stop("Cannot calculate predictions on the response scale ",
						 "with models with different families or link functions")
					ret$se.fit <- ret$se.fit * abs(family(models[[1L]])$mu.eta(ret$fit))
					ret$fit <- family(models[[1L]])$linkinv(ret$fit)
				}
			}

		} else if (all(sapply(pred, is.numeric))) {
			ret <- apply(do.call("cbind", pred), 1L, weighted.mean,
				w = object$summary$Weight)
		} else {
			stop("'predict' method for the component models returned",
				 " a value in unrecognised format")
		}
	}
	return(ret)
}

`fitted.averaging` <-
function (object, ...) predict.averaging(object)

`model.matrix.averaging` <-
function (object, ...) object$x

`summary.averaging` <-
function (object, ...) {
	cf <- object$avg.model
	no.ase <- all(is.na(cf[, 3L]))
	z <- abs(cf[,1] / cf[, if(no.ase) 2L else 3L])
	pval <- 2 * pnorm(z, lower.tail = FALSE)
	coefmat <- cbind(cf[, if(no.ase) 1L:2L else 1L:3L],
					 `z value` = z,
					 `Pr(>|z|)` = zapsmall(pval))
	object$coefmat <- coefmat
	structure(object, class=c("summary.averaging", "averaging"))
}

`confint.averaging` <-
function (object, parm, level = 0.95, ...) {
    cf <- object$coefficients
    pnames <- colnames(cf)
    if (missing(parm)) parm <- pnames
		else if (is.numeric(parm)) parm <- pnames[parm]

    a2 <- 1 - level
    a <- a2 / 2
    se <- object$se
    wts <- object$summary$Weight
    dfs <- object$dfs
    ci <- t(sapply(parm, function(i)
		par.avg(cf[,i], se[,i], wts, object$dfs[, i], alpha = a2)))[, 4L:5L]
    pct <- stats:::format.perc(c(a, 1L - a), 3L)
	
	ci[is.na(object$coef.shrinkage), ] <- NA_real_
    colnames(ci) <- pct
    return(ci)
}

`print.summary.averaging` <-
function (x, digits = max(3L, getOption("digits") - 3L),
    signif.stars = getOption("show.signif.stars"), ...) {

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")

    cat("Component models:\n")
	print(round(as.matrix(x$summary), 2L), na.print = "")

	#cat("\nVariable names abbreviations:\n")
	cat("\nTerm codes:\n")
	print.default(x$term.codes, quote = FALSE)


	cat("\nModel-averaged coefficients: \n")
	if (nnotdef <- sum(is.na(x$coefmat[, 1L])))
		 cat("(", nnotdef, " not defined because of singularities in all ",
			"component models) \n", sep = "")

	#no.ase <- all(is.na(x$avg.model[,3]))

	hasPval <- TRUE

    printCoefmat(x$coefmat, P.values = hasPval, has.Pvalue = hasPval,
		digits = digits, signif.stars = signif.stars)

	#if (no.ase) cat("Confidence intervals are unadjusted \n")

	cat("\nFull model-averaged coefficients (with shrinkage):", "\n")
	printCoefmat(matrix(x$coef.shrink, nrow = 1L,
		dimnames = list("", x$term.names)), P.values = FALSE,
		has.Pvalue = FALSE, cs.ind = seq_along(x$term.names), tst.ind = NULL)

	cat("\nRelative variable importance:\n")
	print(round(x$importance, 2L))
}

`print.averaging` <-
function(x, ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat("Component models:", "\n")
	comp.names <- rownames(x$summary)
	comp.names[comp.names == ""] <- "null"
	cat(format(sQuote(comp.names), justify = "l"), fill = TRUE)
    cat("\nCoefficients:", "\n")
    print.default(x$avg.model[, 1L])
}

#`vcov.averaging` <- function (object, full = FALSE, ...) {
`vcov.averaging` <- function (object, ...) {
	full <- FALSE
	models <- attr(object, "modelList")
	vcovs <- lapply(lapply(models, vcov), as.matrix)
	names.all <- object$term.names
	nvars <- length(names.all)
	nvarseq <- seq(nvars)
	wts <- object$summary$Weight
	wts <- wts / sum(wts) # normalize just in case

	#vcov0 <- matrix(NA, nrow=nvars, ncol = nvars,	dimnames = list(names.all, names.all))
	vcov0 <- matrix(if(full) 0 else NA, nrow = nvars,
		ncol=nvars,	dimnames = list(names.all, names.all))

	vcovs2 <- lapply(vcovs, function(v) {
		i <- match(dimnames(v)[[1L]], names.all)
		vcov0[i,i] <- v
		return(vcov0)
	})
	b1 <- object$coefficients
	if(full) b1[is.na(b1)] <- 0

	avgb <- object$avg.model[, 1L]
	#avgb <- colSums(t(b1) * wts, na.rm=T)

	res <- sapply(nvarseq, function(c1) sapply(nvarseq, function(c2) {
		 weighted.mean(sapply(vcovs2, "[", c1, c2) + (b1[, c1] - avgb[c1]) *
		 (b1[, c2] - avgb[c2]), wts, na.rm = TRUE)
	}))
	dimnames(res) <- list(names.all, names.all)
	return(res)
}

`logLik.averaging` <- function (object, ...)
	return(structure(lapply(attr(object, "modelList"), .getLogLik()),
			  names = rownames(object$summary)))
