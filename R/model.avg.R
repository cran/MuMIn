`model.avg` <-
function (object, ..., revised.var = TRUE) {
	if (isTRUE("method" %in% names(match.call())))
		stop("argument 'method' is no longer accepted")
	UseMethod("model.avg")
}


`model.avg.model.selection` <-
function(object, subset, fit = FALSE, ..., revised.var = TRUE) {

	if(!missing(subset)) {
		cl <- match.call()
		cl[[1L]] <- as.name("subset")
		names(cl)[2L] <- "x"
		object <- eval(cl[1L:3L], parent.frame())
	}
	if(fit || length(list(...))) {
		cl <- match.call()
		cl$fit <- NULL
		arg1 <- names(cl)[-(1L:2L)] %in% names(formals("model.avg.default"))
		cl1 <- cl[c(TRUE, TRUE, !arg1)]
		cl1[[1L]] <- as.name("get.models")
		cl2 <- cl[c(TRUE, TRUE, arg1)]
		cl2[[2L]] <- cl1
		cl2[[1L]] <- as.name("model.avg")
		#message("recreating the model objects")
		return(eval(cl2, parent.frame()))
	}

	ct <- attr(object, "coefTables")
	#coefNames <- if(!is.null(attr(object, "global")))
		#names(coeffs(attr(object, "global"))) else
	coefNames <- fixCoefNames(unique(unlist(lapply(ct, rownames),
		use.names = FALSE)), sort = TRUE)

	cfarr <- coefArray(ct)
	coefNames <- dimnames(cfarr)[[3]]
	nCoef <- length(coefNames)
	nModels <- length(ct)
	weight <- object$weight / sum(object$weight)

	avgcoef <- array(dim = c(nCoef, 5L), dimnames = list(coefNames, c("Estimate",
		"Std. Error", "Adjusted SE", "Lower CI", "Upper CI")))
	for(i in seq_len(nCoef))
		avgcoef[i, ] <- par.avg(cfarr[, 1L, i], cfarr[, 2L, i], weight,
			df = cfarr[, 3L, i], revised.var = revised.var)

	missing.par <- is.na(cfarr[, 1L, ])
	coef.shrink <- avgcoef[, 1L] *
		colSums(array(weight * as.double(!missing.par), dim = c(nModels, nCoef)))

	#allterms1 <- lapply(attr(object, "calls"), function(x)
		#getAllTerms(as.formula(x[[switch(as.character(x[[1L]]),
			#lme=, lme.formula= "fixed", gls= "model", "formula")]])))
	all.terms <- attr(object, "terms")
	all.vterms <- all.terms[!(all.terms %in% attr(all.terms, "interceptLabel")
		| apply(is.na(object[, all.terms]), 2L, all))]
	#allterms1 <- apply(!is.na(object[, all.vterms, drop = FALSE]), 1L, function(x) all.vterms[x])
	allterms1 <- applyrns(!is.na(object[, all.vterms, drop = FALSE]), function(x) all.vterms[x])
	all.model.names <- .modelNames(allTerms = allterms1, uqTerms = all.vterms)

	mstab <- object[, -(seq_len(ncol(object) - 5L))]
	colnames(mstab)[4:5] <- c("Delta", "Weight")
	rownames(mstab) <- all.model.names

	avgcoef[is.nan(avgcoef)] <- NA

	ret <- list(
		summary = as.data.frame(mstab),
		term.codes = attr(all.model.names, "variables"),
		avg.model = avgcoef,
		coef.shrinkage = coef.shrink,
		coefArray = cfarr,
		importance = importance(object),
		beta = attr(object, "beta"),
		term.names = coefNames,
		x = NULL,
		residuals = NULL,
		formula = if(!is.null(attr(object, "global")))
			formula(attr(object, "global")) else NULL,
		call = match.call()
	)

	attr(ret, "beta") <- attr(object, "beta")
	attr(ret, "nobs") <- attr(object, "nobs")
	attr(ret, "revised.var") <- revised.var
	class(ret) <- "averaging"
	return(ret)
}


`model.avg.default` <-
function(object, ..., beta = FALSE,
	rank = NULL, rank.args = NULL, revised.var = TRUE,
	dispersion = NULL, ct.args = NULL) {

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

	nModels <- length(models)
	if(nModels == 1L) stop("only one model supplied. Nothing to do")
	.checkModels(models)

    alpha <- 0.05
	.fnull <- function(...) return(NULL)
	ICname <- deparse(attr(rank, "call")[[1L]])

	allterms1 <- lapply(models, getAllTerms)
	all.terms <- unique(unlist(allterms1, use.names = FALSE))

	# sort by level (main effects first)
	all.terms <- all.terms[order(vapply(gregexpr(":", all.terms),
		function(x) if(x[1L] == -1L) 0L else length(x), numeric(1L)), all.terms)]

	# all.model.names <- modelNames(models, asNumeric = FALSE,
		# withRandomTerms = FALSE, withFamily = FALSE)
	all.model.names <- .modelNames(allTerms = allterms1, uqTerms = all.terms)

	#if(is.null(names(models))) names(models) <- all.model.names
	
	coefTableCall <- as.call(c(alist(coefTable, models[[i]],
		dispersion = dispersion[i]), ct.args))

	# check if models are unique:
	if(!is.null(dispersion)) dispersion <- rep(dispersion, length.out = nModels)
	coefTables <- vector(nModels, mode = "list")
	for(i in seq_len(nModels)) coefTables[[i]] <-
		eval(coefTableCall)
		#coefTable(models[[i]], dispersion = dispersion[i])
	mcoeffs <- lapply(coefTables, "[", , 1L)

	dup <- unique(sapply(mcoeffs, function(i) which(sapply(mcoeffs, identical, i))))
	dup <- dup[sapply(dup, length) > 1L]
	if (length(dup) > 0L) stop("models are not unique. Duplicates: ",
		prettyEnumStr(sapply(dup, paste, sep = "", collapse = " = "),
			quote = "'"))

	# workaround for different behaviour of model.matrix with lme: data argument is required
	if(any(linherits(models, c(lme = TRUE)))) {
		model.matrix.lme <- function(object, data = object$data, ...)
			model.matrix.default(object, data = data, ...)
	}
	
	LL <- .getLik(object)
	logLik <- LL$logLik
	lLName <- LL$name

	ic <- vapply(models, rank, numeric(1L))
	logLiks <- lapply(models, logLik)
	delta <- ic - min(ic)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))
	model.order <- order(weight, decreasing = TRUE)

	# ----!!! From now on, everything MUST BE ORDERED by 'weight' !!!-----------
	mstab <- cbind(df = vapply(logLiks, attr, numeric(1L), "df"),
		logLik = as.numeric(logLiks), IC = ic, Delta = delta, Weight = weight,
		deparse.level = 0L)
	if(!is.null(dispersion)) mstab <- cbind(mstab, Dispersion = dispersion)
	rownames(mstab) <- all.model.names
	mstab <- mstab[model.order, ]
	weight <- mstab[, "Weight"] # has been sorted in table
	models <- models[model.order]
	coefTables <- coefTables[model.order]

	if (beta) {
		response.sd <- sd(model.response(model.frame(object)))
		for(i in seq_along(coefTables))
			coefTables[[i]][, 1L:2L] <-
				coefTables[[i]][, 1L:2L] *
				apply(model.matrix(models[[i]]), 2L, sd) / response.sd
	}

	cfarr <- coefArray(coefTables)
	coefNames <- dimnames(cfarr)[[3]]
	nCoef <- length(coefNames)

	# Benchmark: 3.7x faster
	#system.time(for(i in 1:10000) t(array(unlist(p), dim=c(length(all.terms),length(models)))))
	#system.time(for(i in 1:10000) do.call("rbind", p))

	vpresent <- do.call("rbind", lapply(models, function(x)
		all.terms %in% getAllTerms(x)))

	importance <- apply(weight * vpresent, 2L, sum)
	names(importance) <- all.terms
	importance <- sort(importance, decreasing = TRUE)

	avgcoef <- t(vapply(seq_along(coefNames), function(i) par.avg(
			cfarr[, 1L, i], se = cfarr[, 2L, i], df = cfarr[, 3L, i],
			weight = weight, alpha = alpha,
			revised.var = revised.var), double(5L)))

	avgcoef[is.nan(avgcoef)] <- NA

	dimnames(avgcoef) <- list(coefNames, c("Estimate", "Std. Error",
		"Adjusted SE", "Lower CI", "Upper CI"))

    names(all.terms) <- seq_along(all.terms)
	colnames(mstab)[3L] <- ICname

	missing.par <- is.na(cfarr[, 1L, ])
	coef.shrink <- avgcoef[, 1L] *
		#colSums(weight * !is.na(all.coef))
		colSums(array(weight * as.double(!missing.par), dim = c(nModels, nCoef)))

	#global.mm <- model.matrix(fit)
	#cnmmxx2 <- unique(unlist(lapply(gm, function(x) names(coef(x)))))
	#mmx <- gmm[, cnmmxx[match(colnames(gmm), cnmmxx, nomatch = 0)]]

	mmxs <- tryCatch(cbindDataFrameList(lapply(models, model.matrix)),
					 error = .fnull, warning = .fnull)

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

	ret <- list(
		summary = as.data.frame(mstab),
		term.codes = attr(all.model.names, "variables"),
		avg.model = avgcoef,
		coef.shrinkage = coef.shrink,
		coefArray = cfarr,
		importance = importance,
		beta = beta,
		term.names = coefNames,
		x = mmxs,
		residuals = rsd,
		formula = frm,
		call = match.call()
	)

	attr(ret, "modelList") <- models
	attr(ret, "beta") <- beta
	attr(ret, "nobs") <- nobs(object)
	attr(ret, "revised.var") <- revised.var
	class(ret) <- "averaging"
	return(ret)
}

`coef.averaging` <-
function(object, full = FALSE, ...) if(full) object$coef.shrinkage else
	object$avg.model[, 1L]



 #TODO: predict p -="response" + average on response scale
`predict.averaging` <-
function(object, newdata = NULL, se.fit = FALSE, interval = NULL,
	type = NA, backtransform = FALSE,
	full = TRUE, ...) {

	#if(is.null(type)) type <- "link"
	if (!missing(interval)) .NotYetUsed("interval", error = FALSE)

	models <- attr(object, "modelList")
	if(is.null(models)) stop("cannot predict without a model list")

	# Benchmark: vapply is ~4x faster
	#system.time(for(i in 1:1000) sapply(models, inherits, what="gam")) /
	#system.time(for(i in 1:1000) vapply(models, inherits, logical(1L), what="gam"))

	# If all models inherit from lm:
	if ((missing(se.fit) || !se.fit)
		&& (is.na(type) || type == "link")
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

		cl <- as.list(match.call())
		cl$backtransform <- cl$full <- NULL
		cl[[1L]] <- as.name("predict")
		#if("type" %in% names(cl)) cl$type <- "link"
		if(!missing(newdata)) cl$newdata <- as.name("newdata")
		cl <- as.call(cl)

		#predict.lme <- MuMIn:::predict.lme

		pred <- lapply(models, function(x) {
			cl[[2L]] <- x
			y <- tryCatch({
				y <- eval(cl)
				if(is.numeric(y)) y else structure(as.list(y[c(1L, 2L)]),
					names = c("fit", "se.fit"))
				}, error = function(e) e)
		})

		err <- sapply(pred, inherits, "condition")
		if (any(err)) {
			lapply(pred[err], warning)
			stop(sprintf(ngettext(sum(err), "'predict' for model %s caused error",
				"'predict' for models %s caused errors"),
					paste(sQuote(names(models[err])), collapse = ", ")))
		}

		.untransform <- function(fit, se.fit = NULL, models) {
			fam <- tryCatch(vapply(models, function(z)
				unlist(family(z)[c("family", "link")]), character(2L)),
							error = function(e) NULL)
			if (!is.null(fam)) {
				if(any(fam[, 1L] != fam[, -1L]))
				stop("Cannot calculate prediction on a response scale ",
					 "with models using different families or link functions")
				fam1 <- family(models[[1L]])
				if(is.null(se.fit))
					return(fam1$linkinv(fit))
				else return(list(fit = fam1$linkinv(fit),
					se.fit = se.fit * abs(fam1$mu.eta(fit))
					))
			}
			return(NULL)
		}

		if (all(sapply(pred, is.list))) {
		#if(all(sapply(pred, function(x) c("fit", "se.fit") %in% names(x)))) {
			fit <- do.call("cbind", lapply(pred, "[[", "fit"))
			se.fit <- do.call("cbind", lapply(pred, "[[", "se.fit"))
			revised.var <- attr(object, "revised.var")

			apred <- unname(vapply(seq(nrow(fit)), function(i)
				par.avg(fit[i, ], se.fit[i, ], weight = object$summary$Weight,
					df = NA_integer_, revised.var = revised.var),
					FUN.VALUE = double(5L)))

			# TODO: ase!
			#no.ase <- all(is.na(object$avg.model[,3]))
			# if(no.ase) 2 else 3
			ret <-  if (backtransform)
					.untransform(apred[1L, ], apred[2L, ], models = models)
				else
					list(fit = apred[1L, ], se.fit = apred[2L, ])
		} else {
			#tryCatch({
			i <- !vapply(pred, is.numeric, FALSE)
			if(any(i)) pred[i] <- lapply(pred[i], "[[", 1L)
			ret <- apply(do.call("cbind", pred), 1L, weighted.mean,
				w = object$summary$Weight)
			if (backtransform)
				ret <- .untransform(ret, models = models)
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
	object$coefmat <- cbind(cf[, if(no.ase) 1L:2L else 1L:3L],
					 `z value` = z,
					 `Pr(>|z|)` = zapsmall(pval))
	structure(object, class = c("summary.averaging", "averaging"))
}

`confint.averaging` <-
function (object, parm, level = 0.95, ...) {
    cf <- object$coefArray[, 1L, ]
    pnames <- colnames(cf)
    if (missing(parm)) parm <- pnames
		else if (is.numeric(parm)) parm <- pnames[parm]

    a2 <- 1 - level
    a <- a2 / 2
    se <- object$coefArray[, 2L, ]
    wts <- object$summary$Weight
    dfs <- object$coefArray[, 3L, ]
    ci <- t(sapply(parm, function(i)
		par.avg(cf[,i], se[,i], wts, dfs[, i], alpha = a2)))[, 4L:5L]
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
	printCoefmat(matrix(x$coef.shrinkage, nrow = 1L,
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
	if(is.null(models)) stop("cannot calculate covariance matrix without a model list")

	vcovs <- lapply(lapply(models, vcov), as.matrix)
	names.all <- object$term.names
	nvars <- length(names.all)
	nvarseq <- seq(nvars)
	wts <- object$summary$Weight
	wts <- wts / sum(wts) # normalize just in case

	#vcov0 <- matrix(NA, nrow=nvars, ncol = nvars,	dimnames = list(names.all, names.all))
	vcov0 <- matrix(if(full) 0 else NA, nrow = nvars,
		ncol = nvars, dimnames = list(names.all, names.all))

	v <- vcovs[[1]]

	vcovs2 <- lapply(vcovs, function(v) {
		i <- match(dimnames(v)[[1L]], names.all)
		vcov0[i, i] <- v
		return(vcov0)
	})
	b1 <- object$coefArray[, 1L, ]
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

`logLik.averaging` <- function (object, ...) {
	models <- attr(object, "modelList")
	if(is.null(models)) {
		nobs <- attr(object, "nobs")
		apply(object$summary, 1L, function(x) structure(list(x[2L]),
			df = x[1L], nobs = nobs, class = "logLik"))
	} else {
		structure(lapply(attr(object, "modelList"), .getLogLik()),
			names = rownames(object$summary))
	}
}

`coefTable.averaging` <- function (model, ...)  model$avg.model[, 1:2]
