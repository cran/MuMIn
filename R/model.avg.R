`model.avg` <-
function(m1, ..., beta = FALSE, method = c("0", "NA"), rank = NULL,
	rank.args = NULL, alpha = 0.05, revised.var = TRUE) {

	method <- match.arg(method)

	.fnull <- function(...) return(NULL)

	if (!is.null(rank)) {
	   	rankFn <- match.fun(rank)
		rank.call <- as.call(c(as.name("rankFn"), as.symbol("x"), rank.args))
		rank <- substitute(rank)
	} else if (!is.null(attr(m1, "rank.call"))) {
		rank.call <- attr(m1, "rank.call")
		rank.args <- as.list(attr(m1, "rank.call"))[-(1L:2L)]
		rankFn <- match.fun(rank.call[[1L]])
		rank <- as.character(rank.call[[1L]])
	}

	if (inherits(m1, "list")) {
		models <- m1
		m1 <- models[[1L]]
	} else {
		models <- list(m1, ...)
	}

	if (!is.null(rank)) {
		IC <- function(x) eval(rank.call)
		res <- IC(m1)
  		if (!is.numeric(res) || length(res) != 1L)
			stop("'rank' should return numeric vector of length 1")
	} else {
		IC <- AICc
	}

	if (length(models) == 1L) stop("Only one model supplied. Nothing to do")

	#Try to find if all models are fitted to the same data
	m.resp <- lapply(models, function(x) formula(x)[[2L]])
	if(!all(vapply(m.resp[-1L], "==", logical(1), m.resp[[1L]])))
		stop("Response differs between models")

	m.data <- lapply(models, function(x) (if(mode(x) == "S4") `@` else `$`)
					 (x, "call")$data)
	m.nresid <-	vapply(models, nobs, numeric(1L))
	if(!all(m.data[-1L] == m.data[[1]]) || !all(m.nresid[-1L] == m.nresid[[1L]]))
		stop("Models were not all fitted to the same dataset")

	all.terms <- unique(unlist(lapply(models, getAllTerms)))
	all.terms <- all.terms[order(vapply(gregexpr(":", all.terms),
		function(x) if(x[1L] == -1) 0 else length(x), numeric(1L)), all.terms)]

	all.model.names <- vapply(models,
		function(x) paste(match(getAllTerms(x), all.terms), collapse="+"), character(1L))
	# check if models are unique:
	dup <- duplicated(all.model.names)
	if (any(dup)) {
  		dup <- table(all.model.names)
		dup <- seq(all.model.names)[all.model.names %in% names(dup[dup > 1])]
		stop("Models are not unique. Duplicates: ", paste(dup, collapse=", "))
	}

	# workaround for different behavior of model.matrix with lme: data argument is required
	if(any(vapply(models, inherits, logical(1L), what="lme"))) {
		model.matrix.lme <- function(object, data=object$data, ...)
			model.matrix.default(object, data=data, ...)
	}

	aicc <- vapply(models, IC, numeric(1))
	dev <- if (!is.null(tryCatch(deviance(models[[1L]]), error=.fnull)))
		vapply(models, deviance, numeric(1)) else NA
	delta <- aicc - min(aicc)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))
	model.order <- order(weight, decreasing=TRUE)

	# DEBUG:
	# sapply(sapply(sapply(models[model.order], coef), names), paste, collapse="+")
	# sapply(models, function(x) paste(match(getAllTerms(x), all.terms), collapse="+"))

	# !!! From now on, everything MUST BE SORTED by 'weight' !!!

	selection.table <- data.frame(
		Deviance = dev,	AICc = aicc, Delta = delta, Weight = weight,
		row.names = all.model.names
	)[model.order, ]

	weight <- selection.table$Weight # sorted in table
	models <- models[model.order]

	all.par <- unique(unlist(lapply(models, function(m) names(coeffs(m)))))
	all.par <- all.par[order(vapply(gregexpr(":", all.par),
		function(x) if(x[1L] == -1) 0 else length(x), numeric(1L)), all.par)]
	npar <- length(all.par)
	ac <- rep(0, length = npar)

	if (beta)
		response.sd <- sd(model.response(model.frame(m1)))

	mtable <- t(vapply(models, function(m) {
		m.tTable <- tTable(m)
		m.ncoef <- NROW(m.tTable)
		if(m.ncoef == 0L) { ##
			rep(NA, npar * 3L)
		} else {
			m.se <- m.tTable[,2L]
			# FIXED: below is wrong for Mixed models
			#m.df <- n - length(m.coef)
			m.df <- coefDf(m) # -> double(ncoef)/ NA / NULL
			if(is.null(m.df) || !length(m.df)) m.df <- rep(NA, m.ncoef)

			if (beta) {
				m.tTable[, 1L:2L] <- m.tTable[,1L:2L] * sd(model.matrix(m)) / response.sd
			}

			m.vars <- match(all.par, rownames(m.tTable))
			return(c(m.tTable[m.vars, 1L:2L], m.df[m.vars]))
		}
		}, structure(double(npar * 3L),
		names=c(paste(rep(c("Coef","SE", "DF"), each=npar), all.par, sep=".")))
	)) ### << mtable

	# mtable is already sorted by weigth
	all.coef <- mtable[, seq_len(npar)]
	all.se <- mtable[, npar + seq_len(npar)]
	all.df <- mtable[, 2L * npar + seq_len(npar)]
	##
	rownames(all.se) <- rownames(all.coef) <- rownames(selection.table)

	importance <- apply(weight * t(sapply(models,
		function(x) all.terms %in% getAllTerms(x))), 2L, sum)
	names(importance) <- all.terms
	importance <- sort(importance, decreasing=T)


	if (method == "0") {
		all.coef[is.na(all.coef)] <- 0
		all.se[is.na(all.se)] <- 0
	}

	#avg.model <- t(sapply(seq_along(all.par),
	#	function(i) par.avg(all.coef[,i], all.se[,i], weight, all.df, alpha)))

	avg.model <- t(vapply(seq_along(all.par),
		function(i) par.avg(all.coef[,i], all.se[,i], weight, all.df[,i],
			alpha, revised.var),
		structure(double(5), names=c("Coefficient", "SE", "Adjusted SE",
		"Lower CI", "Upper CI"))))


	all.coef[all.coef == 0] <- NA
	all.se[all.se == 0] <- NA
	colnames(all.coef) <- colnames(all.se) <- rownames(avg.model) <-  all.par
    names(all.terms) <- seq_along(all.terms)

	if (!is.null(rank))
		colnames(selection.table)[2L] <- as.character(rank)

	#global.mm <- model.matrix(fit)
	#cnmmxx2 <- unique(unlist(lapply(gm, function(x) names(coef(x)))))
	#mmx <- gmm[, cnmmxx[match(colnames(gmm), cnmmxx, nomatch = 0)]]

	mmxs <- tryCatch(cbindDataFrameList(lapply(models, model.matrix)),
					 error=.fnull)

	# Far less efficient:
	#mmxs <- lapply(models, model.matrix)
	#mx <- mmxs[[1]];
	#for (i in mmxs[-1])
	#	mx <- cbind(mx, i[,!(colnames(i) %in% colnames(mx)), drop=FALSE])

	# residuals averaged (with brute force)
	rsd <- tryCatch(apply(vapply(models, residuals, residuals(m1)), 1,
		weighted.mean, w=weight), error=.fnull)
	trm <- tryCatch(terms(models[[1]]),
			error=function(e) terms(formula(models[[1L]])))
	frm <- reformulate(all.terms,
				response = attr(trm, "variables")[-1L][[attr(trm, "response")]])

	ret <- list(
		summary = selection.table,
		coefficients = all.coef,
		se = all.se,
		variable.codes = all.terms,
		avg.model = avg.model,
		importance = importance,
		beta = beta,
		term.names = all.par,
		x = mmxs,
		residuals = rsd,
		formula = frm,
		call = match.call(),
		dfs=all.df
	)

	attr(ret, "mList") <- models
	attr(ret, "method") <- method
	attr(ret, "beta") <- beta
	attr(ret, "revised.var") <- revised.var
	class(ret) <- "averaging"
	return(ret)
}


if (!existsFunction("nobs"))
`nobs` <- function(object, ...) UseMethod("nobs")

`nobs.mer` <- function(object, ...) object@dims[["n"]]
`nobs.lme` <- `nobs.gls` <- function(object, ...) object$dims$N
`nobs.glmmML` <- function(object, ...) length(object$coefficients) + object$cluster.null.df
`nobs.default` <- function(object, ...) NROW(resid(object, ...))


`coefDf` <- function(x) UseMethod("coefDf")
`coefDf.lme` <- function(x) x$fixDF$X
`coefDf.mer` <- function(x) rep(NA, x@dims[["p"]])
`coefDf.gls` <- function(x) rep(x$dims$N - x$dims$p, x$dims$p)
`coefDf.default` <- function(x) rep(df.residual(x), length(coef(x)))


`coef.averaging` <-
function(object, ...) object$avg.model[,1]

`predict.averaging` <-
function(object, newdata = NULL, se.fit = NULL, interval = NULL,
	type = c("link", "response"), ...) {

	type <- match.arg(type)

	#if(("type" %in% names(match.call())) && type != "link") {
	#if(!missing("type") && type != "link") {
	#	warning("Only predictions on the link scale are allowed. Argument ",
	#			"'type' ignored")
	#	type <- "link"
	#}
	if (!missing(interval)) .NotYetUsed("interval", error = FALSE)

	models <- attr(object, "mList")

	# vapply is ~4x faster !
	#system.time(for(i in 1:1000) sapply(models, inherits, what="gam")) /
	#system.time(for(i in 1:1000) vapply(models, inherits, logical(1), what="gam"))

	# If all models inherit from lm:
	if ((missing(se.fit) || !se.fit)
		&& all(vapply(models, inherits, logical(1), what="lm"))
		&& !any(vapply(models, inherits, logical(1), what="gam"))
		) {
		coeff <- coef(object)
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

		Xnew <- Xnew[, match(names(coeff), colnames(Xnew), nomatch = 0)]
		ret <- (Xnew %*% coeff)[, 1]

		#if (se.fit) {
		#	covmx <- solve(t(X) %*% X)
		#	se <- sqrt(diag(Xnew %*% covmx %*% t(Xnew))) * sqrt(scale)
		#	return(list(fit = y, se.fit = se))
		#}
	} else {
		# otherwise, use brute force:
		if(attr(object, "method") == "NA")
			warning("Prediction assumes 'method' is \"0\"")

		#pred <- if(!missing(newdata))
		#	lapply(models, predict, newdata = newdata, se.fit = se.fit,...) else
		#	lapply(models, predict, se.fit = se.fit, ...)

		cl <- as.list(match.call())
		cl[[1]] <- as.name("predict")
		if("type" %in% names(cl)) cl$type <- "link"
		if(!missing(newdata)) cl$newdata <- as.name("newdata")
		pred <- do.call("lapply", c(as.name("models"), cl[-2]))

		if(all(sapply(pred, function(x) c("fit", "se.fit") %in% names(x)))) {
			fit <- do.call("cbind", lapply(pred, "[[", "fit"))
			se.fit <- do.call("cbind", lapply(pred, "[[", "se.fit"))
			#npar <- sapply(models, function(x) as.vector(attr(logLik(x), "df")))

			apred <- sapply(seq(nrow(fit)), function(i)
				par.avg(fit[i, ], se.fit[i, ], object$summary$Weight, NA,
						revised.var = attr(object, "revised.var")))

			# TODO: ase!
			#no.ase <- all(is.na(object$avg.model[,3]))
			# if(no.ase) 2 else 3
			ret <- list(fit=apred[1,], se.fit=apred[2, ])

			if (type == "response") {
				fam <- tryCatch(vapply(models, function(z)
					unlist(family(z)[c("family", "link")]), character(2)),
								error=function(e) NULL)

				if(!is.null(fam)) {
					if(any(fam[,1] != fam[,-1]))
					stop("Cannot calculate predictions of response scale",
						 "with models with different families or link functions")
					ret$se.fit <- ret$se.fit * abs(family(models[[1]])$mu.eta(ret$fit))
					ret$fit <- family(models[[1]])$linkinv(ret$fit)
				}
			}


		} else if (all(sapply(pred, is.numeric))) {
			ret <- apply(do.call("cbind", pred), 1, weighted.mean,
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
function (object, ...) print.averaging(object)

`print.averaging` <-
function(x, ...) {
	cat("\nModel summary:\n")
	print(round(as.matrix(x$summary), 2), na.print="")

	cat("\nVariables:\n")
	print(x$variable.codes, quote=F)

	cat("\nAveraged model parameters:\n")

	no.ase <- all(is.na(x$avg.model[,3]))

	print(signif(x$avg.model[, if(no.ase) -3 else TRUE], 3))
	if (no.ase) cat("Confidence intervals are unadjusted \n")

	cat("\nRelative variable importance:\n")
	print(round(x$importance, 2))
}

`vcov.averaging` <- function (object, ...) {
	models <- attr(object, "mList")
	vcovs <- lapply(lapply(models, vcov), as.matrix)
	names.all <- object$term.names
	nvars <- length(names.all)
	nvarseq <- seq(nvars)
	wts <- object$summary$Weight / sum(object$summary$Weight) # normalize just in case

	vcov0 <- matrix(if(attr(object, "method") == "NA") NA else 0, nrow=nvars, ncol=nvars,
		dimnames = list(names.all, names.all))
	vcovs2 <- lapply(vcovs, function(v) {
		i <- match(dimnames(v)[[1]], names.all)
		vcov0[i,i] <- v
		return(vcov0)
	})

	b1 <- object$coefficients
	avgb <- object$avg.model[,1]
	#avgb <- colSums(t(b1) * wts, na.rm=T)

	vcov3 <- sapply(nvarseq, function(c1) sapply(nvarseq, function(c2) {
		 weighted.mean(sapply(vcovs2, "[", c1, c2) + (b1[, c1] - avgb[c1]) *
		 (b1[, c2] - avgb[c2]), wts, na.rm = T)
	}))
	dimnames(vcov3) <- list(names.all, names.all)
	return(vcov3)
}

`logLik.averaging` <- function (object, ...) {
	return(structure(lapply(attr(object, "mList"), logLik),
			  names=rownames(object$summary)))
}
