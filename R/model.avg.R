`model.avg` <-
function(m1, ..., beta = FALSE, method = c("0", "NA"), rank = NULL,
	rank.args = NULL, alpha = 0.05) {

	method <- match.arg(method)

	.fnull <- function(...) return(NULL)

	if (!is.null(rank)) {
	   	rankFn <- match.fun(rank)
		rank.call <- as.call(c(as.name("rankFn"), as.symbol("x"), rank.args))
		rank <- substitute(rank)
	} else if (!is.null(attr(m1, "rank.call"))) {
		rank.call <- attr(m1, "rank.call")
		rank.args <- as.list(attr(m1, "rank.call"))[-(1:2)]
		rankFn <- match.fun(rank.call[[1]])
		rank <- as.character(rank.call[[1]])
	}

	if (inherits(m1, "list")) {
		models <- m1
		m1 <- models[[1]]
	} else {
		models <- list(m1, ...)
	}

	if (!is.null(rank)) {
		IC <- function(x) eval(rank.call)
		res <- IC(m1)
  		if (!is.numeric(res) || length(res) != 1)
			stop("'rank' should return numeric vector of length 1")
	} else {
		IC <- AICc
	}

	if (length(models) == 1) stop("Only one model supplied. Nothing to do")

	#Try to find if all models are fitted to the same data
	m.resp <- lapply(models, function(x) formula(x)[[2]])
	if(!all(sapply(m.resp[-1], "==", m.resp[[1]])))
		stop("Response differs between models")

	m.data <- lapply(models, function(x) (if(mode(x) == "S4") `@` else `$`)
					 (x, "call")$data)
	m.nresid <- sapply(models, function(x) length(resid(x)))
	if(!all(m.data[-1] == m.data[[1]]) || !all(m.nresid[-1] == m.nresid[[1]]))
		stop("Models were not fitted to the same data")

	all.terms <- unique(unlist(lapply(models, getAllTerms)))
	all.terms <- all.terms[order(sapply(gregexpr(":", all.terms),
		function(x) if(x[1] == -1) 0 else length(x)), all.terms)]

	all.model.names <- sapply(models,
		function(x) paste(match(getAllTerms(x), all.terms), collapse="+"))
	# check if models are unique:
	dup <- duplicated(all.model.names)
	if (any(dup)) {
  		dup <- table(all.model.names)
		dup <- seq(all.model.names)[all.model.names %in% names(dup[dup > 1])]
		stop("Models are not unique. Duplicates: ", paste(dup, collapse=", "))
	}

	# workaround for different behavior of model.matrix with lme: data argument is required
	if(any(sapply(models, inherits, "lme"))) {
		model.matrix.lme <- function(object, data=object$data, ...)
			model.matrix.default(object, data=data, ...)
	}

	aicc <- sapply(models, IC)
	dev <- if (!is.null(tryCatch(deviance(models[[1]]), error=.fnull)))
		sapply (models, deviance) else NA
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
	all.par <- all.par[order(sapply(gregexpr(":", all.par),
		function(x) if(x[1] == -1) 0 else length(x)), all.par)]
	npar <- length(all.par)
	ac <- rep(0, length = npar)

	mtable <- t(sapply(models, function(m) {
		m.tTable <- tTable(m)
		n <- NROW(residuals(m))
		m.coef <- m.tTable[,1]
		m.var <- m.tTable[,2]
 		m.df <- n - length(m.coef)
		if (beta) {
			response.sd <- sd(model.frame(m1)[, attr(terms(m1), "response")])
			m.vars.sd <- sd(model.matrix(m))
			bx <- m.vars.sd / response.sd
			m.coef <- m.coef * bx
			m.var <- m.var * bx
		}
		m.vars <- match(all.par, rownames(m.tTable))
		return(c(coef=m.coef[m.vars], var=m.var[m.vars], df=m.df))
	}))

	# mtable is already sorted by weigth
	all.coef <- mtable[, 1:npar]
	all.var <- mtable[, npar + (1:npar)]
	all.df <- mtable[, 2 * npar + 1]
	##
	rownames(all.var) <- rownames(all.coef) <- rownames(selection.table)

	importance <- apply(weight * t(sapply(models,
		function(x) all.terms %in% getAllTerms(x))), 2, sum)
	names(importance) <- all.terms
	importance <- sort(importance, decreasing=T)


	if (method == "0") {
		all.coef[is.na(all.coef)] <- 0
		all.var[is.na(all.var)] <- 0
	}

	avg.model <- t(sapply(seq_along(all.par),
		function(i) par.avg(all.coef[,i], all.var[,i], all.df, weight, alpha)))
	all.coef[all.coef == 0] <- NA
	all.var[all.var == 0] <- NA
	colnames(all.coef) <- colnames(all.var) <- rownames(avg.model) <-  all.par
    names(all.terms) <- seq_along(all.terms)

	if (!is.null(rank))
		colnames(selection.table)[2] <- as.character(rank)

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
	rsd <- tryCatch(apply(sapply(models, residuals), 1, weighted.mean, w=weight),
					error=.fnull)
	trm <- tryCatch(terms(models[[1]]),
			error=function(e) terms(formula(models[[1]])))
	frm <- reformulate(all.terms,
				response = attr(trm, "variables")[-1][[attr(trm, "response")]])

	ret <- list(
		summary = selection.table,
		coefficients = all.coef,
		variable.codes = all.terms,
		variance = all.var,
		avg.model = avg.model,
		relative.importance = importance,
		weights = weight,
		beta = beta,
		term.names = all.par,
		x = mmxs,
		residuals = rsd,
		formula = frm,
		method = method,
		call = match.call()
	)

	attr(ret, "mList") <- models
	class(ret) <- "averaging"
	return(ret)
}

`coef.averaging` <-
function(object, ...) object$avg.model[,1]

`predict.averaging` <-
function(object, newdata = NULL, se.fit = NULL, interval = NULL, type = NULL,
	...) {

	#if(("type" %in% names(match.call())) && type != "link") {
	if(!missing("type") && type != "link") {
		warning("Only predictions on the link scale are allowed. Argument ",
				"'type' ignored")
	}
	if (!missing(se.fit)) .NotYetUsed("se.fit", error = FALSE)
	if (!missing(interval)) .NotYetUsed("interval", error = FALSE)

	models <- attr(object, "mList")

	# If all models inherit from lm:
	if (all(sapply(models, inherits, what="lm"))
		&& !any(sapply(models, inherits, what="gam"))
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

		Xnew <- Xnew[, match(names(coeff),colnames(Xnew), nomatch = 0)]
		ny <- (Xnew %*% coeff)[, 1]

		#if (se.fit) {
		#	covmx <- solve(t(X) %*% X)
		#	se <- sqrt(diag(Xnew %*% covmx %*% t(Xnew))) * sqrt(scale)
		#	return(list(fit = y, se.fit = se))
		#}
	} else {
		# otherwise, use brute force:
		if(object$method == "NA")
			warning("Prediction for this type of model assumes 'method' is \"0\"")

		ny <- if(!missing(newdata))
			sapply(models, predict, newdata = newdata, ...)
		else
			sapply(models, predict, ...)

		ny <- apply(ny, 1, weighted.mean, w = object$weight)
	}

	return(ny)
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
	print(signif(x$avg.model, 3))

	cat("\nRelative variable importance:\n")
	print(round(x$relative.importance, 2))
}
