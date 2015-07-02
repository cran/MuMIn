## EXCLUDE

std_predict <-
function (object, ...) UseMethod("std_predict")

invtransform <-
function (object, y) UseMethod("invtransform")

terms.averaging <-
function (x, ...) {
	terms(formula(x))
}

predict_avg <-
function(x, newdata, type = c("response", "invlink", "link", "terms"),
		 se.fit = FALSE, full =  TRUE, linkinv = NULL,
		 use.lincomb = FALSE, ...) {
	type <- match.arg(type)
	if(missing(newdata)) newdata <- NULL
	models <- getModelList(x)
	
	if((type != "terms" && !full) || use.lincomb) {
		#stop("averaged prediction via linear combination is not implemented yet")
		beta <- coef(x, full = full)
		models <- attr(x, "modelList")
		mf <- mergeMF(models)
		contr <- mergeContrasts(models)
		offwts <- offsetWeights(Weights(x), terms(mf), models)
		y <- lincomb(beta, mf, newdata, contr, offwts)
		if(type == "invlink") {
			if(is.null(linkinv)) stop("argument 'linkinv' is required")
			if(inherits(linkinv, c("family", "link-glm",  "character")))
				linkinv <- glm.link(linkinv)$linkinv
			y <- linkinv(y)
		}
	} else {
		xtype <- if(type == "invlink") "link" else type
		yall <- vector("list", length(models))
		for(i in seq_along(models))
			yall[[i]] <- std_predict(models[[i]], newdata = newdata,
				type = xtype, se.fit = se.fit, ...)
		
		y1 <- yall[[1L]]
		mod1 <- models[[1L]]
		if(!all(vapply(models, inherits, FALSE, class(mod1))))
			stop("cannot average predictions from different component model classes (yet)")
			
		if(type == "invlink") {
			if(!is.null(linkinv)) warning("argument 'linkinv' ignored")
				links <- tryCatch(vapply(models, function(m) glm.link(m)$name, ""),
				error = function(e) NULL)
			if (is.null(links)) stop("link functions of at least some component models not recognized")
			if(any(links[1L] != links[-1L]))
				stop("inverse transformation not possible because component models use different link functions")
		}
		dfs <- NULL
		fn <- if(is.matrix(y1)) {
			#if(all(attr(terms(mod1),"term.labels") == colnames(y1)))
			if(type == "terms") {
				#dfs <- tryCatch(vapply(models, df.residual, 1), error = function(e) NULL)
				avgterms
			} else avgmat
		} else if (is.list(y1) && all(c("fit", "se.fit") %in% names(y1))) {
			#dfs <- tryCatch(vapply(models, df.residual, 1), error = function(e) NULL)
			if(is.matrix(y1$fit)) 	#if(type == "terms" && is.matrix(y1$fit)) {
				avgterms else avgsefit
		} else if(is.numeric(y1)) {
			avgvec
		} else stop("'predict' on component model returned value in unsupported format")
		y <- fn(yall, Weights(x), revised.var = attr(x, "revised.var"),
				full = full, dfs = dfs)
		if(type == "invlink") y <- invtransform(mod1, y)
	}
	y
}

avg.predictions <-
function(yall, w, type = c("atomic", "matrix", "se.fit", "terms"),
		 revised.var = TRUE, full = FALSE, ...) {
	type <- if(is.null(type)) {
		} else match.arg(type)
	switch(type,
		matrix = avgmat(yall, w, ...),
		atomic = avgvec(yall, w, ...),
		se.fit = avgsefit(yall, w, revised.var, full, NULL, ...),
		terms = avgterms(yall, w, revised.var, full, NULL, ...)
		)
}

avgvec <-
function(yall, w, ...) {
	ret <- yall[[1L]]
	ydim <- c(length(ret), length(yall))
	yarr <- array(unlist(yall), dim = ydim)
	for(i in seq_len(ydim[1L]))
		ret[i] <- weighted.mean(yarr[i, ], w = w)
	ret
}

## helper function
toarray <- function(x)
	array(unlist(x), dim = c(dim(x[[1L]]), length(x)),
	 dimnames = c(dimnames(x[[1L]]),  list(names(x))))

issefit <-
function(x)
is.list(x) && (length(x) >= 2L) && all(names(x)[1L:2L] == c("fit", "se.fit"))

predtypes <-
function(cl, types, argname = "type") {
	ntypes <- length(types)
	y <- vector("list", ntypes)
	names(y) <- types <- types
	for(i in seq_len(ntypes)) {
		cl[[argname]] <- types[i]
		y[[i]] <- eval.parent(cl, n = 2L)
	}
	y
}

## averages prediction of each model term (predict.[g]lm):
avgterms <-
function(yall, w, revised.var, full, dfs, ...) {
	fitonly <- vapply(yall, is.atomic, FALSE)
	stopifnot(!any(fitonly[1L] != fitonly[-1L])) ## inconsistent prediction form
	fitonly <- fitonly[1L]

	if(fitonly) {
		all.coefnames <- unique(fixCoefNames(unlist(lapply(yall, colnames))))
		fits <- toarray(lapply(yall, function(x) x[, match(all.coefnames, colnames(x)), drop = FALSE]))
		if(full) fits[is.na(fits)] <- 0
		sefits <- array(0, dim = dim(fits))
	} else {
		all.coefnames <- unique(fixCoefNames(unlist(lapply(yall, function(x) colnames(x$fit)))))
		yall <- lapply(yall, function(x) {	
			j <- match(all.coefnames, colnames(x$fit))
			x$fit <- x$fit[, j]
			x$se.fit <- x$se.fit[, j]
			x
		})
		fits <- toarray(lapply(yall, "[[", "fit"))
		sefits <- toarray(lapply(yall, "[[", "se.fit"))
		
		if(full) {
			sefits[is.na(sefits)] <- 0
			fits[is.na(fits)] <- 0
		}
	}
	
	hasDfs <- !is.null(dfs) && !anyNA(dfs)
	k <- if(hasDfs) 3L else 2L
	d <- dim(fits)
	avgfit <- avgsefit <- array(NA_real_, dim = d[1L:2L], dimnames = dimnames(fits)[1L:2L])
	for (i in 1L:d[1L]) for (j in 1L:d[2L]) {
		res <- par.avg(fits[i, j, ], sefits[i, j, ], weight = w, df = dfs) 
		avgfit[i, j] <- res[1L]
		avgsefit[i, j] <- res[2L]
	}
	colnames(avgfit) <- colnames(avgsefit) <- all.coefnames
	if(fitonly) avgfit else list(fit = avgfit, se.fit = avgsefit)
}

avgsefit <-
function(yall, w, revised.var, full, dfs,  ...) {
	fit <- do.call("cbind", lapply(yall, "[[", "fit"))
	se.fit <- do.call("cbind", lapply(yall, "[[", "se.fit"))
	n <- nrow(fit)
	y <- matrix(0, ncol = 2L, nrow = n)
	hasDfs <- !is.null(dfs) && !anyNA(dfs)
	k <- c(1L, if(hasDfs) 3L else 2L) 
	for(i in 1L:n) y[i, ] <-
		par.avg(fit[i, ], se.fit[i, ], weight = w,
			df = dfs, revised.var = revised.var)[k]
	list(fit = y[, 1L], se.fit = y[, 2L])
}

avgmat <-
function(yall, w, ...) {
	ret <- yall[[1L]]
	ydim <- c(dim(ret), length(yall))
	yarr <- array(unlist(yall), dim = ydim)
	kseq <- seq_len(ydim[2L])
	for(i in seq_len(ydim[1L]))
		for(k in kseq) ret[i, k] <- weighted.mean(yarr[i, k, ], w = w)
	ret
}

std_predict.default <-
function(object, newdata, ...) {
	cl <- match.call()
	if(is.null(newdata)) cl$newdata <- NULL
	cl[[1L]] <- as.name("predict")
	eval.parent(cl)
}

std_predict.multinom <-
std_predict.polr <-
function(object, newdata, type = c("link", "response"), ...) {
	type <- match.arg(type)
	cl <- match.call()
	cl[[1L]] <- as.name("predict")
	cl[['type']] <- "probs"
	rval <- eval.parent(cl)
	if(type == "link") glm.link(object)$linkfun(rval) else rval
}
		
std_predict.merMod <-
function(object, newdata = NULL, type = c("link", "response"),
		 se.fit = FALSE, re.form = NA, ...) {
	if(is.null(newdata)) newdata <- NULL
	predict.merMod(object, newdata, type, se.fit, re.form = re.form, ...)
}

std_predict.lme <-
function(object, newdata = NULL, type = NA, se.fit = FALSE, ...) {
	predict.lme(object, newdata, se.fit = se.fit, ...)
}

std_predict.lm <-
function(object, newdata = NULL, type = c("link", "response", "terms"), se.fit = FALSE, ...) {
	type <- match.arg(type)
	if(!inherits(object, "glm") && type == "link") type <- "response"
	std_predict.default(object, newdata = newdata, type = type, se.fit = se.fit, ...)
}

std_predict.cpglm <-
function(object, newdata = NULL, type = c("link", "response"), ...) {
	type <- match.arg(type)
	as.numeric(predict(object, newdata = newdata, type = type))
	#as.numeric(predict(object, newdata = if(missing(newdata)) NULL else newdata, type = type))
}

invtransform.gls <-
invtransform.lme <-
invtransform.lm <-
function(object, y) y

invtransform.glm <-
invtransform.default <-
function(object, y) {
	link <- glm.link(object)
	if(is.list(y) && all(c("fit", "se.fit") %in% names(y))) {
		y$se.fit <- y$se.fit * abs(link$mu.eta(y$fit))
		y$fit <- link$linkinv(y$fit)
		y
	} else {
		link$linkinv(y)
	}
}

std_predict.zeroinfl <-
function(object, newdata, type = c("link", "response", "prob"), ...) {
	type <- match.arg(type)
	cl <- match.call()
	cl[[1L]] <- as.name("predict")
	if(is.null(cl$newdata)) cl$newdata <- NULL
	if(type == "link") {
		cl$type <- "zero"
		ret <- family(object)$linkfun(eval.parent(cl))
		cl$type <- "count"
		ret <- cbind(ret, log(eval.parent(cl)))
		colnames(ret) <- c("linkphi", "logmu")
		ret
	} else {
		eval.parent(cl)
	}
}

std_predict.hurdle <-
function(object, newdata = NULL, type = c("link", "response", "prob"), ...) {
	type <- match.arg(type)
	cl <- match.call()
	cl[[1L]] <- as.name("predict")
	if(is.null(cl$newdata)) cl$newdata <- NULL
	if(type == "link") {
		cl$type <- "zero"
		ret <- log(eval.parent(cl))
		cl$type <- "count"
		ret <- cbind(ret, log(eval.parent(cl)))
		colnames(ret) <- c("logphi", "logmu")
		ret
	} else {
		eval.parent(cl)
	}
}

std_predict.coxph <-
function (object, newdata, type = NA,
		  predtype = c("lp", "risk", "expected"), ...) {
	predtype <- match.arg(predtype, several.ok = TRUE)
	cl <- match.call()
	cl[[1L]] <- as.name("predict")
	cl$predtype <- NULL
	if(is.null(newdata)) cl$newdata <- NULL
	if(!is.na(type) && type == "terms") return(eval.parent(cl))
	y <- predtypes(cl, predtype)
	if(length(predtype) == 1L) return(y[[1L]])
	if(all(vapply(y, issefit, FALSE))) {
		rval <- vector("list", 2L)
		for(i in 1L:2L) rval[[i]] <- do.call("cbind", lapply(y, "[[", i))
		names(rval) <- c("fit", "se.fit")
	} else {
		rval <- do.call("cbind", y)
		colnames(rval) <- predtype
	}
	return(rval)
}

#{{
# std_predict.coxph <-
# function(object, newdata = NULL, type = NA, se.fit = FALSE, ...) {
	# if(!is.na(type) && type != "response") stop(gettextf("prediction of type '%s' is not implemented", type))
	# cl <- match.call()
	# cl[[1L]] <- as.name("predict")
	# if(missing(newdata)) cl$newdata <- NULL
	# types <- c("lp", "risk", "expected")
	# ntypes <- length(types)
	# y <- vector("list", ntypes)
	# names(y) <- types <- types
	# for(i in seq_len(ntypes)) {
		# cl$type <- types[i]
		# y[[i]] <- eval.parent(cl)
	# }
	# if(se.fit) {
		# rval <- vector("list", ntypes)
		# names(rval) <- nm <- c("fit", "se.fit")
		# for(i in nm) rval[[i]] <- do.call("cbind", lapply(y, "[[", i))
	# } else {
		# rval <- do.call("cbind", y)
	# }
	# rval
	# #rval <- vector("list", 2L)
	# #for(i in 1L:2L) rval[[i]] <- do.call("cbind", lapply(y, "[", , i))
	# #names(rval) <- c("fit", "se.fit")
	# #rval
# }
# std_predict(fm2, lung1[1:4, ], se.fit = T)
# std_predict(fm2, lung1[1:4, ], se.fit = F)
# predict(fm2, type = "risk", se = T)
#}}

std_predict.unmarkedFit <-
function (object, newdata, type = NA, predtype, ...) {
	if(!is.na(type) && type != "response")
		stop(gettextf("prediction of type '%s' is not implemented", type))
	
	predtype <- if(missing(predtype)) names(object@estimates) else
		match.arg(predtype, choices = names(object@estimates), several.ok = TRUE)
	cl <- match.call()
	cl[[1L]] <- as.name("predict")
	if(is.null(newdata)) cl$newdata <- NULL
	y <- predtypes(cl, predtype)
	rval <- vector("list", 2L)
	for(i in 1L:2L) rval[[i]] <- do.call("cbind", lapply(y, "[", , i))
	names(rval) <- c("fit", "se.fit")
	rval
}

std_predict.caic <-
function (object, newdata, type = c("response", "terms"), se.fit = FALSE, interval, ...) {
	cl <- match.call(expand.dots = FALSE)
	if(missing(newdata) || is.null(newdata)) newdata <- caper::caic.table(object)
	do.call("predict", list(object = object$mod, newdata = newdata, type = type, se.fit = se.fit, ...))
}

std_predict.pgls <-
function (object, newdata, type, se.fit, ...) {
	if(missing(newdata)) newdata <- NULL
	predict(object, newdata)[, 1L]
}

std_predict.rq <-
function (object, newdata, type, se.fit, ...) {
	arg <- list(object = object)
	if(!missing(newdata) && !is.null(newdata)) arg$newdata <- newdata
	do.call("predict", arg)
}


invtransform.zeroinfl <-
function(object, y) {
	(1 - object$linkinv(y[, "linkphi"])) * exp(y[, "logmu"])
}

invtransform.hurdle <-
function(object, y) {
	exp(y[, "logphi"] + y[, "logmu"])
}


## merge*(models, ...) functions

mergeContrasts <-
function(models) {
	ctr <- lapply(models, get.contrasts)
	ctrnm <- unique(unlist(lapply(ctr, names)))
	ret <- structure(vector("list", length = length(ctrnm)), names = ctrnm)
	for(x in ctr) {
		for(i in names(x)) {
			if(is.null(ret[[i]])) {
				ret[[i]] <- x[[i]]
			} else 
				if(!identical(ret[[i]], x[[i]]))
					stop(gettextf("inconsistent contrasts in '%s'", i) )
		}
	}
	ret
}

mergeMF <-
function(models) {
	mf <- model.frame(models[[1L]])
	mfNames <- colnames(mf)
	lhs <- asChar(getResponseFormula(mf))

	f <- attr(fixTermsObject(terms(mf)), "term.labels")
	#m <- models[[2]]
	for(m in models[-1L]) {
		mf1 <- model.frame(m)
		if(!identical(lhs, lhs1 <- asChar(getResponseFormula(terms(mf1)))))
			stop("response differs between models: ", sQuote(c(asChar(lhs), lhs1)))
		mf <- cbind(mf, mf1[, !(colnames(mf1) %in% mfNames), drop = FALSE])
		tt1 <- fixTermsObject(terms(mf1))
		f <- c(f, attr(tt1, "term.labels"))
		if(!is.null(attr(tt1,"offset")))
			f <- c(f, sapply(as.list(attr(tt1,"variables")[attr(tt1,"offset") + 1L]), asChar))
	}
	
	f <- reformulate(f[!duplicated(f)])
	f <- as.formula(call("~", parse(text = lhs, keep.source = FALSE)[[1L]], f[[length(f)]]))
	environment(f) <- environment(formula(models[[1L]]))
	tt <- fixTermsObject(terms(f))
	mf <- mf[, rownames(attr(tt, "factors"))]
	attr(tt, "dataClasses") <- vapply(mf, .MFclass, "")
	attr(mf, "terms") <- tt
	mf
}

offsetTermNames <- function(x)
vapply(as.list(attr(x, "variables")[attr(x,"offset") + 1L]), deparse, "", control = NULL)

offsetWeights <-
function(wts, Terms, models) {
	if(is.null(off <- attr(Terms, "offset")))
		return(NULL)
	offnames <- rownames(attr(Terms, "factors"))[off]
	n <- length(offnames)
	v <-  matrix(vapply(models, function(x) {
		offnames %in% offsetTermNames(terms(x))
		}, logical(n)), nrow = n)
	offwts <- vapply(1L:n, function(i) sum(wts[v[i, ]]), 0)
	names(offwts) <- offnames
	offwts
}

## orders terms alphabetically in 'terms' object
fixTermsObject <-
function(x, peel = TRUE) {
	
	peelfun <- function(z) if(is.call(z))
		paste0(as.character(z[-1L]), collapse = " ") else as.character(z)
	
	factors <- attr(x, "factors")
    if (length(factors) != 0L) {
		
		z <- rep(1L, nrow(factors))
		z[attr(x, "response")] <- 0L
		if(hasOff <- !is.null(attr(x, "offset")))
			z[attr(x, "offset")] <- 2L
		charvar <- if (peel) sapply(as.list(attr(x, "variables")[-1L]), peelfun) else
					rownames(factors)
		ov <- order(z, charvar)
		factors <- factors[ov, , drop = FALSE]
		charvar <- charvar[ov]
		v <- rownames(factors)
		
		lab <- lab_ord <- character(ncol(factors))
		lfac <- factors != 0L
		for(i in 1L:ncol(factors)) {
			j <- lfac[, i]
			lab[i] <- paste0(v[j], collapse = ":")
			lab_ord[i] <- paste0(charvar[j], collapse = ":")
		}
		o <- order(attr(x, "order"), lab_ord)
		lab <- lab[o]
		ans <- reformulate(c(lab, offsetTermNames(x)), intercept = attr(x, "intercept"))
		
		if(attr(x,"response") != 0) ans <-
			as.formula(call("~", attr(x, "variables")[[attr(x, "response") + 1L]], ans[[2L]]))
		attributes(ans) <- attributes(x)
		attr(ans, "factors") <- factors[, o, drop = FALSE]
		colnames(attr(ans,"factors")) <- attr(ans, "term.labels") <- lab
		
		if(hasOff) attr(ans, "offset") <- which(z[ov] == 2L)
		for(j in c("variables", "predvars"))
			if(!is.null(attr(ans, j))) attr(ans, j) <- attr(ans, j)[c(1L, ov + 1L)]
		if(!is.null(attr(ans, "dataClasses")))
			attr(ans, "dataClasses") <- attr(ans, "dataClasses")[ov]
		ans		
	} else x
}

get.contrasts <- function(x) UseMethod("get.contrasts")
get.contrasts.lm <- function(x) x$contrasts
get.contrasts.averaging <- function(x) {
	mergeContrasts(getModelList(x))
}


getModelList <-
function(object, error = TRUE) {
	if(is.null(models <- attr(object, "modelList")))
		if(error) stop("component models not included in this 'averaging' object")
	invisible(models)
}


lincomb <-
function(beta, mf, newdata, contr, offwts) {
	if(is.null(newdata)) {
		## TODO
		.NotYetImplemented()
	} else {
		xlev <- .getXlevels(Terms <- terms(mf), mf)
		m <- model.frame(Terms, newdata, na.action = na.fail, xlev = xlev)
		X <- model.matrix(Terms, m, contrasts.arg = contr)
		if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
		if (!is.null(off.num <- attr(Terms, "offset"))) {
			v <- attr(Terms, "variables")
			j <- match(names(offwts), rownames(attr(Terms, "factors")))
			for(i in seq_along(v))
				v[[j[i] + 1L]][[2L]] <- call("*", v[[j[i] + 1L]][[2L]], offwts[i])
			offset <- rep(0, nrow(X))
			for (i in off.num) offset <- offset + eval(v[[i + 1L]], newdata)
		}
	}
	cn <- fixCoefNames(colnames(X), FALSE)
	ans <- (X %*% beta[cn])[, 1L]
	if (!is.null(off.num)) ans <- ans + offset
	ans
}

