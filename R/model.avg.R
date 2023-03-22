`model.avg` <-
function (object, ..., revised.var = TRUE) {
	if (isTRUE("method" %in% names(match.call())))
		stop("argument 'method' is defunct")
	UseMethod("model.avg")
}

.coefarr.avg <-
function(cfarr, weight, revised.var, full, alpha) {	
	weight <- weight / sum(weight)
	nCoef <- dim(cfarr)[3L]
	if(full) {
		nas <- is.na(cfarr[, 1L, ]) & is.na(cfarr[, 2L, ])
		cfarr[, 1L, ][nas] <- cfarr[, 2L, ][nas] <- 0
		#cfarr[, 1L:2L, ][is.na(cfarr[, 1L:2L, ])] <- 0
		if(!all(is.na(cfarr[, 3L, ])))
			cfarr[ ,3L, ][is.na(cfarr[ , 3L, ])] <- Inf
	}
	
	avgcoef <- array(dim = c(nCoef, 5L),
		dimnames = list(dimnames(cfarr)[[3L]], c("Estimate",
		"Std. Error", "Adjusted SE", "Lower CI", "Upper CI")))
	for(i in seq_len(nCoef))
		avgcoef[i, ] <- par.avg(cfarr[, 1L, i], cfarr[, 2L, i], weight,
			df = cfarr[, 3L, i], alpha = alpha, revised.var = revised.var)
		
	avgcoef[is.nan(avgcoef)] <- NA
	return(avgcoef)
}

`model.avg.model.selection` <-
function(object, subset, fit = FALSE, ..., revised.var = TRUE) {

	if(!missing(subset)) {
		cl <- match.call()
		cl[[1L]] <- as.name("subset")
		names(cl)[2L] <- "x"
		object <- eval.parent(cl[1L:3L])
	}
	# TODO: unify refitting conditions in model.avg and model.sel
	
	if(fit || !missing(...)) {
		cl <- match.call()
		cl$fit <- NULL
		arg1 <- names(cl)[-(1L:2L)] %in% names(formals("model.avg.default"))
		cl1 <- cl[c(TRUE, TRUE, !arg1)]
		cl1[[1L]] <- as.name("get.models")
		if(is.null(cl1[["subset"]])) cl1[["subset"]] <- NA
		# TODO: subset = TRUE

		cl2 <- cl[c(TRUE, TRUE, arg1)]
		cl2[[2L]] <- cl1
		cl2[[1L]] <- as.name("model.avg")
		#message("Re-fitting model objects...")
		return(eval(cl2, parent.frame()))
	}

	if(nrow(object) <= 1L) stop("'object' consists of only one model")
	
	ct <- attr(object, "coefTables")
	cfarr <- coefArray(ct)
	weight <- Weights(object)

	cfmat <- as.matrix(cfarr[, 1L, ])
	cfmat[is.na(cfmat)]<- 0
	coefMat <- array(dim = c(2L, ncol(cfmat)),
		dimnames = list(c("full", "subset"), dimnames(cfarr)[[3L]]))
	
	coefMat[1L, ] <- drop(weight %*% cfmat)
	coefMat[2L, ] <- coefMat[1L, ] / colSums(array(weight *
		as.numeric(!is.na(cfarr[, 1L, ])), dim = dim(cfmat)))
	coefMat[is.nan(coefMat)] <- NA_real_

	#allterms1 <- lapply(attr(object, "calls"), function(x)
		#getAllTerms(as.formula(x[[switch(as.character(x[[1L]]),
			#lme=, lme.formula= "fixed", gls= "model", "formula")]])))
	all.terms <- attr(object, "terms") # TERMS
	all.vterms <- all.terms[!(all.terms %in% attr(all.terms, "interceptLabel")
		| apply(is.na(object[, all.terms, drop = FALSE]), 2L, all))]
	#allterms1 <- apply(!is.na(object[, all.vterms, drop = FALSE]), 1L, function(x) all.vterms[x])
	allterms1 <- applyrns(!is.na(object[, all.vterms, drop = FALSE]), function(x) all.vterms[x])
	allmodelnames <- .modelNames(allTerms = allterms1, uqTerms = all.vterms)

	mstab <- itemByType(object, c("df", "loglik", "ic", "delta", "weight"))
	rownames(mstab) <- allmodelnames
	
	.Debug(.Generic <- "model.avg")

	ret <- list(
		msTable = structure(as.data.frame(mstab, stringsAsFactors = TRUE),
			term.codes = attr(allmodelnames, "variables")),
		coefficients = coefMat,
		coefArray = cfarr,
		sw = sw(object),
		x = NULL,
		residuals = NULL,
		formula = if(!is.null(attr(object, "global")))
			formula(attr(object, "global")) else NULL,
		call = {
			cl <- match.call()
			cl[[1L]] <- as.name(.Generic)
			cl
		}
	)

	attr(ret, "rank") <- attr(object, "rank")
	if(is.null(attr(object, "modelList"))) {
		attr(ret, "model.calls") <- attr(object, "model.calls")
		attr(ret, "interceptLabel") <- attr(attr(object, "terms"), "interceptLabel")
	} else {
		attr(ret, "modelList") <- attr(object, "modelList")
	}
	attr(ret, "beta") <- attr(object, "beta")
	attr(ret, "nobs") <- attr(object, "nobs")
	attr(ret, "revised.var") <- revised.var
	class(ret) <- "averaging"
	return(ret)
}

`model.avg.default` <-
function(object, ..., beta = c("none", "sd", "partial.sd"),
		 rank = NULL, rank.args = NULL, revised.var = TRUE,
		 dispersion = NULL, ct.args = NULL) {

	if (is.object(object)) {
		models <- list(object, ...)
        rank <- .getRank(rank, rank.args = rank.args, object = object) 
	} else {
		if(length(object) == 0L) stop("'object' is an empty list")
		models <- object
		object <- object[[1L]]
        if (!is.null(rank) || is.null(rank <- attr(models, "rank"))) {
            rank <- .getRank(rank, rank.args = rank.args, object = object)
      	}
	}
	
	strbeta <- betaMode <- NULL
	eval(.expr_beta_arg)

	nModels <- length(models)
	if(nModels == 1L) stop("only one model supplied. Nothing to do")
	checkIsModelDataIdentical(models)
	
	testSmoothKConsistency(models) # for gam, if any

	ICname <- asChar(attr(rank, "call")[[1L]])

	allterms1 <- lapply(models, getAllTerms)
	all.terms <- unique(unlist(allterms1, use.names = FALSE))

	# sort by level (main effects first)
	all.terms <- all.terms[order(vapply(gregexpr(":", all.terms),
		function(x) if(x[1L] == -1L) 0L else length(x), 1L), all.terms)]

	# allmodelnames <- modelNames(models, asNumeric = FALSE,
		# withRandomTerms = FALSE, withFamily = FALSE)
	allmodelnames <- .modelNames(allTerms = allterms1, uqTerms = all.terms)

	#if(is.null(names(models))) names(models) <- allmodelnames
	
	coefTableCall <- if(betaMode == 2L) 
		 call("std.coef", as.symbol("m"), partial.sd = TRUE)
		else call("coefTable", as.symbol("m"))
	if(!is.null(dispersion))
		coefTableCall[['dispersion']] <- as.symbol("d")
	for(a in names(ct.args))
		coefTableCall[[a]] <- ct.args[[a]]
		

	.DebugPrint(coefTableCall)
    
    # NOTE: first argument in coefTableCall is "m" and "d" for dispersion
	coefTables <-
	    mapply(function(m, d) {
	        rval <- eval(coefTableCall)
	        rownames(rval) <- fixCoefNames(rownames(rval))
	        rval
	    }, m = models, d = if(is.null(dispersion)) NA else dispersion,
			SIMPLIFY = FALSE)
	
    
	#coefTables <- vector(nModels, mode = "list")
    # # NOTE: first argument in coefTableCall is "models[[i]]"
	# for(i in seq_len(nModels)) {
		# coefTables[[i]] <- eval(coefTableCall)
		# rownames(coefTables[[i]]) <- fixCoefNames(rownames(coefTables[[i]]))
	# }
	
	
	# check if models are unique:
	mcoeffs <- lapply(coefTables, "[", , 1L)
	dup <- unique(sapply(mcoeffs, function(i) which(sapply(mcoeffs, identical, i))))
	dup <- dup[sapply(dup, length) > 1L]
	if (length(dup) > 0L) stop("models are not unique. Duplicates: ",
		prettyEnumStr(sapply(dup, paste0, collapse = " = "),
			quote = "'"))

	LL <- .getLik(object)
	logLik <- LL$logLik
	lLName <- LL$name

	ic <- vapply(models, rank, 0)
	logLiks <- lapply(models, logLik)
	delta <- ic - min(ic)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))
	model.order <- order(weight, decreasing = TRUE)

	# ----!!! From now on, everything MUST BE ORDERED by 'weight' !!!-----------
	mstab <- cbind(df = vapply(logLiks, attr, 0, "df"),
		logLik = as.numeric(logLiks), IC = ic, delta = delta, weight = weight,
		deparse.level = 0L)
	if(!is.null(dispersion)) mstab <- cbind(mstab, Dispersion = dispersion)
	rownames(mstab) <- allmodelnames
	mstab <- mstab[model.order, ]
	weight <- mstab[, "weight"] # has been sorted in table
	models <- models[model.order]
	coefTables <- coefTables[model.order]

	if (betaMode == 1L) {
		response.sd <- sd(model.response(model.frame(object)))
		coefTables <-
			mapply(function(m, ct) {
				X <- model.matrix(m)
				ct[, 1L:2L] <- ct[, 1L:2L] *
					apply(X[, match(rownames(ct), colnames(X)), drop = FALSE],
						2L, sd) / response.sd
			}, models, coefTables, SIMPLIFY = FALSE)
	}

	cfarr <- coefArray(coefTables)
	cfmat <- array(cfarr[, 1L, ], dim = dim(cfarr)[-2L], dimnames = dimnames(cfarr)[-2L])
	cfmat[is.na(cfmat)]<- 0
	coefMat <- array(NA_real_, dim = c(2L, ncol(cfmat)),
		dimnames = list(c("full", "subset"), colnames(cfmat)))
	coefMat[1L, ] <- drop(weight %*% cfmat)
	coefMat[2L, ] <- coefMat[1L, ] / colSums(array(weight *
		as.numeric(!is.na(cfarr[, 1L, ])), dim = dim(cfmat)))
	coefMat[is.nan(coefMat)] <- NA_real_
		
    names(all.terms) <- seq_along(all.terms)
	colnames(mstab)[3L] <- ICname

	# Benchmark: 3.7x faster
	#system.time(for(i in 1:10000) t(array(unlist(p), dim=c(length(all.terms),length(models)))))
	#system.time(for(i in 1:10000) do.call("rbind", p))
	
	vpresent <- do.call("rbind", lapply(models, function(x)
		all.terms %in% getAllTerms(x)))
	
	if(all(dim(vpresent) > 0L)) {
		sw <- apply(weight * vpresent, 2L, sum)
		names(sw) <- all.terms
		o <- order(sw, decreasing = TRUE)
		sw <- sw[o]
		attr(sw, "n.models") <- structure(colSums(vpresent)[o], names = all.terms)
		class(sw) <- c("sw", "numeric")
	} else {
		sw <- structure(integer(0L), n.models = integer(0L), class = c("sw", "numeric"))
	}
	
	mmxs <- tryCatch(cbindDataFrameList(lapply(models, model.matrix)),
					 error = return_null, warning = return_null)

	# Far less efficient:
	#mmxs <- lapply(models, model.matrix)
	#mx <- mmxs[[1]];
	#for (i in mmxs[-1])
	#	mx <- cbind(mx, i[,!(colnames(i) %in% colnames(mx)), drop=FALSE])

	# residuals averaged (with brute force)
	#rsd <- tryCatch(apply(vapply(models, residuals, residuals(object)), 1L,
		#weighted.mean, w = weight), error = return_null)
	#rsd <- NULL
	## XXX: how to calc residuals ?
	
	modelClasses <- lapply(models, class)
	frm <-
	if(all(vapply(modelClasses[-1L], identical, FALSE, modelClasses[[1L]]))) {
		trm <- tryCatch(terms(models[[1L]]),
				error = function(e) terms(formula(models[[1L]])))
		response <- attr(trm, "response")
		m1 <- models[[1L]]
		makeArgs(m1, all.terms, opt = list(
			response = if(response > 0L) attr(trm, "variables")[[response + 1L]] else NULL,
			gmFormulaEnv = environment(formula(m1)),
			intercept = ! identical(unique(unlist(lapply(allterms1, attr, "intercept"))), 0),
			interceptLabel = unique(unlist(lapply(allterms1, attr, "interceptLabel"))),
			#	random = attr(allTerms0, "random"),
			gmCall = get_call(m1),
			gmEnv = parent.frame(),
			allTerms = all.terms,
			random = . ~ .
			))[[1L]]
	} else NA
	
	.Debug(.Generic <- "model.avg")
	
	ret <- list(
		msTable = structure(as.data.frame(mstab, stringsAsFactors = TRUE),
			term.codes = attr(allmodelnames, "variables")),
		coefficients = coefMat,
		coefArray = cfarr,
		sw = sw,
		x = mmxs,
		residuals = NULL, # no residuals, as they can be calculated in several ways
		formula = frm,
		call = {
			cl <- match.call()
			cl[[1L]] <- as.name(.Generic)
			cl
		}
	)

	attr(ret, "rank") <- rank
	attr(ret, "modelList") <- models
	attr(ret, "beta") <- strbeta
	attr(ret, "nobs") <- nobs(object)
	attr(ret, "revised.var") <- revised.var
	class(ret) <- "averaging"
	return(ret)
}

.checkFull <-
function(object, full, warn = TRUE) {
	if(isTRUE(attr(object, "arm")) && !full) {
		if(warn) cry(-1L, "'subset' averaged coefficients are not available with ARM algorithm",
					 warn = TRUE)
		return(TRUE)
	} else return(full)
}

## XXX: backward compatibility (< 0.15.0):
upgrade_averaging_object <-
function(x) {
	if(is.matrix(x$coefficients)) return(x)
	if(all(c("coefTable", "coef.shrinkage") %in% names(x))) {
		x$coefficients <- rbind(full = x$coef.shrinkage, subset = x$coefTable[, 1L])
		x$coefTable <- NULL
	} else stop("'averaging' object is corrupt")
	x
}
