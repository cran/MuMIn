#TODO: checking if models are fitted to the same dataset <- model.avg

`model.sel` <-
function (object, ...) UseMethod("model.sel")

`model.sel.averaging` <-
function (object, rank = attr(object, "rank"),
    rank.args = NULL, ...,  beta = c("none", "sd", "partial.sd"), extra) {
    model.sel(get.models(object, NA), rank = rank, rank.args = rank.args, ...,
       beta = beta, extra)
}

`model.sel.model.selection` <-
function (object, rank = NULL, rank.args = NULL, fit = NA, ...,
		  beta = c("none", "sd", "partial.sd"),
		  extra) {

	strbeta <- betaMode <- NULL
	eval(.expr_beta_arg)
	
	reFit <- !missing(extra) || (strbeta != attr(object, "beta"))
	if(!is.null(rank.args) && !identical(fit, FALSE)) reFit <- TRUE
	
	if(!isTRUE(fit) && !is.null(rank)) {
		rank <- .getRank(rank, rank.args = rank.args)
		ic <- tryCatch(sapply(logLik(object), rank), error = identity)
		if(inherits(ic, "error") || !is.numeric(ic)) {
			#message("'rank' cannot be applied to 'logLik' object. Re-fitting model objects.")
			reFit <- TRUE
		}
	} #else rank <- .getRank(attr(object, "rank"))

    # XXX: sometimes the message is shown but no refitting is needed or done.
	if(reFit && !isTRUE(fit)) {
		if(anyNA(fit)) message("Re-fitting models...")
			else stop("cannot proceed without re-fitting models ('fit' is FALSE)")
	}
	
	if(isTRUE(fit) || reFit) {
		#message("to compute 'extras' or beta-weights, need to re-fit model objects.")
		cl <- match.call()
		ss <- if(is.null(cl$subset)) TRUE else cl$subset
		models <- do.call("get.models", list(object, subset = ss), envir = parent.frame())
		cl$subset <- NULL
		cl$object <- models
		rval <- do.call("model.sel", as.list(cl), envir = parent.frame())
	} else if(!is.null(rank)) {
		newRankName <- as.character(attr(rank, "call")[[1L]])
		message(gettextf("New rank '%s' applied to logLik objects", newRankName))
		k <- type2col(object, "ic")
		attr(object, "names")[k] <- names(attr(object, "column.types"))[k] <-
			newRankName
		itemByType(object, "ic") <- ic
		itemByType(object, "delta") <- ic - min(ic, na.rm = TRUE)
		itemByType(object, "weight") <- Weights(ic)
		rval <- object[order(ic), ]
		attr(rval, "rank") <- rank
	} else rval <- object
	return(rval)
}

`model.sel.default` <-
function(object, ..., rank = NULL, rank.args = NULL,
		 beta = c("none", "sd", "partial.sd"),
		 extra) {
	.makemnames <- function(cl) {
		cl[c("rank", "rank.args", "beta", "extra")] <- NULL
		unlist(.makeListNames(cl[-1L]))
	}
	
	strbeta <- betaMode <- NULL
	eval(.expr_beta_arg)

	if (missing(object) && length(models <- list(...)) > 0L) {
		object <- models[[1L]]
		names(models) <- .makemnames(sys.call())
	} else if (is.list(object) && !is.object(object)) {
		if(length(object) ==  0L) stop("at least one model must be given")
		models <- object
		object <- models[[1L]]
		names(models) <- unlist(.makeListNames(models))
	} else {
		models <- list(object, ...)
		if(length(models) > 1L) {
			names(models) <- .makemnames(sys.call())
		} else {
			names(models)[1L] <- unlist(.makeListNames(list(substitute(object))))
		}
	}

	if(length(models) == 0L) stop("at least one model must be given")

	checkIsModelDataIdentical(models, FALSE)

	if(is.null(names(models)) || anyNA(names(models)))
		names(models) <- seq_along(models)
	names(models) <- make.unique(names(models), sep = "")

    if(is.null(rank) && !missing(rank.args))
        warning("'rank.args' ignored with no 'rank' given")
        
	rank <- .getRank(rank, rank.args = rank.args, object = object)
	ICname <- asChar(attr(rank, "call")[[1L]])
	allTermsList <- lapply(models, getAllTerms, intercept = TRUE)
	random.terms <- lapply(allTermsList, attr, "random.terms")
	all.terms <- unique(unlist(allTermsList, use.names = FALSE))
    
    lapply(models, function(fit) {
        if(any(dup <- duplicated(cfn <- names(coeffs(fit)))))
        cry(-2L, "models cannot have duplicated coefficient names: %s",
             prettyEnumStr(cfn[dup]))
    })

    
	all.coef <- fixCoefNames(unique(unlist(lapply(lapply(models, coeffs), names),
		use.names = FALSE)))

	## TODO: case when models belong to different classes using logLik or qLik 
	## - give error
	LL <- .getLik(models[[1L]])
	logLik <- LL$logLik
	lLName <- LL$name

	j <- !(all.terms %in% all.coef)
	#d <- as.data.frame(t(sapply(models, matchCoef, all.terms = all.terms)))

	coefTables <- lapply(models, matchCoef, all.terms = all.terms,
						allCoef = TRUE, beta = betaMode)
	d <- as.data.frame(do.call("rbind", coefTables))
	coefTables <- lapply(coefTables, attr, "coefTable")

	d[,j] <- lapply(d[,j, drop = FALSE], function(x) factor(is.nan(x),
		levels = TRUE, labels = "+"))
	
	rval <- vapply(models, function(x) {
		ll <- logLik(x)
		ic <- tryCatch(rank(x), error = function(e) e)
		if(inherits(ic, "error")) {
			ic$call <- sys.call(sys.nframe() - 4L)
			ic$message <- gettextf("evaluating 'rank' failed with message: %s",
				ic$message)
			stop(ic)
		}
		c(attr(ll, "df"), ll, ic)
		}, structure(double(3L), names = c("df", lLName, ICname)))
	rval <- as.data.frame(t(rval), stringsAsFactors = TRUE)
	rval <- cbind(d, rval)
    o <- order(rval[, ICname], decreasing = FALSE)
	rval[, "delta"] <- rval[, ICname] - rval[o[1L], ICname]
	rval[, "weight"] <- Weights(rval[, ICname])
	mode(rval[, "df"]) <- "integer"
	
	
	descrf <- modelDescr(models)
	descrf$model <- NULL
	if(nlevels(descrf$family) == 1L) descrf$family <- NULL
	if(ncol(descrf)) {
		i <- seq_len(length(all.terms))
		rval <- cbind(rval[, i, drop = FALSE], descrf, rval[, -i, drop = FALSE],
			deparse.level = 0L)
	}
	
	if(!missing(extra) && length(extra) != 0L) {
		# a cumbersome way of evaluating a non-exported function in a parent frame:
		#extra <- eval.parent(call(".get.extras", substitute(extra)))
		extra <- eval.parent(as.call(list(call("get", ".get.extras",
			envir = call("asNamespace",
			.packageName), inherits = FALSE), substitute(extra),
			r2nullfit = TRUE)))
		res <- lapply(models, function(x) unlist(lapply(extra, function(f) f(x))))
		extraResultNames <- unique(unlist(lapply(res, names)))
		nextra <- length(extraResultNames)
		i <- seq_len(length(all.terms))
		rval <- cbind(rval[, i, drop = FALSE], do.call("rbind", lapply(res, function(x) {
			if(length(x) < nextra) {
				tmp <- rep(NA_real_, nextra)
				tmp[match(names(x), extraResultNames)] <- x
				tmp
			} else x
		})), rval[, -i, drop = FALSE])
	} else nextra <- 0L
	row.names(rval) <- names(models)
	
	rval <- structure(
		rval[o, , drop = FALSE],
		 # TERMS
		terms = structure(all.terms, interceptLabel =
			unique(unlist(lapply(allTermsList, attr, "interceptLabel")))),
		model.calls = lapply(models, get_call)[o],
		model.family = lapply(models, function(x) tryCatch(family(x), error = function(e) NULL)),
		modelList = models[o],
		order = o,
		rank = rank,
		beta = strbeta,
		call = match.call(),
		nobs = nobs(models[[1L]]),
		coefTables = coefTables[o],
		vCols = colnames(descrf),
		column.types = {
			colTypes <- c(terms = length(all.terms), varying = ncol(descrf), 
				extra = nextra, df = 1L, loglik = 1L, ic = 1L, delta = 1L,
				weight = 1L)
			column.types <- rep(1L:length(colTypes), colTypes)
			names(column.types) <- colnames(rval)
			lv <- 1L:length(colTypes)
			factor(column.types, levels = lv, labels = names(colTypes)[lv])
		},
		class = c("model.selection", "data.frame")
	)
	if(!("class" %in% colnames(rval)))
		attr(rval, "model.class") <- class(models[[1L]])[1L]
	
	if (!all(sapply(random.terms, is.null)))
		attr(rval, "random.terms") <- random.terms[o]

	rval
}
