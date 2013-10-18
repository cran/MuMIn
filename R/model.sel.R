#TODO: checking if models are fitted to the same dataset <- model.avg

`model.sel` <-
`mod.sel` <-
function (object, ...) UseMethod("model.sel")

`model.sel.model.selection` <-
function (object, rank = NULL, rank.args = NULL, ...) {
	if(!is.null(rank)) {
		rank <- .getRank(rank, rank.args = rank.args)
		ic <- tryCatch(sapply(logLik(object), rank), error = function(e) e)
		if(!inherits(ic, "error") && is.numeric(ic)) {
			oldRankCol <- as.character(attr(attr(object, "rank"), "call")[[1L]])
			rankCol <- as.character(attr(rank, "call")[[1L]])
			colnames(object)[colnames(object) == oldRankCol] <- rankCol
			object[, rankCol] <- ic
			object$delta <- ic - min(ic)
			object$weight <- Weights(ic)
			ret <- object[order(ic), ]
			#attr(ret, "order") <- o
			attr(ret, "rank") <- rank
			attr(ret, "rank.call") <- attr(rank, "call")
		} else {
			message("'rank' cannot be applied to 'logLik' object. Recreating model fits.")
			models <- get.models(object, seq.int(nrow(object)))
			ret <- model.sel.default(models, rank = rank)
		}
		return(ret)
	} else return(object)
}


`model.sel.default` <-
function(object, ..., rank = NULL, rank.args = NULL) {
	.makemnames <- function(cl) {
		cl[1L] <- cl$rank <- cl$rank.args <- NULL
		unlist(.makeListNames(cl))
	}

	if (missing(object) && length(models <- list(...)) > 0L) {
		object <- models[[1L]]
		names(models) <- .makemnames(sys.call())
	} else if (inherits(object, "list")) {
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

	.checkModels(models, FALSE)

	if(is.null(names(models)) || any(is.na(names(models))))
		names(models) <- seq_along(models)
	names(models) <- make.unique(names(models), sep = "")

	rank <- .getRank(rank, rank.args = rank.args, object = object)
	ICname <- deparse(attr(rank, "call")[[1L]])
	allTermsList <- lapply(models, getAllTerms, intercept = TRUE)
	random.terms <- lapply(allTermsList, attr, "random.terms")
	all.terms <- unique(unlist(allTermsList, use.names = FALSE))
	all.coef <- fixCoefNames(unique(unlist(lapply(lapply(models, coeffs), names),
		use.names = FALSE)))

	## TODO: case when models belong to different classes using logLik or qLik 
	## - give error
	LL <- .getLik(models[[1L]])
	logLik <- LL$logLik
	lLName <- LL$name

	j <- !(all.terms %in% all.coef)
	#d <- as.data.frame(t(sapply(models, matchCoef, all.terms = all.terms)))

	mcoeflist <- lapply(models, matchCoef, all.terms = all.terms, allCoef = TRUE)
	d <- as.data.frame(do.call("rbind", mcoeflist))

	retCoefTable <-	lapply(mcoeflist, attr, "coefTable")

	d[,j] <- lapply(d[,j, drop = FALSE], function(x) factor(is.nan(x),
		levels = TRUE, labels = "+"))

	ret <- vapply(models, function(x) {
		ll <- logLik(x)
		ic <- tryCatch(rank(x), error = function(e) e)
		if(inherits(ic, "error")) {
			ic$call <- sys.call(sys.nframe() - 4L)
			ic$message <- gettextf("evaluating 'rank' for an object failed with message: %s", ic$message)
			stop(ic)
		}
		c(attr(ll, "df"), ll, ic)
		}, structure(double(3L), names=c("df", lLName, ICname)))
	ret <- as.data.frame(t(ret))

	ret <- cbind(d, ret)
	ret[, "delta"] <- ret[, ICname] - min(ret[, ICname])
	ret[, "weight"] <- Weights(ret[,ICname])
	o <- order(ret[, "delta"], decreasing = FALSE)

	descrf <- modelDescr(models)
	descrf$model <- NULL
	if(nlevels(descrf$family) == 1L) descrf$family <- NULL
	if(ncol(descrf)) {
		i <- seq_len(length(all.terms))
		ret <- cbind(ret[, i], descrf, ret[, -i])
	}

	rownames(ret) <- names(models)

	ret <- structure(
		ret[o, ],
		terms = structure(all.terms, interceptLabel =
			unique(unlist(lapply(allTermsList, attr, "interceptLabel")))),
		calls = lapply(models, .getCall)[o],
		order = o,
		rank = rank,
		rank.call = attr(rank, "call"),
		call = match.call(expand.dots = TRUE),
		nobs = nobs(models[[1L]]),
		coefTables = retCoefTable[o],
		vCols = colnames(descrf),
		class = c("model.selection", "data.frame")
	)

	if (!all(sapply(random.terms, is.null)))
		attr(ret, "random.terms") <- random.terms[o]

	ret
}
