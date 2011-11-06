#TODO: checking if models are fitted to the same dataset <- model.avg

`model.sel` <-
`mod.sel` <-
function (object, ...) UseMethod("mod.sel")


`mod.sel.model.selection` <-
function (object, rank = NULL, rank.args = NULL, ...) {
	#if(!is.null(rank)) .NotYetUsed("rank")
	if(!is.null(rank)) {
		models <- get.models(object, seq.int(nrow(object)))
		ret <- mod.sel.default(models, rank = .getRank(rank,
			rank.args = rank.args, object = models[[1L]]))
		return(ret)
	} else {
		return(object)
	}
}


`mod.sel.default` <-
function(object, ..., rank = NULL, rank.args = NULL) {

	if (missing(object) && length(models <- list(...)) > 0L) {
		object <- models[[1L]]
	} else if (inherits(object, "list")) {
		if(length(object) ==  0L) stop("at least one model must be given")
		models <- object
		object <- models[[1L]]
	} else {
		models <- list(object, ...)
		names(models)[1L] <- deparse(substitute(object))
	}
	if(length(models) == 0L) stop("at least one model must be given")

	.checkModels(models, FALSE)

	if(is.null(names(models)) || any(is.na(names(models))))
		names(models) <- seq_along(models)
	names(models) <- make.unique(names(models), sep = "")

	rank <- .getRank(rank, rank.args = rank.args, object = object)
	ICname <- deparse(attr(rank, "call")[[1L]])
	all.terms <- unique(unlist(lapply(models, getAllTerms, intercept = TRUE)))
	all.coef <- fixCoefNames(unique(unlist(lapply(lapply(models, coeffs), names))))

	logLik <- .getLogLik()

	j <- !(all.terms %in% all.coef)
	d <- as.data.frame(t(sapply(models, matchCoef, all.terms=all.terms)))
	d[,j] <- lapply(d[,j, drop=FALSE], function(x) factor(is.nan(x),
		levels=c(F, T), labels=c("", "+")))

	ret <- as.data.frame(t(vapply(models, function(x) {
		ll <- logLik(x)
		c(attr(ll, "df"), ll, rank(x))
		}, structure(double(3L), names=c("df", "logLik", ICname)))))

	ret <- cbind(d, ret)
	ret[, "delta"] <- ret[, ICname] - min(ret[, ICname])
	ret[, "weight"] <- Weights(ret[,ICname])
	o <- order(ret[, "delta"], decreasing = FALSE)

	rownames(ret) <- names(models)
	ret <- ret[o, ]

	attr(ret, "terms") <- all.terms
	attr(ret, "calls") <- lapply(models, .getCall)[o]
	attr(ret, "order") <- o
	attr(ret, "rank") <- rank
	attr(ret, "rank.call") <- attr(rank, "call")
	attr(ret, "call") <- match.call(expand.dots = TRUE)

	class(ret) <- c("model.selection", "data.frame")
	ret
}
