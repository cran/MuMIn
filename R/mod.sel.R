`mod.sel` <-
function(object, ..., rank = AICc, rank.args = NULL) {
	if (missing(object) && length(models) > 0) {
		object <- models[[1L]]
	} else if (inherits(object, "list")) {
		if(length(object) ==  0) stop("At least one model must be given")
		models <- object
		object <- models[[1L]]
	} else {
		models <- list(object, ...)
		names(models)[1] <- deparse(substitute(object))
	}
	if(length(models) == 0) stop("At least one model must be given")

	if(is.null(names(models)) || any(is.na(names(models))))
		names(models) <- seq_along(models)
	names(models) <- make.unique(names(models), sep="")

	rank <- .getRank(rank, rank.args = rank.args, object = object)

	ICname <- deparse(attr(rank,"call")[[1]])
	
	

	all.terms <- c("(Intercept)", unique(unlist(lapply(models, getAllTerms))))
	#browser()
	all.coef <- fixCoefNames(unique(unlist(lapply(lapply(models, coeffs), names))))
	
	logLik <- .getLogLik()
	
	
	j <- !(all.terms %in% all.coef)
	d <- as.data.frame(t(sapply(models, matchCoef, all.terms=all.terms)))
	d[,j] <- lapply(d[,j, drop=FALSE], function(x) factor(is.nan(x), levels=c(F, T), labels=c("", "+")))
	
	ret <- as.data.frame(t(vapply(models, function(x) { 
		ll <- logLik(x)
		c(attr(ll, "df"), attr(ll, "nobs"), ll, rank(x)) 
		}, structure(double(4), names=c("df", "nobs", "logLik", ICname)))))
	
	ret <- cbind(d, ret)
	ret[,"delta"] <- ret[, ICname] - min(ret[, ICname])
	ret[,"weight"] <- Weights(ret[,ICname])
	
	ret <- ret[order(ret[,"delta"], decreasing = FALSE), ]

	attr(ret, "terms") <- all.terms
	rownames(ret) <- names(models)
	class(ret) <- c("model.selection", "data.frame")
	ret
}	
	
	# data(Cement)
	# Cement$X1 <- cut(Cement$X1, 3)
	# Cement$X2 <- cut(Cement$X2, 2)
	# fm1 <- glm(formula = y ~ X + X1 + X2 * X3, data = Cement)
	# dd <- dredge(fm1)
	# mod <- get.models(dd, 1:20)
	
	# names(coeffs(fm1))

# }
