
`predict.averaging` <-
function(object, newdata = NULL, se.fit = FALSE, interval = NULL,
	type = NA, backtransform = FALSE,
	full = TRUE, ...) {
	## XXX: backward compatibility:
	    
    object <- upgrade_averaging_object(object)
	full <- .checkFull(object, full) 

	if (!missing(interval)) .NotYetUsed("interval", error = FALSE)
	
	if(backtransform && !is.na(type) && type == "response")
		warning("back-transforming predictions already on response scale")

	models <- attr(object, "modelList")
	if(is.null(models))
        stop("can predict only from 'averaging' object containing model list")
		
	if(is.null(names(models)))
		names(models) <- seq.int(length(models))

	# If all models inherit from lm:
	if (
        (missing(se.fit) || !se.fit) &&
        (is.na(type) || type == "link") &&
        !backtransform &&
        all(linherits(models, c(gam = FALSE, lm = TRUE))) &&
        !anyNA(object$coefficients[1L, ])
		) {
		
		coeff <- coef(object, full = full)

		X <- model.matrix(object)
		if (missing(newdata) || is.null(newdata)) {
			Xnew <- X
		} else {
			tt <- delete.response(terms(formula(object)))
			xlev <- unlist(unname(lapply(models, "[[", "xlevels")),
						   recursive = FALSE, use.names = TRUE)
			Xnew <- model.matrix(tt, data = newdata, xlev = xlev)
		}
		
		colnames(Xnew) <- fixCoefNames(colnames(Xnew))
		Xnew <- Xnew[, match(names(coeff), colnames(Xnew), nomatch = 0L)]
		ret <- (Xnew %*% coeff)[, 1L]
		#if (se.fit) {
		#	scale <- 1
		#	covmx <- solve(t(X) %*% X)
		#	se <- sqrt(diag(Xnew %*% covmx %*% t(Xnew))) * sqrt(scale) ## TODO: use matmult
		#	return(list(fit = y, se.fit = se))
		#}
	} else {
		# otherwise, use brute force:
		if(isFALSE(full)) warning("argument 'full' ignored")

		cl <- as.list(match.call())
		cl$backtransform <- cl$full <- NULL
		cl[[1L]] <- as.name("predict")
		cl <- as.call(cl)

		pred <- lapply(models, function(x, pfr) {
			cl[[2L]] <- x
			y <- tryCatch({
				y <- eval(cl, pfr)
				if(is.numeric(y)) y else structure(as.list(y[c(1L, 2L)]),
					names = c("fit", "se.fit"))
				}, error = function(e) e)
		}, pfr = parent.frame())
		
		err <- sapply(pred, inherits, "condition")
		if (any(err)) {
			lapply(pred[err], warning)
			stop(sprintf(ngettext(sum(err), "'predict' for model %s caused error",
				"'predict' for models %s caused errors"),
				prettyEnumStr(names(models[err]), quote = "'")
				))
		}

        .untransform <- 
        function(fit, se.fit = NULL, models) {
			links <- tryCatch(vapply(models, function(m) family(m)[["link"]], ""),
							error = function(e) NULL)
			if (!is.null(links)) {
				if(any(links[1L] != links[-1L]))
					cry(-1L, "cannot inverse-transform averaged prediction of models using different link functions")
				fam1 <- family(models[[1L]])
				if(is.null(se.fit)) return(fam1$linkinv(fit)) else
					return(list(fit = fam1$linkinv(fit),
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
				par.avg(fit[i, ], se.fit[i, ], weight = Weights(object),
					df = NA_integer_, revised.var = revised.var),
					FUN.VALUE = double(5L)))

			# TODO: ase!
			#no.ase <- all(is.na(object$coefTable[,3]))
			# if(no.ase) 2 else 3
			ret <-  if (backtransform)
					.untransform(apred[1L, ], apred[2L, ], models = models) else
						list(fit = apred[1L, ], se.fit = apred[2L, ])
		} else {
			#tryCatch({
			i <- !vapply(pred, is.numeric, FALSE)
			if(any(i)) pred[i] <- lapply(pred[i], "[[", 1L)
			ret <- apply(do.call("cbind", pred), 1L, weighted.mean,
				w = Weights(object))
			if (backtransform)
				ret <- .untransform(ret, models = models)
		}

	}
	return(ret)
}
