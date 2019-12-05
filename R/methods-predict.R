## Predict methods for objects for which they are not available in their
## original packages, or replacements.

# Add 'se.fit' argument for predict:
# https://stat.ethz.ch/pipermail/r-help/2004-April/050144.html
# http://web.archiveorange.com/archive/v/rOz2zbtjRgntPMuIDoIl

# based on the original 'predict.gls' in package 'nlme'

`predict.gls` <-
function (object, newdata, se.fit = FALSE, na.action = na.fail, ...) {
    if (missing(newdata) && !se.fit) {
        return(fitted(object))
    }
    form <- getFrom("nlme", "getCovariateFormula")(object)
    mfArgs <- list(formula = form, data = newdata, na.action = na.action)
    mfArgs$drop.unused.levels <- TRUE
    dataMod <- do.call(model.frame, mfArgs)
    contr <- object$contrasts
    for (i in names(dataMod)) {
        if (inherits(dataMod[, i], "factor") && !is.null(contr[[i]])) {
            levs <- levels(dataMod[, i])
            levsC <- dimnames(contr[[i]])[[1]]
            if (any(wch <- is.na(match(levs, levsC)))) {
                stop(sprintf(ngettext(sum(wch), "level %s not allowed for %s", 
                  "levels %s not allowed for %s"), paste(levs[wch], collapse = ",")), domain = NA)
            }
            attr(dataMod[, i], "contrasts") <- contr[[i]][levs, , drop = FALSE]
        }
    }
    N <- nrow(dataMod)
    if (length(all.vars(form)) > 0) {
        X <- model.matrix(form, dataMod)
    } else {
        X <- array(1, c(N, 1), list(row.names(dataMod), "(Intercept)"))
    }
    cf <- coef(object)
    val <- c(X[, names(cf), drop = FALSE] %*% cf)
	if(se.fit) {
		se <- sqrt(matmultdiag(X %*% vcov(object), ty = X))
		val <- list(fit = val, se.fit = unname(se))
	}	
    lab <- "Predicted values"
    if (!is.null(aux <- attr(object, "units")$y)) {
        lab <- paste(lab, aux)
    }
    structure(val, label = lab)
}



`predict.lme` <-
function (object, newdata, level, asList = FALSE,
	na.action = na.fail, se.fit = FALSE, ...) {
	cl <- match.call()
	cl$se.fit <- NULL
	cl[[1L]] <- call("get", "predict.lme", asNamespace("nlme"))	
	res <- eval.parent(cl)
	
	if(se.fit && (missing(level) || any(level > 0)))
		warning("cannot calculate standard errors for level > 0")
	if(se.fit && !missing(level) && length(level) == 1L && all(level == 0)) {
		if (missing(newdata) || is.null(newdata)) {
			X <- model.matrix(object, data = object$data)
		} else {
			tt <- delete.response(terms(formula(object)))
			xlev <- .getXlevels(tt, model.frame(object, data = object$data))
			X <- model.matrix(tt, data = newdata, contrasts.arg =
							  object$contrasts, xlev = xlev)
		}
		se <- sqrt(matmultdiag(X %*% vcov(object), ty = X))
		# se <- sqrt(rowSums((X %*% vcov(object)) * X))
		# se <- sqrt(diag(X %*% vcov(object) %*% t(X))) ## TODO: use matmult
		names(se) <- names(res)
		list(fit = c(res), se.fit = se)
	} else res
}

.predict_glm <-
function (object, newdata, type = c("link", "response"), se.fit = FALSE,
	trms, coeff, offset, ...) {
	
    type <- match.arg(type)
    if (!missing(newdata) && !is.null(newdata)) {
        xlev <- .getXlevels(trms, model.frame(trms, data = newdata))
        X <- model.matrix(trms, data = newdata, contrasts.arg = attr(model.matrix(object), 
            "contrasts"), xlev = xlev)
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(trms, "offset"))) 
            for (i in off.num) offset <- offset + eval(attr(trms, 
                "variables")[[i + 1L]], newdata)
		
		cl <- get_call(object)
        if (!is.null(cl$offset)) 
            offset <- offset + eval(cl$offset, newdata)
    } else {
        X <- model.matrix(object)
        if (!length(offset))  offset <- NULL
    }
    y <- (X %*% coeff)[, 1L]
    if (!is.null(offset)) 
        y <- y + offset
    fam <- family(object)
    if (se.fit) {
        # covmat <- as.matrix(vcov(object))
		se <- sqrt(matmultdiag(X %*% as.matrix(vcov(object)), ty = X))
		# se <- sqrt(rowSums((X %*% covmat) * X))
        # se <- sqrt(diag(X %*% covmat %*% t(X)))
        if (type == "response" && inherits(fam, "family")) 
            list(fit = fam$linkinv(y), se.fit = se * abs(fam$mu.eta(y)))
        else list(fit = y, se.fit = se)
    } else {
        if (type == "response" && inherits(fam, "family")) 
            fam$linkinv(y)
        else y
    }
}

`predict.gamm` <- 
function (object, ...) mgcv::predict.gam(object[['gam']], ...)
