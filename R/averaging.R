`coef.averaging` <-
function(object, full = FALSE, ...) {
	## XXX: backward compatibility:
	object <- upgrade_averaging_object(object)
	full <- .checkFull(object, full) 
	object$coefficients[if(full) 1L else 2L, ]
}

`fitted.averaging` <-
function (object, ...) .NotYetImplemented()
# predict.averaging(object, backtransform = TRUE, type = "link")

model.frame.averaging <-
function (formula, ...) {
	mergeMF(getModelList(formula))
}

model.matrix.averaging <-
function (object, ...) {
	if(j <- match("x", names(object), nomatch = 0L)) return(object[[j]])
	
	mf <- model.frame(object)
	do.call("model.matrix", list(object = terms(mf), data = mf,
		contrasts.arg = get.contrasts(object)))
}

`summary.averaging` <-
function (object, ...) {
	## XXX: backward compatibility:
	object <- upgrade_averaging_object(object)

	.makecoefmat <- function(cf) {
		no.ase <- all(is.na(cf[, 3L]))
		z <- abs(cf[, 1L] / cf[, if(no.ase) 2L else 3L])
		pval <- 2 * pnorm(z, lower.tail = FALSE)
		cbind(cf[, if(no.ase) 1L:2L else 1L:3L, drop = FALSE],
			`z value` = z, `Pr(>|z|)` = zapsmall(pval))
	}
	
	is.arm <- ncol(object$msTable) == 6L && (colnames(object$msTable)[6L] == "ARM weight")

	weight <- object$msTable[, if(is.arm) 6L else 5L]
	
	object$coefmat.full <- .makecoefmat(.coefarr.avg(object$coefArray, weight,
			attr(object, "revised.var"), TRUE, 0.05))

	if(!is.arm) object$coefmat.subset <-
		.makecoefmat(.coefarr.avg(object$coefArray, weight,
			attr(object, "revised.var"), FALSE, 0.05))
	
	object$coef.nmod <- colSums(!is.na(object$coefArray[, 1L, , drop = FALSE]))

	structure(object, ARM = is.arm, class = c("summary.averaging", "averaging"))
}

`confint.averaging` <-
function (object, parm, level = 0.95, full = FALSE, ...) {
	## XXX: backward compatibility:
	object <- upgrade_averaging_object(object)
	full <- .checkFull(object, full) 


    a2 <- 1 - level
    a <- a2 / 2
    cf <- object$coefArray[, 1L, ]
    pnames <- colnames(cf)
    if (missing(parm)) parm <- pnames
		else if (is.numeric(parm)) parm <- pnames[parm]
	missing.par <- is.na(cf)
    se <- object$coefArray[, 2L, ]
    dfs <- object$coefArray[, 3L, ]
	if(full) {
		se[missing.par] <- cf[missing.par] <- 0
		if(!all(is.na(dfs))) dfs[missing.par] <- Inf
	}
    wts <- Weights(object) ## XXX: !
    ci <- t(sapply(parm, function(i)
		par.avg(cf[,i], se[,i], wts, dfs[, i], alpha = a2)))[, 4L:5L]
    
	
	ci[is.na(object$coefficients[1L, ]), ] <- NA_real_
    colnames(ci) <- getFrom("stats", "format.perc")(c(a, 1L - a), 3L)
    return(ci)
}

`print.summary.averaging` <-
function (x, digits = max(3L, getOption("digits") - 3L),
    signif.stars = getOption("show.signif.stars"), ...) {

    cat("\nCall:\n", paste(asChar(x$call, nlines = -1L), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
		
	comcallstr <- 
	if(!is.null(attr(x, "model.calls"))) {
		commonCallStr(calls = attr(x, "model.calls"))
	} else if(!is.null(attr(x, "modelList"))) {
		commonCallStr(attr(x, "modelList"))
	} else NA
	if(!is.na(comcallstr)) {
		cat("Component model call: \n")
		cat(strwrap(comcallstr), sep = " \n     ")
	}		
		
    cat("\nComponent models: \n")
	msTable <- x$msTable
	wi <- ncol(msTable)
	if(!isTRUE(attr(x, "ARM")) && names(msTable)[wi] != "weight")
		msTable <- msTable[, c(1L, wi), drop = FALSE]
		
	print(round(as.matrix(msTable), 2L), na.print = "")

	if(!is.null(attr(x$msTable, "term.codes"))) {
		cat("\nTerm codes: \n")
		print.default(attr(x$msTable, "term.codes"), quote = FALSE)
	}
	cat("\nModel-averaged coefficients: ")
	if (nnotdef <- sum(is.na(x$coefmat.full[, 1L]))) {
		 msg <- paste0("\n(", nnotdef, " not defined because of singularities in all ",
			"component models)", collapse = "")
		cat(strwrap(msg, exdent = 4L), sep = "\n")
	}

		 
	hasPval <- TRUE
	coefTitles <- if(isTRUE(attr(x, "ARM")))
		c(coefmat.full = "(ARM average)") else
		c(coefmat.full = "(full average)",
		  coefmat.subset = "(conditional average)")
		
	n <- length(coefTitles)	
	for (i in seq.int(n)) {
	    iname <- names(coefTitles[i])
	    if (is.null(x[[i]])) next
	    cat(" \n", coefTitles[i], " \n", sep = "")
	    printCoefmat(x[[iname]],
	        P.values = hasPval, has.Pvalue = hasPval,
	        digits = digits, signif.stars = signif.stars,
	        signif.legend = i == n
	    )
	}
	
	nose <- apply(x$coefArray[, 2L, ], 1L, function(x) all(is.na(x)))
	msg <- if(all(nose)) "Standard errors cannot be calculated because no component models provide them \n" else
		if(any(nose)) "Standard errors cannot be calculated because some component models do not provide them \n"
	cat(strwrap(msg, exdent = 4L), sep = "\n")
	
	#if (no.ase) cat("Confidence intervals are unadjusted \n")
	#printCoefmat(matrix(x$coef  .shrinkage, nrow = 1L,
		#dimnames = list("", x$term.names)), P.values = FALSE,
		#has.Pvalue = FALSE, cs.ind = seq_along(x$term.names), tst.ind = NULL)

	# cat("\nSum of weights: \n")
	# print(round(x$sw, 2L))
}

`print.averaging` <-
function(x, ...) {
	## XXX: backward compatibility:
	x <- upgrade_averaging_object(x)
    cat("\nCall:\n", paste(asChar(x$call, nlines = -1L), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat("Component models:", "\n")
	comp.names <- rownames(x$msTable)
	comp.names[comp.names == ""] <- "null"
	cat(format(sQuote(comp.names), justify = "l"), fill = TRUE)
	cat("\nCoefficients:", "\n")
	print(x$coefficients[!is.na(x$coefficients[,1L]), , drop = FALSE])
    x
}

`vcov.averaging` <- 
function (object, full = FALSE, ...) {
	## XXX: backward compatibility:
	object <- upgrade_averaging_object(object)
	full <- .checkFull(object, full) 

	full <- as.logical(full)[1L]

	models <- attr(object, "modelList")
	if(is.null(models)) stop("cannot calculate covariance matrix from ",
							 "'averaging' object without component models")

	vcovs <- lapply(lapply(models, vcov), as.matrix)
	names.all <- dimnames(object$coefArray)[[3L]]
	nvars <- length(names.all)
	nvarseq <- seq(nvars)
	wts <- Weights(object)
	wts <- wts / sum(wts) # normalize just in case

	vcov0 <- matrix(if(full) 0 else NA_real_, nrow = nvars,
		ncol = nvars, dimnames = list(names.all, names.all))

	vcovs2 <- lapply(vcovs, function(v) {
		i <- match(fixCoefNames(dimnames(v)[[1L]]), names.all)
		vcov0[i, i] <- v
		return(vcov0)
	})
	b1 <- object$coefArray[, 1L, ]
	if(full) b1[is.na(b1)] <- 0
	avgb <- object$coefficients[2L - full, ]
	
	res <- sapply(nvarseq, function(c1) sapply(nvarseq, function(c2) {
		 weighted.mean(sapply(vcovs2, "[", c1, c2) + (b1[, c1] - avgb[c1]) *
			(b1[, c2] - avgb[c2]), wts, na.rm = TRUE)
	}))
	dimnames(res) <- list(names.all, names.all)
	return(res)
}

`logLik.averaging` <- function (object, ...) {
	models <- attr(object, "modelList")
	if(is.null(models)) {
		nobs <- attr(object, "nobs")
		apply(object$msTable, 1L, function(x) structure(list(x[2L]),
			df = x[1L], nobs = nobs, class = "logLik"))
	} else {
		structure(lapply(attr(object, "modelList"), logLik),
			names = rownames(object$msTable))
	}
}

`coefTable.averaging` <-
function (model, full = FALSE, adjust.se = TRUE, ...) {
	full <- .checkFull(model, full) 
	
    no.ase <- any(is.na(model$coefArray[,3L,]) & !is.na(model$coefArray[,1L,]))
	if(!missing(adjust.se) && adjust.se && no.ase) 
        warning("adjusted std. error not available for this type of model")
		
	weight <- model$msTable[, ncol(model$msTable)]

    cols <- c(1L, if(!adjust.se || no.ase) 2L else 3L)
	ct <- .coefarr.avg(model$coefArray, weight, TRUE, full, .05)[, cols, drop = FALSE] 
	.makeCoefTable(ct[,1L], ct[,2L], NA, rownames(ct))
}
