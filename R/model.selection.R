`coefTable.model.selection` <-
function (model, ...) {
	#structure(attr(model, "coefTables"), names = rownames(model))
	ret <- attr(model, "coefTables")
	names(ret) <- rownames(model)
	ret
}

`coef.model.selection` <-
function (object, ...) {
	ct <- attr(object, "coefTables")
	n <- length(ct)
	allcf <- unique(unlist(lapply(ct, rownames)))
	ret <- matrix(NA_real_, nrow = n, ncol = length(allcf),
		dimnames = list(rownames(object), allcf))
	for(i in seq_len(n))
		ret[i, match(rownames(ct[[i]]), allcf)] <- ct[[i]][, 1L]
	ret
}

`coeffs.model.selection` <-
function (model) coef.model.selection(model)

`coefArray` <- function(object) {
	coefNames <- fixCoefNames(unique(unlist(lapply(object, rownames),
		use.names = FALSE)), sort = TRUE)
	nCoef <- length(coefNames)
	nModels <- length(object)
	ret <- array(NA_real_, dim = c(nModels, 3L, nCoef),
		dimnames = list(names(object), c("Estimate", "Std. Error", "df"), coefNames))
	for(i in seq_along(object)) {
		z <- object[[i]]
		ret[i, 1:3, ]
		ret[i, seq_len(ncol(z)), ] <- t(z[match(coefNames, fixCoefNames(rownames(z))), ])
	}
	ret
}


`subset.model.selection` <-
function(x, subset, select, recalc.weights = TRUE, ...) {
	if (missing(select)) {
		if(missing(subset)) return(x)
		e <- .substHas(substitute(subset))
		i <- eval(e, x, parent.frame())
		return(`[.model.selection`(x, i, recalc.weights = recalc.weights, ...))
	} else {
		cl <- match.call(expand.dots = FALSE)
	    cl <- cl[c(1L, match(names(formals("subset.data.frame")), names(cl), 0L))]
	    cl[[1L]] <- as.name("subset.data.frame")
		ret <- eval(cl, parent.frame())
		if(recalc.weights && ("weight" %in% colnames(ret)))
			ret[, 'weight'] <- ret[, 'weight'] / sum(ret[, 'weight'])
	    return(ret)
	}
}

`[.model.selection` <-
function (x, i, j, recalc.weights = TRUE, ...) {
	ret <- `[.data.frame`(x, i, j, ...)
	if (missing(j)) {
		s <- c("row.names", "calls", "coefTables", "random.terms", "order")
		k <- match(dimnames(ret)[[1L]], dimnames(x)[[1L]])
		attrib <- attributes(x)
		attrib[s] <- lapply(attrib[s], `[`, k)
		attributes(ret) <- attrib
		if(recalc.weights)
			ret$weight <- ret$weight / sum(ret$weight)
			#ret[, 'weight'] <- ret[, 'weight'] / sum(ret[, 'weight'])

		if(!is.null(warningList <- attr(ret, "warnings")))
			attr(ret, "warnings") <- warningList[sapply(warningList, attr, "id") %in% rownames(ret)]
	} else {
		cls <- class(ret)
		class(ret) <- cls[cls != "model.selection"] # numeric or data.frame
	}
	return(ret)
}

`print.model.selection` <-
function(x, abbrev.names = TRUE, warnings = getOption("warn") != -1L, ...) {
	orig.x <- x
	if(!is.null(x$weight)) x$weight <- round(x$weight, 3L)
	xterms <- attr(x, "terms")
	if(is.null(xterms)) {
		print.data.frame(x, ...)
	} else {
		if(abbrev.names) xterms <- abbreviateTerms(xterms, 3L)

		colnames(x)[seq_along(xterms)] <-  xterms
		globcl <- attr(x, "global.call")
		if(!is.null(globcl)) {
			cat("Global model call: ")
			print(globcl)
			cat("---\n")
			random.terms <- attr(getAllTerms(attr(x, "global")), "random.terms")
			if(!is.null(random.terms)) random.terms <- list(random.terms)
		} else random.terms <- attr(x, "random.terms")

		cat("Model selection table \n")
		dig <- c("R^2" = 4L, df = 0L, logLik = 3L, AICc = 1L, AICc = 1L,
			AIC = 1L, BIC = 1L, QAIC = 1L, QAICc = 1L, ICOMP = 1L, Cp = 1L,
			delta = 2L,	weight = 3L)

		j <- match(colnames(x), names(dig), nomatch = 0L)
		i <- sapply(x, is.numeric) & (j == 0L)
		x[, i] <- signif(x[, i], 4L)
		for(i in names(dig)[j]) x[, i] <- round(x[, i], digits = dig[i])

		vLegend <- character(0L)
		if(abbrev.names) {
			vCols <- attr(x, "vCols")
			vlen <- nchar(vCols)
			if(!is.null(vCols)) {
				for(i in vCols) {
					z <- x[, i]
					lev <- levels(z)
					lev <- lev[!(lev %in% c("", "NULL"))]
					z <- factor(z, levels = lev)
					n <- nchar(i)
					for(k in seq.int(n, 1L)) {
						shlev <- abbreviateTerms(lev, k, deflate = TRUE)
						if(all(nchar(shlev) <= n)) break;
					}
					x[, i] <- factor(z, labels = shlev)
					if(any(j <- shlev != lev)) vLegend <- c(vLegend, paste(i,
						": ", paste(shlev[j], "=", sQuote(lev[j]),
						collapse = ", "), sep = ""))
				}
			}
		}

		uqran <- unique(unlist(random.terms, use.names = FALSE))
		abbran <- abbreviateTerms(gsub("1 | ", "", uqran, fixed = TRUE), 1L,
			deflate = TRUE)
		colran <- vapply(random.terms, function(s) paste(abbran[match(s, uqran)],
			collapse = "+"), "")

		if(addrandcol <- length(unique(colran)) > 1L) {
			k <- which(colnames(x) == "df")[1L]
			x <- cbind(x[, 1L:(k - 1L)], random = colran, x[, k:ncol(x)])
		}

		print.default(as.matrix(x)[, !sapply(x, function(.x) all(is.na(.x))),
			drop = FALSE], na.print = "", quote = FALSE)

		if(abbrev.names && length(vLegend))
			cat("Abbreviations:", vLegend, sep = "\n")

		if(!is.null(random.terms)) {
			if(addrandcol) {
				cat("Random terms: \n")
				cat(paste(abbran, "=", sQuote(uqran)), sep = "\n")
			} else {
				cat("Random terms (all models): \n")
				cat(paste(sQuote(uqran)), sep = ", ")
				cat("\n")
			}

		}
		if (warnings && !is.null(attr(x, "warnings"))) {
			cat("\n"); print.warnings(attr(x, "warnings"))
		}
	}
	invisible(orig.x)
}

`update.model.selection` <- function (object, global.model, ..., evaluate = TRUE) {
    cl <- attr(object, "call")
    if (is.null(cl)) stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...

	if(!missing(global.model))
		extras <- c(list(global.model = substitute(global.model)), extras)

    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(cl)))
        for (a in names(extras)[existing]) cl[a] <- extras[a]
        if (any(!existing)) {
            cl <- c(as.list(cl), extras[!existing])
            cl <- as.call(cl)
        }
    }
    return(if (evaluate) eval(cl, parent.frame()) else cl)
}

`coef.model.selection` <- function (object, ...)
	object[, attr(object, "terms")]

`logLik.model.selection` <- function (object, ...) {
	nobs <- attr(object, "nobs")
	n <- nrow(object)
	ret <- vector(n, mode = "list")
	for(i in 1:n) ret[[i]] <-
		structure(object[i, "logLik"], df = object[i, "df"], nobs = nobs,
			class = "logLik")
	ret
}
