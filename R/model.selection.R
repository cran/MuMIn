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
		use.names = FALSE)))
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

`getCall.model.selection`  <-
function (x, i = NULL, ...) {
	if(is.null(i)) return(attr(x, "call", exact = TRUE))
	if(length(i) == 1L) return(attr(x, "model.calls", exact = TRUE)[[i]])
	return(attr(x, "model.calls", exact = TRUE)[i])
}

`subset.model.selection` <-
function(x, subset, select, recalc.weights = TRUE, recalc.delta = FALSE, ...) {
	expr.sub.expand <- expression(
		.substFunc(.substFunc(.substFunc(
			substitute(subset),
			"dc", .sub_dc_has, as.name(".subset_vdc")),
			c("{", "Term"), .sub_Term),
			"has", .sub_has))
					
	subst <- function(cl, ...) eval(call("substitute", cl, list(...)))

	if (missing(select)) {
		if(missing(subset)) return(x)
		e <- eval(expr.sub.expand)
		e <- subst(e, . = x)
		DebugPrint(e)
		i <- eval(e, x, parent.frame())
		return(`[.model.selection`(x, i, recalc.weights = recalc.weights, 
			recalc.delta = recalc.delta, ...))
	} else {
		cl <- match.call(expand.dots = FALSE)
		if(!missing(subset)) cl$subset <- 
			subst(eval(expr.sub.expand), . = substitute(x))
	    cl <- cl[c(1L, match(names(formals("subset.data.frame")), names(cl), 0L))]
	    cl[[1L]] <- as.name("subset.data.frame")
		DebugPrint(cl)
		ret <- eval(cl, parent.frame())
		if(recalc.weights && ("weight" %in% colnames(ret)))
			ret[, 'weight'] <- ret[, 'weight'] / sum(ret[, 'weight'])
		if(recalc.delta && ("delta" %in% colnames(ret)))
			ret[, 'delta'] <- ret[, 'delta'] - min(ret[, 'delta'])
	    return(ret)
	}
}

`[.model.selection` <-
function (x, i, j, recalc.weights = TRUE, recalc.delta = FALSE, ...) {
	ret <- `[.data.frame`(x, i, j, ...)
	if (missing(j)) {
		s <- c("row.names", "model.calls", "coefTables", "random.terms", "order")
		k <- match(dimnames(ret)[[1L]], dimnames(x)[[1L]])
		attrib <- attributes(x)
		attrib[s] <- lapply(attrib[s], `[`, k)
		attributes(ret) <- attrib
		if(recalc.weights) {
			ret[, 'weight'] <- Weights(`[.data.frame`(ret, , 
				which(names(ret) == "delta") - 1L))
		}
		if(recalc.delta) {
			ic <- `[.data.frame`(ret, , which(names(ret) == "delta") - 1L)
            ret[, "delta"] <- ic - min(ic)
		}
		if(!is.null(warningList <- attr(ret, "warnings")))
			attr(ret, "warnings") <- warningList[sapply(warningList, attr, "id") %in% rownames(ret)]
	} else {
		cls <- class(ret)
		class(ret) <- cls[cls != "model.selection"] # numeric or data.frame
	}
	return(ret)
}

getModelClass <-
function(x) {
	if(inherits(x, "model.selection")) {
		if(!is.null(attr(x, "global"))) return(class(attr(x, "global"))[1L])
		if("class" %in% colnames(x)) return(as.character(x[, "class"]))
		if(!is.null(attr(x, "model.class"))) return(attr(x, "model.class"))
	}
	return(NULL)
}

# TODO: named arguments...
`rbind.model.selection` <- 
function (..., deparse.level = 1) {
	allargs <- list(...)
	n <- length(allargs) 
	if(n == 1L) {
		return(allargs[[1L]])
	} else if(n > 1L) {		
		idx <- seq(1L, n)
		nm <- split(make.unique(unlist(lapply(allargs, row.names))),
					rep(idx, sapply(allargs, nrow)))		
		for(i in idx)
			row.names(allargs[[i]]) <- nm[[i]]	
		res <- allargs[[1L]]
		for(i in seq(2L, n))
			res <- merge(res, allargs[[i]], suffixes = NULL)
		return(res)
	} else return(NULL)
}

`merge.model.selection` <-
function (x, y, suffixes = c(".x",".y"), ...)  {
	
	a1 <- attributes(x)
	a2 <- attributes(y)
	if(!identical(a1$rank.call, a2$rank.call))
		stop("models not ranked by the same IC")
	if(!identical(a1$nobs, a2$nobs))
		stop("models fitted to different number of observations")
	c1 <- c(a1$terms, a1$vCols)
	c2 <- c(a2$terms, a2$vCols)
	res <- cbind(rbindDataFrameList(list(x[, c1, drop = FALSE], y[, c2, drop = FALSE])),
				 rbindDataFrameList(list(x[, !(colnames(x) %in% c1), drop = FALSE],
										 y[, !(colnames(y) %in% c2), drop = FALSE])))
	
	if(!is.null(suffixes))
		row.names(res) <- c(paste0(row.names(x), suffixes[1L]),
			 paste0(rownames(y), suffixes[2L]))

	nm <- rownames(res)
	
	newattr <- list()
	for(i in c("model.calls", "coefTables"))
		newattr[[i]] <- structure(c(a1[[i]], a2[[i]]), names = nm)
	for(i in c("rank", "rank.call", "nobs", "class"))
		newattr[[i]] <- a1[[i]]
	for(i in c("terms", "vCols"))
		newattr[[i]] <- unique(c(a1[[i]], a2[[i]]))
	attr(newattr[["terms"]], "interceptLabel") <-
		unique(c(attr(a1$terms,"interceptLabel"),
				 attr(a2$terms,"interceptLabel")))
		
	mclsx <- getModelClass(x)
	mclsy <- getModelClass(y)
	if(length(mclsx) != 1L || !identical(mclsx, mclsy)) {
		if(!("class" %in% colnames(res))) {
			res[, "class"] <- NA_integer_
			pos <- length(c(newattr$terms, newattr$vCols))
			nc <- ncol(res)
			res <- res[, c(1L:pos, nc, (pos + 1L):(nc - 1L))]
		}
		res[, "class"] <- as.factor(c(rep(mclsx, length.out = nrow(x)),
							rep(mclsy, length.out = nrow(y))))
	} else {
		newattr[["model.class"]] <- mclsx
	}
	
	for(i in names(newattr)) attr(res, i) <- newattr[[i]]
	class(res) <- c("model.selection", "data.frame")
		
	o <- order(res[, which(colnames(res) == "delta") - 1L])
	res <- res[o, recalc.delta = TRUE]
	res
}

`print.model.selection` <-
function(x, abbrev.names = TRUE, warnings = getOption("warn") != -1L, ...) {
	orig.x <- x
	if(!is.null(x$weight)) x$weight <- round(x$weight, 3L)
	xterms <- attr(x, "terms")
	if(is.null(xterms) || !all(xterms %in% colnames(x)[seq_along(xterms)])) {
		print.data.frame(x, ...)
	} else {
		if(abbrev.names) xterms <- abbreviateTerms(xterms, 6L, 3L, deflate = TRUE)

		colnames(x)[seq_along(xterms)] <- xterms
		globcl <- attr(x, "global.call")
		if(!is.null(globcl)) {
			cat("Global model call: ")
			print(globcl)
			cat("---\n")
			random.terms <- attr(getAllTerms(attr(x, "global")), "random.terms")
			if(!is.null(random.terms)) random.terms <- list(random.terms)
		} else random.terms <- attr(x, "random.terms")
		cat("Model selection table \n")
		dig <- c(AnyIC = 1L, "R^2" = 4L, df = 0L, logLik = 3L,
			delta = 2L,	weight = 3L)

		j <- match(colnames(x), names(dig), nomatch = 0L)
		iic <- length(j) - 2L
		j[iic] <- 1L # AnyIC
		names(dig)[1L] <- colnames(x)[iic]
		i <- sapply(x, is.numeric) & (j == 0L)
		
		x[, i] <- signif(x[, i], 4L)
		for(i in names(dig)[j]) x[, i] <- round(x[, i], digits = dig[i])

		vLegend <- NULL
		if(abbrev.names) {
			vCols <- attr(x, "vCols")
			vCols <- vCols[(vCols %in% colnames(x)) & !(vCols %in% c("class"))]
			vlen <- nchar(vCols)
			vLegend <- vector(length(vCols), mode = "list")
			names(vLegend) <- vCols

			if(!is.null(vCols)) {
				for(i in vCols) {
					lev <- levels(x[, i])
					lev <- lev[!(lev %in% c("", "NULL"))]
					shlev <- abbreviateTerms(lev, nchar(i), deflate = TRUE)
					x[, i] <- factor(x[, i], levels = lev, labels = shlev)
					if(any(j <- shlev != lev)) vLegend[[i]] <-
						paste(shlev[j], "=", sQuote(lev[j]))
				}
				vLegend <- vLegend[!vapply(vLegend, is.null, TRUE)]
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

		if(abbrev.names && length(vLegend)) {
			cat("Abbreviations:", sep = "\n")
			for(i in names(vLegend)) {
				cat(vLegend[[i]], sep = ", ", fill = TRUE, labels =
					c(paste0(i, ":"), rep(paste(rep(" ", nchar(i) + 1L),
					collapse = ""), length(vLegend[[i]]) - 1L)))
			}
		}
		
		cat("Models ranked by", deparse(attr(attr(x, 'rank'), "call"), control = NULL), "\n")

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

`logLik.model.selection` <- function (object, ...) {
	nobs <- attr(object, "nobs")
	n <- nrow(object)
	ret <- vector(n, mode = "list")
	for(i in 1:n) ret[[i]] <-
		structure(object[i, "logLik"], df = object[i, "df"], nobs = nobs,
			class = "logLik")
	ret
}

`$<-.model.selection` <- function (x, name, value) {
	ret <- base::`$<-.data.frame`(x, name, value)
	if(name %in% attr(x, "terms")) class(ret) <- "data.frame"
	ret
}

`row.names<-.model.selection` <- 
function (x, value) {
	x <- `row.names<-.data.frame`(x, value)
	names(attr(x, "coefTables")) <-
	names(attr(x, "model.calls")) <- value
	x
}

`family.model.selection` <-
function (object, ...) {
	if(!is.null(attr(object, "global"))) {
		model.calls <- attr(object, "model.calls")
		if(!is.null(model.calls[[1L]][["family"]])) {
			fam <- lapply(model.calls, "[[", "family")
			fam1 <- unique(fam)
			ret <- lapply(unique(fam), eval)[
				as.integer(as.factor(vapply(fam, deparse, "", control = NULL)))
				]
			names(ret) <- rownames(object)
			#index <- split(seq_along(fam), vapply(fam, deparse, "", control = NULL))
			#for(i in seq_along(fam1)) fam1[[i]] <- list(family = eval(fam1[[i]]), index = index[[i]])
			#fam <- family(dd1)
			#index <- lapply(fam, "[[", "index")
			#ret <- rep(lapply(fam, "[[", "family"), vapply(index, length, 1L))[order(unlist(index))]
			return(ret)
		} else return(family(attr(object, "global")))
	} else {
		return(attr(object, "model.family"))
	}
}

#### XXX
.argTable <-
function(cl, family = NULL, class = NULL,
		 args.omit = NULL, different.only = FALSE) {
	haveNoCall <-  vapply(cl, is.null, FALSE)
	cl[haveNoCall] <- lapply(cl[haveNoCall], function(x) call("<unknown>", formula = NA))
	arg <- lapply(cl, function(x) sapply(x, function(argval)
		switch(mode(argval), character = , logical = argval,
		numeric = signif(argval, 3L), deparse(argval, nlines = 1L))))
	arg <- rbindDataFrameList(lapply(lapply(arg, t), as.data.frame))
	if(!is.null(args.omit)) arg <- arg[, !(colnames(arg) %in% args.omit)]

	arg[] <- lapply(arg, as.factor)
	
	if(!is.null(family)) {
		.getFam <- function(x) unlist(x[c("family",	"link")])
		fam <-  if(inherits(family, "family"))
			matrix(.getFam(family), dimnames = list(c("family", "link"), NULL)) else
			sapply(family, .getFam)
		f <- fam[1L, ]
		f[is.na(f)] <- ""
		f <- vapply(strsplit(f, "(", fixed = TRUE), "[", "", 1L)
		f[f == "Negative Binomial"] <- "negative.binomial"
		fam[2L, fam[2L, ] == vapply(unique(f), function(x) if(is.na(x))
									NA_character_ else formals(get(x))$link,
									FUN.VALUE = "")[f]] <- NA_character_
		j <- !is.na(fam[2L,])
		famname <- fam[1L, j]
		famname <- ifelse(substring(famname, nchar(famname)) != ")",
			paste0(famname, "("), paste0(substring(famname, 1L, nchar(famname) - 1L),
				", "))
		fam[1L, j] <- paste0(famname, fam[2L, j], ")")
		arg <- cbind(arg, t(fam))
	}
	if(!is.null(class)) arg[, "class"] <- rep(class, length.out = nrow(arg))

	#arg <- as.matrix(arg)
	#arg[is.na(arg) | arg == "NULL"] <- ""
	colnames(arg)[1L] <- "FUN"
	for (i in seq_len(ncol(arg))) {
		v <- arg[, i]
		if(any(j <- is.na(v) | v == "NULL" | v == "")) {
			levels(v) <- c(levels(v), "")
			v[j] <- ""
			arg[, i] <- v
		}
	}
	
	if(different.only)
		arg <- arg[, vapply(arg, nlevels, 1L) != 1L, drop = FALSE]

	#if(ncol(arg) != 0L) arg <- gsub("([\"'\\s]+|\\w+ *=)","", arg, perl = TRUE)
	arg
}



