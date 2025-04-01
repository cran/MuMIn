
.MuMInEnv <- new.env(parent = baseenv())

testStart <- 
function(...) {
    p <- c(...)
	res <- length(find.package(p, quiet = TRUE)) == length(p)
    if(res) {
		for(a in p) suppressPackageStartupMessages(
			library(a, character.only = TRUE, quietly = TRUE))
		require(.packageName, character.only = TRUE)
		options(na.action = "na.fail")
	}
	res
}

warnonce <- 
function(..., show.instance = 0L) {
	id <- make.names(deparse1(match.call(expand.dots = FALSE)$...))
	count <- get0(flag <- paste0("warned.", as.character(id)[1L]), .MuMInEnv,
					ifnotfound = 0L)
	if(count <= show.instance)
		assign(flag, count + 1L, envir = .MuMInEnv)
	if(count == show.instance) {
		cl <- match.call()
		cl$show.instance <- NULL
		cl[[1L]] <- as.name("warning")
		eval.parent(cl)
	}
}

`cry` <-
function(Call = NA, Message, ..., warn = FALSE, domain = paste0("R-", .packageName)) {
	if (is.character(Call)) {
		Call <- call(Call[1L], sys.call(-1L)[[1L]])
	} else if(is.numeric(Call)) {
		Call <- sys.call(Call - 1L)
	} else if (!is.call(Call) && !is.null(Call))
		Call <- sys.call(-1L)
	if(warn) warning(simpleWarning(gettextf(Message, ..., domain = domain), Call)) else
		stop(simpleError(gettextf(Message, ..., domain = domain), Call))
}


`getElement` <- function (object, name) {
    if (isS4(object))
		if (.hasSlot(object, name)) slot(object, name) else NULL
    else object[[name, exact = TRUE]]
}

# cbind list of data.frames omitting duplicated column (names)
`cbindDataFrameList` <-
function(x) {
	dfnames <- unlist(lapply(x, colnames), use.names = FALSE)
	uq <- !duplicated(dfnames)
	res <- do.call("cbind", x)[,uq]
	colnames(res) <- dfnames[uq]
	return(res)
}

# same for rbind, check colnames and add NA's when any are missing
`rbindDataFrameList` <-
function(x) {
	all.colnames <- unique(unlist(lapply(x, colnames), use.names = FALSE))
	x <- lapply(x, function(y) {
		y[all.colnames[!(all.colnames %in% colnames(y))]] <- NA
		return(y[all.colnames])
	})
	return(do.call("rbind", x))
}

`videntical` <-
function(x) all(vapply(x[-1L], identical, logical(1L), x[[1L]]))

# Check class for each object in a list
`linherits` <- function(x, whats) {
	as.logical(vapply(x, inherits, integer(length(whats)), names(whats),
		which = TRUE)) == whats
}

# tries to make a list of element names
`.makeListNames` <- function(x) {
	nm <- names(x)
	lapply(seq_along(x), function(i) {
		if(is.null(nm) || !nzchar(nm[i])) {
			switch(mode(x[[i]]),
				call = {
						v <- asChar(x[[i]], width.cutoff = 20L)
						if(length(v) != 1L) v <- sprintf("%s...", v[1L])
						v },
				symbol =, name = as.character(x[[i]]),
				NULL =, logical =, numeric =, complex =, character = x[[i]], i
				)
		} else nm[i]
	})
}

# test if dependency chain is satisfied: x[n] can be TRUE only if x[n+1] are also TRUE
`.subset_dc` <- function(...) {
	n <- length(x <- c(...))
	if(n > 1L) all(x[-n] >= x[-1L]) else TRUE
}

# vectorized version of .subset_do (used within subset.model.selection)
`.subset_vdc` <- function(...) apply(cbind(..., deparse.level = 0L), 1L, .subset_dc)

`prettyEnumStr` <- function(x, sep = ", ", sep.last = gettext(" and "), quote = TRUE) {
	n <- length(x)
	if(is.function(quote))
		x <- quote(x) else {
			if(identical(quote, TRUE)) quote <- '"'
			if(is.character(quote)) x <- paste0(quote, x, quote)
		}
	paste0(x, if(n > 1L) c(rep(sep, n - 2L), sep.last, "") else NULL,
		collapse = "")
}

# `splitList` <- function (x, k) {
    # n <- length(x)
    # ret <- unname(split.default(x, findInterval(seq_len(n), seq(0L, n +
        # 1L, length = k + 1L))))
	# if(k > n) ret <- c(ret, vector(k - n, mode = "list"))
	# ret
# }


`.parallelPkgCheck` <- function(quiet = FALSE) {
	# all this is to trick the R-check
	if(!("snow" %in% loadedNamespaces())) {
		if(getRversion() < "2.14.0") {
			if(length(find.package("snow", quiet = TRUE)))
				do.call("require", list("snow"))
		} else if(length(find.package("parallel", quiet = TRUE)))
			do.call("require", list("parallel", quiet = TRUE))
	}
	if(!exists("clusterCall", mode = "function")) {
		if(quiet) return(FALSE) else
			stop("cannot find function 'clusterCall'")
	} else return(TRUE)
}

`clusterVExport` <- local({
   `getv` <- function(obj, env = as.environment(1L))
		for (i in names(obj)) assign(i, obj[[i]], envir = env)
	function(cluster, ...) {
		Call <- match.call()
		Call$cluster <- NULL
		Call <- Call[-1L]
		vars <- list(...)
		vnames <- names(vars)
		if (is.null(vnames)) {
			names(vars) <- vapply(Call, asChar, "")
		} else if (any(!nzchar(vnames))) {
			names(vars) <- ifelse(!nzchar(vnames), vapply(Call, asChar, ""), vnames)
		}
		get("clusterCall")(cluster, getv, vars)
		# clusterCall(cluster, getv, vars)
	}
})

# test if 'x' can be updated (in current environment or on a cluster)
# level is 0/FALSE - no checking, 1 - check if variables and functions exist,
# >1 - reevaluate x and compare with original 
`testUpdatedObj` <- function(cluster = NA, x, call = get_call(x),
	level = 1L, exclude = "subset") {
	
	if(isTRUE(level)) level <- 2L

	if (level > 0L) {
		xname <- asChar(substitute(x))
		doParallel <- inherits(cluster, "cluster")
		if(doParallel) {
			clusterCall <- get("clusterCall")
			whereStr <- gettext(" in the cluster nodes' environment")
			csapply <- function(...) clusterCall(cluster, "sapply", ...)
		} else {
			whereStr <- ""
			csapply <- function(...) sapply(...)
		}
		if(is.null(call)) stop(gettextf("'%s' has no call component", xname))
		call.orig <- call
		if(!is.null(call$data)) {
			# get rid of formulas, as they are evaluated within 'data'
			call <- call[!sapply(call, function(x) "~" %in% all.names(x))]
			call[exclude] <- NULL
		}	
		
		v <- all.vars(call, functions = FALSE)
		if(!all(z <- unlist(csapply(v, "exists", where = 1L)))) {
			z <- unique(names(z[!z]))
			stop(sprintf(ngettext(length(z), "variable %s not found%s",
				"variables %s not found%s"), prettyEnumStr(z, quote = "'"), whereStr))
			}
		vfun <- all.vars(call, functions = TRUE)
		if(!all(z <- unlist(csapply(vfun[!(vfun %in% v)], "exists",
			mode = "function", where = 1L)))) {
			zz <- unique(names(z[!z]))
			stop(sprintf(ngettext(length(zz), "function %s not found%s",
				"functions %s not found%s"), prettyEnumStr(zz, quote = "'"), whereStr))
			}
		if(level > 1L && !missing(x)) {
			if(doParallel) {
				# XXX: Import: clusterCall
				if(!all(vapply(lapply(clusterCall(cluster, eval, call.orig), all.equal, x), isTRUE, TRUE)))
					stop(gettextf("'%s' evaluated on the cluster nodes differs from the original one",
				xname))
			} else if (!isTRUE(all.equal(x, update(x))))
				stop(gettextf("updated '%s' differ(s) from the original one", xname))
		}
	}
}

`tryCatchWE` <- function (expr) {
	Warnings <- NULL
	list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
		warning = function(w) {
			Warnings <<- c(Warnings, list(w))
			invokeRestart("muffleWarning")
		}), warnings = Warnings)
}

# like apply(, 2) but returns a list (does not do any checking)
`applyrns` <- function (X, FUN, ...) {
	n <- nrow(X)
	ret <- vector(n, mode = "list")
	for(i in seq_len(n)) if(!is.null(z <- FUN(X[i, ], ...))) ret[[i]] <- z
	ret
}


## from stats:::format_perc
`format_perc` <-
function (probs, digits) 
paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), 
    "%")

return_null <-
function(...) NULL 

## Cheating RCheck:
getFrom <-
function(pkg, name)
get(name, envir = asNamespace(pkg), inherits = FALSE)

# used by 'model.sel' and 'dredge' with argument 'extra'
.get.extras <-
function(extra, r2nullfit = NULL) {
	
	extraExpr <- substitute(extra)
	
	if(!is.vector(extra)) {
		extraExpr <- call("alist", extraExpr)
		extra <- list(extra)
	}
	
	isfun <- vapply(extra, is.function, NA)
	if(any(isfun)) {
		if(all(isfun) && !is.null(names(extra)) && all(nzchar(names(extra))))
			return(extra)
		
		extraExpr[[1L]] <- as.name("alist")
		extra <- eval.parent(extraExpr)
	}
	
	extraNames <- names(extra) %||% character(length(extra))
	emptynames <- !nzchar(extraNames)
	if(any(emptynames)) {
		extraNames[emptynames] <-
			do.call("c", .mapply(function(x, i) switch(mode(x),
				call = asChar(x[[1L]]), name = asChar(x),
					`function` = paste0("function", i),
					character = , as.character(x)[1L]),
				list(x = extra[emptynames], i = which(emptynames)), MoreArgs = list()))
	}
	
	extra <- as.list(extra)
	names(extra) <- extraNames
	if(anyDuplicated(extra)) {
		ok <- !duplicated(extra)
		extra <- extra[ok]
	}
	if(any(i <- vapply(extra, is.language, TRUE)))
		extra[i] <- lapply(extra[i], eval.parent)
    pos <- match(c("R^2", "adjR^2"), extra, nomatch = 0L)
	if(any(pos != 0L)) {
        nullfit_ <- NULL # to pass R check
		if(!is.null(r2nullfit)) {
            r2env <- new.env()
            assign("nullfit_", r2nullfit, envir = r2env)  
            if(pos[1L] != 0L) {
                f <-  function(x) r.squaredLR(x, null = nullfit_)
                environment(f) <- r2env
                extra[["R^2"]] <- f
            }
            if(pos[2L] != 0L) {
                f <-  function(x) attr(r.squaredLR(x, null = nullfit_), "adj.r.squared")
                environment(f) <- r2env
                extra[["adjR^2"]] <- f
            }
		} else {
            if(pos[1L] != 0L) extra[["R^2"]] <- function(x) r.squaredLR(x)
            if(pos[2L] != 0L) extra[["adjR^2"]] <- function(x) attr(r.squaredLR(x), "adj.r.squared")
		}
	}
	sapply(extra, match.fun, simplify = FALSE)
}


.applyExtras <- 
function(model, extra) {
    rval <-.mapply(function(f, a, model) tryCatch({
        z <- f(model)
        mode(z) <- "numeric"
        z
    }, error = function(err) {
                err$call  <- call(a, quote(submodel))
                err$message <- sprintf("while evaluating \"extra\" function '%s': '%s'", a, err$message)
                stop(err)
            }), list(f = extra, a = names(extra)), MoreArgs = list(model = model))
    names(rval) <- names(extra)
    unlist(rval)
}

## matrix multiplication with option of calculating the diagonal only
## It is more memory efficient and faster than `crossprod` for large matrices
matmult <-
function(x, y, diag.only = FALSE) {
	if(ncol(x) != nrow(y)) stop('non-conformable arguments')
	n1 <- nrow(x)
	n2 <- ncol(y)
	if(diag.only) {
		if(n1 != n2) stop('non-conformable arguments')
		## >2x faster:
		return(rowSums(x * t(y)))
		#res <- numeric(n1)
		#for(i in seq.int(n1)) res[i] <- sum(x[i, ] * y[, i])
	} else {
		res <- matrix(nrow = n1, ncol = n2)
		for(i in seq.int(n1)) for(j in seq.int(n2))
			res[i, j] <- sum(x[i, ] * y[, j])
		return(res)
	}
}

## matmultdiag(x, ty = y) == diag(x %*% t(y))
matmultdiag <-
function(x, y, ty = t(y)) {
	if(ncol(x) != ncol(ty)) stop('non-conformable arguments')
	if(nrow(x) != nrow(ty)) stop('result is not a square matrix')
	return(rowSums(x * ty))
}


tmpvarname <-
function(envir, n = 8L) {
	while(exists(x <- paste0(c("*", sample(letters, n), "*"),
		collapse = ""), envir)) {}
	x
}

.lab2expr <-
function(x) {
	x <- gsub(":", "%*%", x, perl = TRUE)
	x <- gsub("\\B_?(\\d+)(?![\\w\\._])", "[\\1]", x, perl = TRUE)
    x <- gsub("((?<=[,=]) +(?=[\\w\"'])|(?<=[\\w\"']) +(?==))", "", x, perl = TRUE)
    x <- gsub("[ _]", "~~", x)
	x <- str2expression(x)
	x[] <- lapply(x, function(x)
		if(is.call(x) && x[[1L]] == "I"  && length(x) == 2L)
			x[[2]] else x)
	x
}


coefmatch <-
function(x, y) {
    match(fixCoefNames(names(coeffs(x))),
    fixCoefNames(names(coeffs(y))))
}

