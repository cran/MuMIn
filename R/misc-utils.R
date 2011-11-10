`DebugPrint` <- function(x) { cat(deparse(substitute(x)), "= \n") ; print(x) }
`srcc` <- function() {
	ret <- eval(expression(source("clipboard", local = TRUE)), .GlobalEnv)
	return(if(ret$visible) ret$value else invisible(ret$value))
}


#if (!exists("getElement", mode = "function", where = "package:base", inherits = FALSE)) {
getElement <- function (object, name) {
    if (isS4(object))
		if (.hasSlot(object, name)) slot(object, name) else NULL
    else object[[name, exact = TRUE]]
}
#}

# cbind list of data.frames omitting duplicated column (names)
`cbindDataFrameList` <-
function(x) {
	dfnames <- unlist(lapply(x, colnames))
	uq <- !duplicated(dfnames)
	res <- do.call("cbind", x)[,uq]
	colnames(res) <- dfnames[uq]
	return(res)
}

# same for rbind, check colnames and add NA's when any are missing
`rbindDataFrameList` <-
function(x) {
	all.colnames <- unique(unlist(lapply(x, colnames)))
	x <- lapply(x, function(y) {
		y[all.colnames[!(all.colnames %in% colnames(y))]] <- NA
		return(y[all.colnames])
	})
	return(do.call("rbind", x))
}

`videntical` <-
function(x) all(vapply(x[-1L], identical, logical(1), x[[1L]]))

# Check class for each object in a list
`linherits` <- function(x, whats) {
	as.logical(vapply(x, inherits, integer(length(whats)), names(whats),
		which=TRUE)) == whats
}

# substitute has(a, !b, ...) for !is.na(a) & is.na(b) ..., in expression
`.substHas` <- function(e) {
	if(is.expression(e)) e <- e[[1L]]
	n <- length(e)
	if(n == 1L) return(e)
	if(e[[1L]] != "has") {
		for(i in 1:n) e[[i]] <- .substHas(e[[i]])
		return(e)
	}
	res <- NULL
	for(i in seq.int(2L, n)) {
		ex <- if(length(e[[i]]) == 2L && e[[i]][[1L]] == "!")
			call("is.na", e[[i]][[2L]]) else
			call("!", call("is.na", e[[i]]))
		res <- if(i == 2L) ex else call("&", res, ex)
	}
	res <- call("(", res)
	return(res)
}
