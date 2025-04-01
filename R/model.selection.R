`coefTable.model.selection` <-
function (model, ...) {
	rval <- attr(model, "coefTables")
	names(rval) <- rownames(model)
	rval
}

`coef.model.selection` <-
function (object, ...) {
	ct <- attr(object, "coefTables")
	n <- length(ct)
	allcf <- unique(unlist(lapply(ct, rownames)))
	rval <- matrix(NA_real_, nrow = n, ncol = length(allcf),
		dimnames = list(rownames(object), allcf))
	for(i in seq_len(n))
		rval[i, match(rownames(ct[[i]]), allcf)] <- ct[[i]][, 1L]
	rval
}

`coeffs.model.selection` <-
function (model) coef.model.selection(model)

`coefArray` <- 
function(object) {
	coefNames <- fixCoefNames(unique(unlist(lapply(object, rownames),
		use.names = FALSE)))
	nCoef <- length(coefNames)
	nModels <- length(object)
	rval <- array(NA_real_, dim = c(nModels, 3L, nCoef),
		dimnames = list(names(object), c("Estimate", "Std. Error", "df"), coefNames))
	for(i in seq_along(object)) {
		z <- object[[i]]
		rval[i, seq_len(ncol(z)), ] <- t(z[match(coefNames, fixCoefNames(rownames(z))), ])
	}
	rval
}

`getCall.model.selection` <-
function (x, i = NULL, ...) {
	if(is.null(i)) return(attr(x, "call", exact = TRUE))
	if(length(i) == 1L) return(attr(x, "model.calls", exact = TRUE)[[i]])
	return(attr(x, "model.calls", exact = TRUE)[i])
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
    if (evaluate) eval.parent(cl) else cl
}

`logLik.model.selection` <- function (object, ...) {
	nobs <- attr(object, "nobs")
	n <- nrow(object)
	rval <- vector(n, mode = "list")
	for(i in 1L:n) rval[[i]] <-
		structure(object[i, "logLik"], df = object[i, "df"], nobs = nobs,
			class = "logLik")
	rval
}

`family.model.selection` <-
function (object, ...) {
	if(!is.null(attr(object, "global"))) {
		model.calls <- attr(object, "model.calls")
		if(!is.null(model.calls[[1L]][["family"]])) {
			fam <- lapply(model.calls, "[[", "family")
			rval <- lapply(unique(fam), eval)[
				as.integer(as.factor(vapply(fam, asChar, "")))
				]
			names(rval) <- rownames(object)
			return(rval)
		} else return(family(attr(object, "global")))
	} else {
		attr(object, "model.family")
	}
}

`nobs.model.selection` <-
function (object, ...)
attr(object, "nobs")

## internal: translate column type to column indices
type2col <-
function (x, type) {
    if (inherits(x, "model.selection")) 
        x <- attr(x, "column.types")
	k <- match(x, type, nomatch = 0L)
	i <- k != 0
	which(i)[order(k[i])]
}

## internal: translate column type to column names
type2colname <-
function(x, type)
names(x)[type2col(x, type)] 


`item<-` <- function(x, name, i, value)
`[<-.data.frame`(x, i, name, value)

`item` <- function(x, name, i, ...)
`[.data.frame`(x, i, name, ...)

`itemByType` <- function(x, type, i, ...) 
`[.data.frame`(x, i, type2col(x, type), ...)

`itemByType<-` <- function(x, type, i, value)
`[<-.data.frame`(x, i, type2col(x, type), value)


duplicated.model.selection <-
function (x, incomparables = FALSE, fromLast = FALSE, ...) {
    duplicated.data.frame(x[, type2col(x, c("loglik", "terms"))],
        incomparables = incomparables, fromLast = fromLast, ...)
}
