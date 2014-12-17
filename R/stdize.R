isFALSE <- function(x) identical(FALSE, x)

stdize <-
function(x, ...) UseMethod("stdize")

rootmeansq <- function(v) {
	v <- as.numeric(v[!is.na(v)])
	sqrt(sum(v^2) / max(1, length(v) - 1L))
}

stdize.default <-
stdize.numeric <-
function(x, center = TRUE, scale = TRUE, ...) {
	if(is.function(scale)) {
		scaleFunc <- scale
		scale <- TRUE
	} else scaleFunc <- function(x) sd(x, na.rm = TRUE)
	
	#if(length(list(...))) warning("additional arguments ignored")
	if(length(scale) != 1L) warning("only first element of 'center' is used")
	if(length(center) != 1L) warning("only first element of 'scale' is used")
	scale <- scale[1L]
	center <- center[1L]
	if(is.logical(scale)) scale <- if(scale) scaleFunc(x) else 1
	if(is.logical(center)) center <- if(center) mean(x, na.rm = TRUE) else 0
	x <- (x - center) / scale
	attr(x, "scaled:center") <-  center
	attr(x, "scaled:scale") <- scale
	x
}

stdize.matrix <-
function(x, center = TRUE, scale = TRUE, ...) {
	if(!is.numeric(x)) return(x)
	#if(length(list(...))) warning("additional arguments ignored")
	if(is.function(scale)) {
		scaleFunc <- scale
		scale <- TRUE
	} else scaleFunc <- function(x) sd(x, na.rm = TRUE)
	if(is.logical(scale)) scale <- if(scale) apply(x, 2L, scaleFunc) else 1
	if(is.logical(center)) center <- if(center) colMeans(x, na.rm = TRUE) else 0
	nc <- ncol(x)
	center <- rep(center, length.out = nc)
	scale <- rep(scale, length.out = nc)
	for(i in 1L:nc) x[, i] <- (x[, i] - center[i]) / scale[i]
	attr(x, "scaled:center") <- center
	attr(x, "scaled:scale") <- scale
	x
}

stdize.factor <-
function(x, binary = c("center", "scale", "binary", "half", "omit"),
center = TRUE, scale = FALSE, 
...) {
	#if(length(list(...))) warning("additional arguments ignored")
	if(nlevels(x) == 2L) {
		stdize.logical(as.numeric(x) - 1, binary, center, scale)
	} else x
}

stdize.logical <-
function(x, binary = c("center", "scale", "binary", "half", "omit"),
center = TRUE, scale = FALSE, 
...) {
	#if(length(list(...))) warning("additional arguments ignored")
	binary <- if(is.null(binary) || is.na(binary)) "" else match.arg(binary)
	switch(binary,
			   center = stdize.numeric(x, center = TRUE, scale = 1),
			   scale = stdize.numeric(x, center = TRUE, scale = TRUE),
			   half = stdize.numeric(x, center = 0.5, scale = 1),
			   binary = stdize.numeric(x, center = 0, scale = 1),
			   omit = x,
			   stdize.numeric(as.numeric(x), center = center, scale = scale)
			   )
}

stdize.data.frame <-
function(x, binary = c("center", "scale", "binary", "half", "omit"),
	center = TRUE, scale = TRUE,
	omit.cols = NULL,
	source = NULL, prefix = TRUE, ...) {

	if(is.function(scale)) {
		scaleFunc <- scale
		scale <- TRUE
	} else scaleFunc <- function(x) sd(x, na.rm = TRUE)

	if(!is.null(source)) {
		if(!missing(center) || !missing(scale) || !missing(binary))
			warning("arguments 'center', 'scale' and 'binary' ignored if 'source' is given")
		j <- match(colnames(x), attr(source, "orig.names"))
		if(any(is.na(j))) stop("some columns in 'x' are missing from 'source'")
		center <- attr(source, "scaled:center")[j]
		scale <- attr(source, "scaled:scale")[j]
		binary <- ""
		if(is.null(center) || is.null(scale)) stop("invalid 'source' object")
	} else
		binary <- if(is.null(binary) || is.na(binary)) "" else match.arg(binary)
	
	dataClasses <- vapply(x, function(x) {
		if (is.logical(x)) return("logical")
		if (is.factor(x)) if(nlevels(x) == 2L) return("factor2") else 
			return("other")
		if (is.matrix(x) && is.numeric(x)) return("nmatrix")
		if (is.numeric(x)) return("numeric")
		 return("other")
	}, "")
	
	if(is.character(omit.cols)) dataClasses[colnames(x) %in% omit.cols] <- "omit"
		else if(is.numeric(omit.cols)) dataClasses[omit.cols] <- "omit"
		
	numData <- dataClasses == "numeric"

	if(binary == "omit")  {
		binaryData <- FALSE
	} else {
		binaryData <- dataClasses == "factor2" | dataClasses == "logical"
		for (i in which(binaryData))
			x[, i] <- as.numeric(x[, i]) - if(dataClasses[i] == "factor2") 1 else 0
	}
		
	nc <- ncol(x)

	f <- function(x, bin) {
		if(is.numeric(x)) {
			calc <- isTRUE(bin) & binaryData
			do <- numData | calc | (binaryData & (is.na(bin) | !isFALSE(bin)))
			x[!calc & do & binaryData & !is.na(bin)] <-  bin
			return(list(x, calc, do))
		} else {
			calc <- (numData & x) | (binaryData & ((is.na(bin) & x) | isTRUE(bin)))
			do <- calc | (binaryData & (!is.na(bin) & !isFALSE(bin)))
			num <- numeric(nc)
			num[!calc & do & binaryData & !is.na(bin)] <-  bin
			return(list(num, calc, do))
		}
	}
	
	ctr <- f(center, switch(binary, center = TRUE, scale = TRUE, half = .5, 
		binary = 0, omit = FALSE, NA))
	scl <- f(scale, switch(binary, center = FALSE, scale = TRUE, half = 1, 
		binary = FALSE, omit = FALSE, NA))
	
	center <- ctr[[1L]]
	scale <- scl[[1L]]
	center[ctr[[2L]]] <- colMeans(x[, ctr[[2L]], drop = FALSE], na.rm = TRUE)
	scale[scl[[2L]]] <- apply(x[, scl[[2L]], drop = FALSE], 2L, scaleFunc)
	
	jTransformed <- ctr[[3L]] | scl[[3L]]
	center[jTransformed & !ctr[[3L]]] <- 0
	scale[jTransformed & !scl[[3L]]] <- 1
	for (i in which(jTransformed)) x[, i] <- (x[, i] - center[i]) / scale[i]
	
	attr(x, "scaled:center") <- ifelse(jTransformed, center, NA)
	attr(x, "scaled:scale") <- ifelse(jTransformed, scale, NA)		
	attr(x, "orig.names") <- colnames(x)
	doprefix <-  FALSE
	if(is.character(prefix) ||
	   (doprefix <- (is.logical(prefix) && isTRUE(prefix)))) {
		prefix <- if(doprefix) c("z.", "c.") else rep(prefix, length.out = 2L)
		colnames(x)[jTransformed] <-
			paste0(prefix[jTransformed + (ctr[[3L]] & !scl[[3L]])], 
				colnames(x)[jTransformed])
	}
	x
}

stdize.formula <-
function(x, data = NULL, response = FALSE,
binary = c("center", "scale", "binary", "half", "omit"),
center = TRUE, scale = TRUE,
omit.cols = NULL,
prefix = TRUE, ...) {
	mf <- model.frame(x, data = data, drop.unused.levels = TRUE, ...)
	if(!is.null(omit.cols)) 
		omit.cols <- if(is.character(omit.cols)) 
			which(colnames(mf) %in% omit.cols) else
			stop("'omit.cols' must be a character vector")
	if(!response) omit.cols <- unique(c(omit.cols, 1L))
	attr(mf, "terms") <- NULL
	mf <- stdize.data.frame(mf, center = center, scale = scale, 
		omit.cols = omit.cols, binary = binary, prefix = prefix)
	mf
}
