# Calculate Akaike weights
`Weights` <-
function(x)  UseMethod("Weights")

`Weights.model.selection` <-
function(x) {
	i <- type2col(x, "weight")
	structure(item(x, i) / sum(item(x, i)),	names = row.names(x),
			  wt.type = colnames(x)[type2col(x, "ic")],
			  class = c("model.weights", "numeric"))
}

`Weights.averaging` <-
function(x) {
	rval <- x$msTable[, ncol(x$msTable)]
	class(rval) <- c("model.weights", "numeric")
	attr(rval, "wt.type") <- 
		if(!is.null(attr(x, "model.weights"))) 
			attr(x, "model.weights") else
			asChar(attr(attr(x, "rank"), "call")[[1L]])
	rval
}

`Weights.data.frame` <-
function(x) {
	if(ncol(x) == 2L && colnames(x)[1L] == "df"	&& is.numeric(x[, 2L]))
		return(Weights(x[, 2L]))
	if(ncol(x) == 1L && is.numeric(x[, 1L]))
		return(Weights(x[, 1L]))
	return(NA)
}

`Weights.numeric` <-
function(x) {
	x <- x - min(x)
	d <- exp(-x / 2)
	structure(d / sum(d), class = c("model.weights", "numeric"))
}

`Weights.default` <-
function(x) {
    cry(, "cannot use \"%s\" as 'x'", class(x)[1L])
}

`Weights<-` <-
function(x, value)  UseMethod("Weights<-")


`Weights<-.default` <-
function(x, value) {
	stop("'Weights' can assign weights only to an \"averaging\" object")
}

`Weights<-.averaging` <-
function(x, value) {

	wi <- ncol(x$msTable)
	if(is.null(value)) {
		wts <- Weights(x$msTable[, wi - 1L])
		x$msTable[, wi] <- wts
		colnames(x$msTable)[wi] <- "weight"
		attr(x, "model.weights") <- NULL
	} else {
		x$msTable[, wi] <- value
		wts <- x$msTable[, wi]
		wts <- wts / sum(wts)
		x$msTable[, wi] <- wts
		
		colnames(x$msTable)[wi] <-
			if(inherits(value, "model.weights") && is.character(attr(value, "wt.type")[1L])) {
				paste0(attr(value, "wt.type")[1L], " weight")
			} else 	"[weight]"
			
		attr(x, "model.weights") <-
			if(is.null(attr(value, "wt.type")))
				"unknown" else
				attr(value, "wt.type")
	}

	rv <-  attr(x, "revised.var")
	for(i in 1L:nrow(x$coefficients)) {
		  full <- rownames(x$coefficients)[i] == "full"
		  x$coefficients[i, ] <- .coefarr.avg(x$coefArray, wts, full = full, alpha = 0.05, revised.var = rv)[, 1L]
	}
	
	o <- order(wts, decreasing = TRUE)
	x$msTable <- x$msTable[o, ]
	x$coefArray <- x$coefArray[o,,]
	if(!is.null(attr(x, "modelList"))) attr(x, "modelList") <- attr(x, "modelList")[o]
	x
}



`[.model.weights` <-
function (x, ...) {
	wt.type <- attr(x, "wt.type")
	x <- NextMethod()
	attr(x, "wt.type") <- wt.type
	class(x) <- c("model.weights", class(x))
	x
}

print.model.weights <-
function (x, ...) {
	cat(attr(x, "wt.type"), "model weights", "\n")
	print(format(round(x, 3L), scientific = FALSE), quote = FALSE, right = TRUE)
	invisible(x)
}

