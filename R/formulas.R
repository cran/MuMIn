# test for marginality constraints
`formulaMargChk` <-
function(frm, except = NULL) {
	
	if(isTRUE(except)) return(TRUE)
	factors <- attr(terms.formula(frm, simplify = FALSE), "factors")
	if(length(factors) == 0L) return(TRUE)
	
	#benchmark({
	#X <- factors
	#n <- nrow(X)
	#res <- vector(mode = "logical", n)
	#for(i in 1L:n) res[i] <- any(X[i, ] > 1L)
	#})
	#benchmark(rowSums(factors > 1L) != 0L)
	#benchmark(apply(factors > 1L, 1L, any))
	#benchmark(apply(factors, 1L, function(x) any(x > 1L)))
	
	ex <- dimnames(factors)[[1L]][rowSums(factors > 1L) != 0L]
	if(is.character(except))
		factors <- factors[!(dimnames(factors)[[1L]] %in% except), ]
	ret <- all(factors < 2L)
	attr(ret, "marg.ex") <- ex
	return(ret)
}

# slightly faster than stats::reformulate
# response must be a character string
#Reformulate <- function(termlabels, response = NULL, intercept = TRUE, envir = parent.frame()) {
#	res <- parse(text = paste(if(!is.null(response)) 'Y', "~", paste(termlabels, collapse = "+"), collapse = ""))[[1L]]
#	class(res) <- "formula"
#	environment(res) <- envir
#	res
#}

`.formulaEnv` <- 
function(object, env = .GlobalEnv) {
	res <- formula(object, env = baseenv())
	environment(res) <- env
	res
}


`simplify.formula` <- 
function(x) {
	tt <- terms(as.formula(x))
	fac <- attr(tt, "factors")
	if(length(fac) == 0L) {
		x[[length(x)]] <- if(attr(tt, "intercept")) 1 else -1
		return(x)
	}
	if(ncol(fac) == 1L) {
		tnm <- attr(tt, "term.labels")
	} else {
		ord <- attr(tt, "order")
		k <- seq_along(colnames(fac))
		names(k) <- colnames(fac)
		k <- k[order(ord, decreasing = TRUE)]
		ret <- sapply(k, function(i) sapply(k, function(j)
			if(ord[j] >= ord[i]) NA else !any(!(fac[, i] == 1L) & fac[, j])
		))
		i <- (!apply(ret, 1L, function(x) any(x, na.rm = TRUE)))
		j <- i & apply(fac[, k], 2L, function(x) all(x < 2L)) &
			 ord[k] > 1
		tnm <- rownames(ret)
		tnm[j] <- gsub(":", "*", tnm[j])
		tnm <- tnm[i][order(ord[k][i])]
	}
	x[[length(x)]]  <- reformulate(tnm, intercept = attr(tt, "intercept"))[[2L]]
	return(x)
}

`expand.formula` <- 
function(x) {
	x <- formula(x)
	env <- environment(x)
	tt <- terms(x)
	x[[length(x)]] <- reformulate(attr(tt, "term.labels"),
		intercept = attr(tt,"intercept"))[[2L]]
	environment(x) <- env
	x
}


