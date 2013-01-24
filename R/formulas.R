# test for marginality constraints
`formulaMargChk` <-
function(frm, except = NULL) {
	if(isTRUE(except)) return(TRUE)
	factors <- attr(terms(frm), "factors")
	if(length(factors) == 0L) return(TRUE)
	ex <- rownames(factors)[apply(factors, 1L, function(x) any(x > 1L))]
	if(is.character(except))
		factors <- factors[!(rownames(factors) %in% except), ]
	ret <- all(factors < 2L)
	attr(ret, "marg.ex") <- ex
	return(ret)
}

`.formulaEnv` <- function(object, env = .GlobalEnv) {
	res <- as.formula(object, env = env)
	environment(res) <- env
	res
}


`simplify.formula` <- function(x) {
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

`expand.formula` <- function(x) {
	x <- formula(x)
	tt <- terms(x)
	x[[length(x)]] <- reformulate(attr(tt, "term.labels"),
		intercept = attr(tt,"intercept"))[[2L]]
	x
}

