`beta.weights` <-
function(model) {
	.Deprecated("std.coef")
	std.coef(model, FALSE)
}

.vif <-
function(x) {
	v <- .vcov(x)
	nam <- dimnames(v)[[1L]]
	if(dim(v)[1L] < 2L) return(structure(rep_len(1, dim(v)[1L]),
										 names = dimnames(v)[[1L]]))
	if ((ndef <- sum(is.na(coef(x)))) > 0L) 
        stop(sprintf(ngettext(ndef, "one coefficient is not defined",
			"%d coefficients are not defined"), ndef))
	o <- attr(model.matrix(x), "assign")
	if (any(int <-(o == 0))) {
		v <- v[!int, !int, drop = FALSE]
	} else warning("no intercept: VIFs may not be sensible")
	d <- sqrt(diag(v))
	rval <- numeric(length(nam))
	names(rval) <- nam
	rval[!int] <- diag(solve(v / (d %o% d)))
	rval[int] <- 1
	rval
}

.partialsd <-
function(x, sd, vif, n, p = length(x) - 1) {
	sd * sqrt(1 / vif) * sqrt((n - 1) / (n - p))
}

partial.sd <-
function(x) {
	b <- coef(x)
	mm <- model.matrix(x)
	mm <- mm[, match(names(b), colnames(mm)), drop = FALSE]
	colnames(mm) <- names(b)
	m <- ncol(mm)
	.partialsd(b, apply(mm, 2L, sd), .vif(x), nobs(x),
			   sum(attr(mm, "assign") != 0))
}

`std.coef` <-
function(x, partial.sd, ...) {
	#b <- coefTable(x, ...)[, 1L:2L, drop = FALSE]
	b <- coefTable(x, ...)
	mm <- model.matrix(x)
	mm <- mm[, match(rownames(b), colnames(mm)), drop = FALSE]
	colnames(mm) <- names(b)
	#b <- b[colnames(mm), ]
	if(partial.sd) {
		bx <- .partialsd(b[, 1L], apply(mm, 2L, sd),
				.vif(x), nobs(x), sum(attr(mm, "assign") != 0))
	} else {
		response.sd <- sd(model.response(model.frame(x)))
		bx <- apply(mm, 2L, sd) / response.sd
	}
	b[, 1L:2L] <- b[, 1L:2L] * bx
	colnames(b)[1L:2L] <- c("Estimate*", "Std. Error*")
	return(b)
}
