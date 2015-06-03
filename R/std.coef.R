`beta.weights` <-
function(model) {
	.Deprecated("std.coef")
	std.coef(model, FALSE)
}

.vif <-
function(x) {
	v <- vcov(x)
	if(dim(v)[1L] < 2L) return(structure(rep_len(1, dim(v)[1L]),
										 names = dimnames(v)[[1L]]))
	if ((ndef <- sum(is.na(coef(x)))) > 0L) 
        stop(sprintf(ngettext(ndef, "one coefficient is not defined",
			"%d coefficients are not defined"), ndef))
	o <- attr(model.matrix(x), "assign")
	if (any(int <-(o == 0))) {
		v <- v[!int, !int, drop = FALSE]
	} else warning("no intercept: VIFs may not be sensible")
	R <- cov2cor(v)
	detR <- det(R)
	vfs <- numeric(m <- dim(v)[1L])
	for (j in 1L:m)
		vfs[j] <- det(as.matrix(R[j, j])) * det(as.matrix(R[-j,-j])) / detR
	
	rval <- numeric(length(int))
	names(rval) <- dimnames(vcov(x))[[1L]]
	rval[!int] <- vfs
	rval[int] <- 1
	rval
}

.partialsd <-
function(x, sd, vif, n, p = length(x) - 1) {
	sd * sqrt(1 / vif) * sqrt((n - 1)/(n - p))
}

partial.sd <- function(x) {
	mm <- model.matrix(x)
	.partialsd(coef(x), apply(mm, 2L, sd), .vif(x), nobs(x), sum(attr(mm, "assign") != 0))
}

`std.coef` <-
function(x, partial.sd, ...) {
	#coefmat <- coefTable(x, ...)[, 1L:2L, drop = FALSE]
	coefmat <- coefTable(x, ...)
	mm <- model.matrix(x)
	if(partial.sd) {
		bx <- .partialsd(coefmat[, 1L], apply(mm, 2L, sd),
				.vif(x), nobs(x), sum(attr(mm, "assign") != 0))
	} else {
		response.sd <- sd(model.response(model.frame(x)))
		bx <- apply(mm, 2L, sd) / response.sd
	}
	coefmat[, 1L:2L] <- coefmat[, 1L:2L] * bx
	colnames(coefmat)[1L:2L] <- c("Estimate*", "Std. Error*")
	return (coefmat)
}


