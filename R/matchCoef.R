`matchCoef` <-
function(m1, m2,
	all.terms = getAllTerms(m2, intercept = TRUE),
	beta = 0L,
	terms1 = getAllTerms(m1, intercept = TRUE),
	coef1 = NULL,
	allCoef = FALSE,
	...
	) {
	
	if(is.null(coef1)) {
		ct <- if (beta != 0L) std.coef(m1, beta == 2L, ...) else coefTable(m1, ...)
		coef1 <- ct[, 1L]
		names(coef1) <- rownames(ct)
	} else if(allCoef) stop("'coef1' is given and 'allCoef' is not FALSE")

	if(any((terms1 %in% all.terms) == FALSE)) stop("'m1' is not nested within 'm2'")
	row <- structure(rep(NA_real_, length(all.terms)), names = all.terms)

	fxdCoefNames <- fixCoefNames(names(coef1))
	row[terms1] <- NaN
	pos <- match(terms1, fxdCoefNames, nomatch = 0L)
	row[fxdCoefNames[pos]] <- coef1[pos]
	if(allCoef) {
		i <- match(names(coef1), rownames(ct))
		j <- !is.na(i)
		rownames(ct)[i[j]] <- fxdCoefNames[j]
		attr(row, "coefTable") <- ct
	}
	row
}