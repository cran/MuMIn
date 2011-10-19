# Mallow's Cp

`Cp` <- 
function(object, dispersion = NULL) {
	rss <- deviance(object)
	df.r <- df.residual(object)
	scale <- if (!is.null(dispersion)) dispersion else if 
		(family(object)$family %in% c("poisson", "binomial")) 1 else if
		(df.r > 0) rss / df.r else 
		NaN
	rss + 2 * scale  * (nobs(object) - df.r)
}