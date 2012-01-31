# support for unmarked

`logLik.unmarkedFit` <- function(object, ...) {
	ret <- -object@negLogLike
	attr(ret, "df") <- length(object@opt$par)
	attr(ret, "nobs") <- unmarked::sampleSize(object)
	class(ret) <- "logLik"
	return(ret)
}
#setMethod("logLik", "unmarkedFit", logLik.unmarkedFit)

`formula.unmarkedFit` <- function (x, ...) x@formula

`getAllTerms.unmarkedFit` <- function (x, intercept = FALSE, ...)  {
	f <- formula(x)
	ret <- list()
	while(is.call(f) && f[[1L]] == "~") {
		ret <- c(ret, f[c(1L, length(f))])
		f <- f[[2L]]
	}
	ret <- lapply(ret, `environment<-`, NULL)
	names(ret) <- sapply(x@estimates@estimates, slot, "short.name")
	#ret <- lapply(ret, function(z) getAllTerms(call("~", z), intercept=FALSE))
	ret <- lapply(ret, getAllTerms.formula, intercept = FALSE)
	attrInt <- sapply(ret, attr, "intercept")
	#ret <- unlist(lapply(names(ret), function(i) sprintf("%s(%s)", i, ret[[i]])))
	ret <- unlist(lapply(names(ret), function(i) if(length(ret[[i]])) paste(i, "(", ret[[i]], ")",
		sep = "") else character(0L)))
	Ints <- paste(names(attrInt[attrInt != 0L]), "(Int)", sep = "")
	if(intercept) ret <- c(Ints, ret)
	attr(ret, "intercept") <- attrInt
	attr(ret, "interceptLabel") <-  Ints
	return(ret)
}

## tweak for 'distsamp' models: prefix the detection "p(...)" terms with 'sigma'
getAllTerms.unmarkedFitDS  <- function (x, intercept = FALSE, ...)  {
	tt <- getAllTerms.unmarkedFit(x, intercept = FALSE)
	ret <- gsub("^p\\(", "p(sigma", c(tt))
	intLab <- attr(tt, "interceptLabel")
	intLab[intLab == "p(Int)"] <- "p(sigma(Intercept))"
	if(intercept) ret <- c(intLab, ret)
	mostattributes(ret) <- attributes(tt)
	attr(ret, "interceptLabel") <- intLab
	ret
}

`coefTable.unmarkedFit` <- function (model, ...)
	.makeCoefTable(coef(model), sqrt(diag(vcov(model))))

`nobs.unmarkedFit` <- function(object, ...) unmarked::sampleSize(object)

coeffs.unmarkedFit <- function(model) {
	ret <- lapply(model@estimates@estimates, coef, altNames = FALSE)
	pfx <- rep(vapply(model@estimates@estimates, slot, "", "short.name"),
		vapply(ret, length, 1L))
	ret <- unlist(unname(ret))
	Ints <- which(names(ret) == "Int")
	names(ret) <- paste(pfx, "(", names(ret), ")", sep = "")
	attr(ret, "Intercept") <- Ints
	ret
}
