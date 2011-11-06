# support for unmarked

`logLik.unmarkedFit` <- function(object, ...) {
  ll <- -object@negLogLike
  attr(ll, "df") <- length(object@opt$par)
  attr(ll, "nobs") <- unmarked::sampleSize(object)
  class(ll) <- "logLik"
  ll
}
#setMethod("logLik", "unmarkedFit", logLik.unmarkedFit)

`formula.unmarkedFit` <- function (x, ...) x@formula

`getAllTerms.unmarkedFit` <- function (x, intercept = FALSE, ...)  {
	f <- formula(x)
	res <- list()
	while(is.call(f) && f[[1L]] == "~") {
		res <- c(f[c(1L, length(f))], res)
		f <- f[[2L]]
	}
	names(res) <- rev(sapply(x@estimates@estimates, slot, "short.name"))

	#res <- lapply(res, function(z) getAllTerms(call("~", z), intercept=FALSE))
	res <- lapply(res, getAllTerms, intercept = FALSE)
	attrInt <- sapply(res, attr, "intercept")
	res <- unlist(lapply(names(res), function(i) sprintf("%s(%s)", i, res[[i]])))
	Ints <- paste(names(attrInt[attrInt != 0L]), "(Int)", sep="")
	if(intercept) res <- c(Ints, res)
	attr(res, "intercept") <- attrInt
	attr(res, "interceptLabel") <-  Ints
	res
}

`tTable.unmarkedFit` <- function (model, ...) {
  do.call("rbind", lapply(model@estimates@estimates, function(y) {
    ret <- cbind(Estimate=y@estimates, SE = sqrt(diag(y@covMat)))
    rn <- rownames(ret)
    rn[rn == "(Intercept)"] <- "Int"
    rownames(ret) <- paste(y@short.name, "(", rn, ")", sep="")
    ret
  }))
}

`coefDf.unmarkedFit` <- function(x) rep(NA, length(coef(x)))
`nobs.unmarkedFit` <- function(object, ...) unmarked::sampleSize(object)