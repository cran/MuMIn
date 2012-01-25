# This is merely to get rid of the annoying behaviour in summary.glmML.
# it does not do anything except for printing the model output.
`summary.glmmML` <- function(object, ...) object

# family
`family.default` <- function (object, ...)  {
	cl <- getElement(object, "call")
	if(is.null(cl)) return(NULL)
	fam <- cl$family
	if(is.null(fam)) fam <- formals(match.fun(cl[[1L]]))$family
	if(is.null(fam)) return(gaussian())
	switch(mode(fam), call = eval(fam), name =, character = match.fun(fam)())
}


`family.gls` <-
`family.lme` <-
stats:::family.lm

`nobs.rq` <-
function (object, ...) length(object$y)


# coxme:

`logLik.coxme` <-
function(object, type = c("integrated", "penalized"), ...) {
	type <- match.arg(type)
	i <- which(type == c("integrated", "penalized"))[1L]
	ret <- object$loglik[[i + 1L]]
	attr(ret, "df") <- object$df[i]
	attr(ret, "nobs") <- object$n[1L]
	class(ret) <- "logLik"
	ret
}

`logLik.lmekin` <-
function(object, ...) {
	ret <- object$loglik
	attr(ret, "nobs") <- object$n
	attr(ret, "df") <- length(object$coefficients$fixed) +
		length(object$coefficients$random) + 1L
	class(ret) <- "logLik"
	ret
}

`nobs.coxme` <-
`nobs.lmekin` <-
function (object, ...) object$n[1L]

`getAllTerms.coxme` <-
function(x, ...)  {
	ret <- MuMIn:::getAllTerms.terms(terms(x))
	random <- x$formulaList$random
	attr(ret, "random.terms") <- as.character(random)
	f <- as.name(".")
	for(f1 in random) f <- call("+", f, f1)
	attr(ret, "random") <- call("~", as.name("."), f)
	attr(ret, "intercept") <- 0L
	attr(ret, "interceptLabel") <- NULL
	ret
}

`formula.coxme` <-
function(x, ...)  {
	ret <- x$formulaList$fixed
	f <- ret[[3L]]
	for(f1 in x$formulaList$random) f <- call("+", f, f1)
	ret[[3L]] <- f
	ret
}

`formula.lmekin` <-
function(x, ...) eval(x$call$formula, parent.frame())


`coeffs.coxme` <-
`coeffs.lmekin` <-
function(model) {
	# for class coxme:
	ret <- model$coefficients
	# for class lmekin and older coxme
	if(is.list(ret) && !is.null(ret$fixed)) return(ret$fixed)
	ret
}


`makeArgs.coxme` <-
`makeArgs.lmekin` <-
function(obj, termNames, comb, opt, ...) {
	ret <- makeArgs.default(obj, termNames, comb, opt)
	f <- .Internal(update.formula(as.formula(ret$formula), as.formula(. ~ . + 1)))
	ret$formula <- update.formula(f, opt$random)
	ret
}
