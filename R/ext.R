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

`coefTable.rq` <- function(model, ...)
	.makeCoefTable(model$coefficients, rep(NA_real_, length(model$coefficients)))


# Classes 'coxme' and 'lmekin' from package 'coxme':

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
function (object, ...) object$n[2L]

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


## Classes 'hurdle' and 'zeroinfl' from package 'pscl':

`nobs.hurdle` <-
`nobs.zeroinfl` <- `nobs.lmekin`

`getAllTerms.hurdle` <- function(x, intercept = FALSE, ...) {
	f <- as.formula(formula(x))
	# to deal with a dot in formula (other classes seem to expand it)
	if("." %in% all.vars(f))
		getAllTerms.terms(terms.formula(f, data = eval(x$call$data, envir =
			environment(f))), intercept = intercept)
	else getAllTerms.formula(f, intercept = intercept)
}

`getAllTerms.zeroinfl` <- function(x, intercept = FALSE, ...) {
	f <- formula(x)
	if(length(f[[3L]]) != 1L && f[[3L]][[1L]] == "|"){
		f1 <- call("~", f[[2L]], f[[3L]][[2L]])
		f2 <- call("~", f[[3L]][[3L]])
	} else {
		f1 <- f
		f2 <- NULL
	}
	fs <- lapply(lapply(c(f1, f2), terms.formula, data = eval(x$call$data)),
		formula)
	z <- lapply(fs, getAllTerms, intercept = TRUE)

	ord <- unlist(lapply(z, attr, "order"))
	n <- sapply(z, length)
	if(length(z) > 1L) ord[-j] <- ord[-(j <- seq_len(n[1L]))] + n[1L]
	zz <- unlist(z)
	Ints <- which(zz == "(Intercept)")
	#zz[Ints] <- "1"
	#zz <- paste(rep(c("count", "zero")[seq_along(z)], sapply(z, length)),
		#"(", zz, ")", sep = "")
	zz <- paste(rep(c("count", "zero")[seq_along(z)], sapply(z, length)),
		"_", zz, sep = "")
	ret <- if(!intercept) zz[-Ints] else zz
	attr(ret, "intercept") <- pmin(Ints, 1)
	attr(ret, "interceptLabel") <- zz[Ints]
	attr(ret, "response") <- attr(z[[1L]], "response")
	attr(ret, "order") <- if(!intercept) order(ord[-Ints]) else ord
	ret
}

`coefTable.zeroinfl` <-
function(model, ...)
	.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))

`coefTable.hurdle` <- function(model, ...) {
	cts <- summary(model)$coefficients
	ct <- do.call("rbind", unname(cts))
	cfnames <- paste(rep(names(cts), vapply(cts, nrow, 1L)), "_", rownames(ct),
		sep = "")
	.makeCoefTable(ct[, 1L], ct[, 2L], coefNames = cfnames)
	#.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))
}
