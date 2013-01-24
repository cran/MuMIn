`getAllTerms.default` <-
#function(x, ...) getAllTerms.formula(as.formula(formula(x)), ...)
function(x, ...) getAllTerms.terms(terms(as.formula(formula(x))), ...)

`getAllTerms.gam` <-
function(x, intercept = FALSE, ...)
	getAllTerms.terms(terms(formula(x), ...), intercept = intercept)

`getAllTerms.lm` <-
function(x, intercept = FALSE, ...)
	getAllTerms.terms(terms(x, ...), intercept = intercept)

`getAllTerms.terms` <-
function(x, offset = TRUE, intercept = FALSE, ...) {

	interceptLabel <- "(Intercept)"
	variables <- attr(x, "variables")[-1L]

	if (!is.null(attr(x, "offset"))){
		offs <- sapply(variables[attr(x, "offset")], deparse)
	} else offs <- NULL

	ret <- attr(x, "term.labels")

	# Get term names, with higher order term components arranged alphabetically
	if (length(ret) > 0) {
		factors <- attr(x, "factors")
		factors <- factors[order(rownames(factors)), , drop = FALSE]
		v <- rownames(factors)
		ret <- apply(factors != 0L, 2L, function(x) paste(v[x], collapse = ":"))
	}

	# Leave out random terms (lmer type)
	#ran <- attr(x, "variables")[-1][-c(attr(x, "offset"), attr(x, "response"))]
	ran <- variables
	ran <- as.character(ran[sapply(ran,
		function(x) length(x) == 3L && x[[1L]] == "|")])
	ifx <- !(ret %in% ran)

	ret <- ret[ifx] # ifx - indexes of fixed terms
	#retUnsorted <- ret

	# finally, sort by order and then alphabetically
	#ret <- unname(ret[order(attr(x, "order")[ifx], ret)])
	ord <- order(attr(x, "order")[ifx], gsub("I\\((.*)\\)", "\\1", ret))
	ret <- unname(ret[ord])

	if(intercept && attr(x, "intercept")) {
		ret <- c(interceptLabel, ret)
		ord <- c(1, ord + 1L)
	}

	if (!is.null(offs[1L])) {
		if (offset) {
			ret <- c(ret, offs)
			ord <- c(ord, length(ord) + 1L)
		}
		attr(ret, "offset") <- offs
	}
	attr(ret, "intercept") <- attr(x, "intercept")
	attr(ret, "interceptLabel") <- interceptLabel

	if (length(ran) > 0L) {
		attr(ret, "random.terms") <- ran
		attr(ret, "random") <- reformulate(c(".", paste("(", ran, ")",
			sep = "")), response = ".")
	}

	response <- attr(x, "response")
	response <- if(response == 0L) NULL else variables[[response]]
	attr(ret, "response") <- response
	attr(ret, "order") <- order(ord)

	return(ret)
}

`getAllTerms.formula` <-
function(x, ...) getAllTerms.terms(terms.formula(x), ...)

`getAllTerms.lme` <-
function(x, ...) {
	ret <- getAllTerms.terms(terms(x), ...)
	attr(ret, "random") <- . ~ .

	# Code from nlme:::print.reStruct, modified slightly
	reStruct <- x$modelStruct$reStruct
	nobj <- length(reStruct)
	if (is.null(namx <- names(reStruct)))
		names(reStruct) <- nobj:1L
	aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
	aux[lower.tri(aux)] <- ""
	reStruct[] <- rev(reStruct)
	aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
	aux[lower.tri(aux)] <- ""
	attr(ret, "random.terms") <- paste(lapply(lapply(reStruct, attr, "formula"),
		"[[", 2L), "|",
		rev(apply(aux, 1L, function(z) paste(z[z != ""], collapse = " %in% "))))

	return(ret)
}

# `getAllTerms.glmer` <- # For backwards compatibility
# `getAllTerms.lmer` <-  # with older versions of lme4
`getAllTerms.mer` <-
function(x, ...) getAllTerms(lme4::formula(x), ...)

# Apparently there is no (explicit) intercept in coxph, but 'terms' gives
# attr(,"intercept") == 1.
`getAllTerms.coxph` <- function (x, ...) {
	ret <- getAllTerms.default(x, ...)
	attr(ret, "intercept") <- 0L
	attr(ret, "interceptLabel") <- NULL
	return(ret)
}

`getAllTerms.glmmML` <- function (x, ...) {
	ret <- getAllTerms.terms(terms(x), ...)
	#attr(ret, "random.terms") <-  deparse(call("|", 1,  x$call$cluster))
	attr(ret, "random.terms") <-  paste("1 |",  x$call$cluster)
	return(ret)
}

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

`getAllTerms.glimML` <- function(x, intercept = FALSE, ...) {
	ret <- getAllTerms.default(x, intercept = intercept, ...)
	ttran <- terms.formula(x@random)
	ran <- attr(ttran, "term.labels")
	if(length(ran)) attr(ret, "random.terms") <- paste("1 |", ran)
	ret
}

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
`getAllTerms.unmarkedFitDS` <- function (x, intercept = FALSE, ...)  {
	tt <- getAllTerms.unmarkedFit(x, intercept = FALSE)
	ret <- gsub("^p\\(", "p(sigma", c(tt))
	intLab <- attr(tt, "interceptLabel")
	intLab[intLab == "p(Int)"] <- "p(sigma(Intercept))"
	if(intercept) ret <- c(intLab, ret)
	mostattributes(ret) <- attributes(tt)
	attr(ret, "interceptLabel") <- intLab
	ret
}

`getAllTerms.MCMCglmm` <- 
function (x, ...) {
	res <- MuMIn:::getAllTerms.default(x, ...)
	attr(res, "random") <- .formulaEnv(.~., environment(formula(x)))
	attr(res, "random.terms") <- deparse(x$Random$formula, control = NULL)[1]
	res
}

`getAllTerms.gamm` <-
function (x, ...) getAllTerms(x$gam, ...)

`getAllTerms.mark` <- 
function (x, intercept = FALSE, ...) {
	
	f <- formula(x, expand = FALSE)[[2L]]
	ret <- list()
	while(length(f) == 3L && f[[1L]] == "+") {
		ret <- append(f[[3L]], ret)
		f <- f[[2L]]
	}
	ret <- append(f, ret)
	res <- lapply(ret, function(x) {
		func <- deparse(x[[1L]], control = NULL)
		tt <- terms(eval(call("~", x[[2L]])))
		tlab <- attr(tt, "term.labels")
		torder <- attr(tt, "order")
		if(attr(tt, "intercept")) {
			tlab <- append("(Intercept)", tlab)
			torder <- c(0L, torder)
		}
		res1 <- lapply(fixCoefNames(tlab), function(z) paste(func, "(", z, ")", sep = ""))
		attr(res1, "order") <- torder
		res1
	})
	
	ord <- order(rep(seq_along(res), sapply(res, length)),
		unlist(lapply(res, attr, "order")))
	res <- unlist(res, recursive = TRUE)[ord]
	ints <- grep("((Intercept))", res, fixed = TRUE)
	attr(res, "intercept") <- as.numeric(ints != 0L)
	attr(res, "interceptLabel") <- res[ints]
	if(!intercept) {
		res <- do.call("structure", c(list(res[-ints]), attributes(res)))
		attr(res, "order") <- order(ord[-ints])
	} else {
		attr(res, "order") <- order(ord)
	}
	
	res
}


`getAllTerms` <-
function(x, ...) UseMethod("getAllTerms")

print.allTerms <- function(x, ...) {
	cat("Model terms: \n")
	if(!length(x)) {
		cat("<None>\n")
	} else {
	print.default(as.vector(x),quote = TRUE)
	}
	ints <- attr(x, "interceptLabel")
	if(!is.null(ints)) {
		cat(ngettext(n = length(ints), "Intercept:", "Intercepts:"), "\n")
		print.default(ints,quote = TRUE)
	}
}
