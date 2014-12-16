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
	if (length(ret) > 0L) {
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

	# finally, sort by term order and then alphabetically
	#ret <- unname(ret[order(attr(x, "order")[ifx], ret)])
	ord <- order(attr(x, "order")[ifx], gsub("I\\((.*)\\)", "\\1", ret))
	ret <- unname(ret[ord])

	deps <- if (length(ret) > 0L) termdepmat(reformulate(ret)) else
		matrix(FALSE, 0L, 0L)
		
	dimnames(deps) <- list(ret, ret)
	diag(deps) <- NA
	
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
		f.random <- reformulate(c(".", paste0("(", ran, ")")), response = ".")
		environment(f.random) <- environment(x)
		attr(ret, "random") <- f.random
	}

	response <- attr(x, "response")
	response <- if(response == 0L) NULL else variables[[response]]
	attr(ret, "response") <- response
	attr(ret, "order") <- order(ord)
	attr(ret, "deps") <- deps


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
function(x, ...) getAllTerms(.xget("lme4", "formula")(x), ...)

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

#`getAllTerms.hurdle` <- function(x, intercept = FALSE, ...) {
#	f <- as.formula(formula(x))
#	# to deal with a dot in formula (other classes seem to expand it)
#	if("." %in% all.vars(f))
#		getAllTerms.terms(terms.formula(f, data = eval(x$call$data, envir = environment(f)))
#			
#			, intercept = intercept)
#	else getAllTerms.formula(f, intercept = intercept)
#}

`getAllTerms.hurdle` <- 
`getAllTerms.zeroinfl` <-
function(x, intercept = FALSE, ...) {
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
	
	deps <- termdepmat_combine(lapply(z, attr, "deps"))

	ord <- unlist(lapply(z, attr, "order"))
	n <- sapply(z, length)
	if(length(z) > 1L) ord[-j] <- ord[-(j <- seq_len(n[1L]))] + n[1L]
	zz <- unlist(z)
	Ints <- which(zz == "(Intercept)")
	#zz[Ints] <- "1"
	#zz <- paste0(rep(c("count", "zero")[seq_along(z)], sapply(z, length)),
		#"(", zz, ")")
	zz <- paste0(rep(c("count", "zero")[seq_along(z)], sapply(z, length)),
				 "_", zz)
	
	dimnames(deps) <- list(zz[-Ints], zz[-Ints])
	
	ret <- if(!intercept) zz[-Ints] else zz
	attr(ret, "intercept") <- pmin(Ints, 1)
	attr(ret, "interceptLabel") <- zz[Ints]
	attr(ret, "response") <- attr(z[[1L]], "response")
	attr(ret, "order") <- if(!intercept) order(ord[-Ints]) else ord
	attr(ret, "deps") <- deps
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
	ret <- getAllTerms.terms(terms(x))
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
		ret <- c(ret, as.formula(f[c(1L, length(f))]))
		f <- f[[2L]]
	}
	ret <- lapply(ret, `environment<-`, NULL)
	names(ret) <- sapply(x@estimates@estimates, slot, "short.name")[seq_along(ret)]
	ret <- lapply(ret, getAllTerms.formula, intercept = FALSE)
	
	deps <- termdepmat_combine(lapply(ret, attr, "deps"))

	attrInt <- sapply(ret, attr, "intercept")
	#ret <- unlist(lapply(names(ret), function(i) sprintf("%s(%s)", i, ret[[i]])))
	ret <- unlist(lapply(names(ret), function(i) if(length(ret[[i]]))
						 paste0(i, "(", ret[[i]], ")") else character(0L)))

	dimnames(deps) <- list(ret, ret)

	Ints <- paste0(names(attrInt[attrInt != 0L]), "(Int)")
	if(intercept) ret <- c(Ints, ret)
	attr(ret, "intercept") <- attrInt
	attr(ret, "interceptLabel") <- Ints
	attr(ret, "deps") <- deps
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
	deps <- attr(ret, "deps")
	dn <- gsub("^p\\(", "p(sigma", rownames(deps))
	dimnames(deps) <- list(dn, dn)
	attr(ret, "deps") <- deps
	ret
}

`getAllTerms.MCMCglmm` <- 
function (x, ...) {
	res <- getAllTerms.default(x, ...) 
	attr(res, "random") <- .formulaEnv(.~., environment(formula(x)))
	attr(res, "random.terms") <- deparse(x$Random$formula, control = NULL)[1]
	res
}

`getAllTerms.gamm` <-
function (x, ...) getAllTerms(x$gam, ...)


`getAllTerms.mark` <- 
function (x, intercept = FALSE, ...) {
	
	f <- formula(x, expand = FALSE)[[2L]]
	formlist <- list()
	while(length(f) == 3L && f[[1L]] == "+") {
		formlist <- c(f[[3L]], formlist)
		f <- f[[2L]]
	}
	formlist <- append(f, formlist)
	
	wrapfunc <- function(x, func) if(length(x) == 0L) x else paste0(func, "(", x, ")")

	alltermlist <- lapply(formlist, function(x, intercept) {
		func <- deparse(x[[1L]], control = NULL)
		at <- getAllTerms(terms(eval(call("~", x[[2L]]))), intercept = intercept)
		at[] <- wrapfunc(at, func)
		dn <- wrapfunc(rownames(attr(at, "deps")), func)
		attr(at, "interceptLabel") <- wrapfunc(attr(at, "interceptLabel"), func)
		dimnames(attr(at, "deps")) <- list(dn, dn)
		at
	}, intercept)
	
	retval <- unlist(alltermlist, recursive = TRUE)
	for(a in c("intercept", "interceptLabel")) {
		attr(retval, a) <-	unlist(sapply(alltermlist, attr, a))
	}
	attr(retval, "order") <- order(rep(seq_along(alltermlist), vapply(alltermlist, length, 1L)),
		unlist(lapply(alltermlist, attr, "order")))
	attr(retval, "deps") <- termdepmat_combine(lapply(alltermlist, attr, "deps"))
	retval
}

`getAllTerms.betareg` <-
function(x, intercept = FALSE, ...) {
	f <- formula(x)
	if(length(f[[3L]]) != 1L && f[[3L]][[1L]] == "|"){
		f1 <- call("~", f[[2L]], f[[3L]][[2L]])
		f2 <- call("~", f[[3L]][[3L]])
	} else {
		f1 <- f
		f2 <- NULL
	}
	fs <- lapply(lapply(c(f1, f2), terms.formula, data = model.frame(x)),
		formula)
	z <- lapply(fs, getAllTerms, intercept = TRUE)
	
	deps <- termdepmat_combine(lapply(z, attr, "deps"))
	
	ord <- unlist(lapply(z, attr, "order"))
	n <- sapply(z, length)
	if(length(z) > 1L) ord[-j] <- ord[-(j <- seq_len(n[1L]))] + n[1L]
	zz <- unlist(z)
	Ints <- which(zz == "(Intercept)")

	if(!is.null(f2) && n[2L] != 0L) {
		i.phi <- -seq.int(n[1L])
		zz[i.phi] <- paste("(phi)", zz[i.phi], sep = "_")
	}
	ret <- if(!intercept) zz[-Ints] else zz
	dimnames(deps) <- list(zz[-Ints], zz[-Ints])

	attr(ret, "intercept") <- pmin(Ints, 1)
	attr(ret, "interceptLabel") <- zz[Ints]
	attr(ret, "response") <- attr(z[[1L]], "response")
	attr(ret, "order") <- if(!intercept) order(ord[-Ints]) else ord
	
	attr(ret, "deps") <- deps
	ret
}

`getAllTerms.asreml`  <-
function(x, intercept = FALSE, ...)
getAllTerms.terms(terms(formula(x), ...), intercept = intercept)


`getAllTerms.cpglmm` <-
function (x, intercept = FALSE, ...) 
getAllTerms(x@formula, intercept = intercept)


`getAllTerms` <-
function(x, ...)
UseMethod("getAllTerms")

# TODO: return object of class 'allTerms'
print.allTerms <-
function(x, ...) {
	cat("Model terms: \n")
	if(!length(x)) {
		cat("<None> \n")
	} else {
		print.default(as.vector(x), quote = TRUE)
	}
	ints <- attr(x, "interceptLabel")
	if(!is.null(ints)) {
		cat(ngettext(n = length(ints), "Intercept:", "Intercepts:"), "\n")
		print.default(ints,quote = TRUE)
	}
}
