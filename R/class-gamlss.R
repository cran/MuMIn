
coeffs.gamlss <-
function(model) {
    cf <- model[c('mu.coefficients', 'sigma.coefficients', 'nu.coefficients', 'tau.coefficients')]
    cf <- lapply(cf, function(x) x[!is.na(x)])
	n <- vapply(cf, length, 0L)
    nm <- unlist(lapply(cf, names), recursive = FALSE, use.names = FALSE)
    nm[nm == "(Intercept)"] <- "(Int)"
    rval <- unlist(cf, use.names = FALSE, recursive = FALSE)
    pfx <- c("mu", "sigma", "nu", "tau")
    names(rval) <- paste0(rep(pfx, n), "(", nm, ")")
    rval
}

coefTable.gamlss <-
function(model, ...)  {
	cf <- coeffs(model)
    .makeCoefTable(cf, vcov(model, type = "se"), coefNames = names(cf))
}

`makeArgs.gamlss` <- 
function(obj, termNames, opt, ...) {
	zarg <- umf_terms2formulalist(termNames, opt)
	formulanames <- c(mu = "formula", sigma = "sigma.formula",
		nu = "nu.formula", tau = "tau.formula")[
			attr(opt$allTerms, "term.kind")]
	names(zarg) <- formulanames
	f <- zarg[[1L]][c(1L, NA, 2L)]
	f[[2L]] <- opt$response
	zarg[[1L]] <- f
	zarg
}


getAllTerms.gamlss <-
function(x, intercept = FALSE, offset = TRUE, ...) {
	
	formlist <- list(mu = x$mu.formula, sigma = x$sigma.formula, 
        nu = x$nu.formula, tau = x$tau.formula)
	
	formlist <- formlist[!vapply(formlist, is.null, logical(1L))]
	
	allterms <- lapply(formlist, getAllTerms.formula, intercept = FALSE, offset = offset, ...)
	attrint <- vapply(allterms, attr, 0L, "intercept")
	
    term.prefix <- names(allterms)
	n <- length(allterms)
	rval <- vector("list", n)
	for(i in which(sapply(allterms, length) != 0L)) {
		rval[[i]] <- paste0(term.prefix[i], "(", allterms[[i]], ")")
	}
	rval <- unlist(rval)

	attrint <- vapply(allterms, attr, 0L, "intercept")
	names(attrint) <- term.prefix[match(names(attrint), term.prefix)]

	ints <- paste0(names(attrint[attrint != 0L]), "(",
		unlist(lapply(allterms, "attr", "interceptLabel")), ")")
	ints <- sub("((Intercept))", "(Int)", ints, fixed = TRUE)
    
	depslist <- lapply(allterms, attr, "deps")
	deps <- termdepmat_combine(depslist)
	if(ncol(deps) != 0L)
		colnames(deps) <- rownames(deps) <-
			paste0(rep(term.prefix, sapply(depslist, ncol)),
				"(", colnames(deps), ")")

	#dimnames(deps) <- list(rval, rval)
	if(intercept) rval <- c(ints, rval)
	mode(rval) <- "character"		
	attr(rval, "intercept") <- attrint
	attr(rval, "interceptLabel") <- ints
	attr(rval, "response") <- attr(allterms$mu, "response")
	attr(rval, "term.kind") <- names(formlist)
	if(intercept) attr(rval, "interceptIdx") <- seq_along(ints)
	attr(rval, "deps") <- deps
	return(rval)
}
