
coeffs.gamlss <-
function(model) {
    cf <- model[c('mu.coefficients', 'sigma.coefficients', 'nu.coefficients', 'tau.coefficients')]
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
    capture.output(sm <- summary(model))
    ct <- sm[, c(1L, 2L), drop = FALSE]
    .makeCoefTable(ct[, 1L], ct[, 2L], coefNames = names(coeffs(model)))
}

`makeArgs.gamlss` <- 
function(obj, termNames, opt, ...) {
	formulanames <- c("formula", "sigma.formula", "nu.formula", "tau.formula")
	zarg <- umf_terms2formulalist(termNames, opt)
	names(zarg) <- formulanames
	f <- zarg[[1L]][c(1L, NA, 2L)]
	f[[2L]] <- opt$response
	zarg[[1L]] <- f
	zarg
}

getAllTerms.gamlss <-
function(x, intercept = FALSE, ...) {
	formlist <- list(mu = x$mu.formula, sigma = x$sigma.formula, 
        nu = x$nu.formula, tau = x$tau.formula)
	allterms <- lapply(formlist, getAllTerms.formula, intercept = FALSE)
	attrint <- vapply(allterms, attr, 0L, "intercept")
	
    term.prefix <- names(allterms)
	n <- length(allterms)
	rval <- vector("list", n)
	for(i in which(sapply(allterms, length) != 0L))
		rval[[i]] <- paste0(term.prefix[i], "(", allterms[[i]], ")")
	rval <- unlist(rval)
    
	attrint <- vapply(allterms, attr, 0L, "intercept")
	names(attrint) <- term.prefix[match(names(attrint), term.prefix)]

	ints <- paste0(names(attrint[attrint != 0L]), "(",
		unlist(lapply(allterms, "attr", "interceptLabel")), ")")
	ints <- sub("((Intercept))", "(Int)", ints, fixed = TRUE)
    
	deps <- termdepmat_combine(lapply(allterms, attr, "deps"))
	dimnames(deps) <- list(rval, rval)
	if(intercept) rval <- c(ints, rval)
	mode(rval) <- "character"		
	attr(rval, "intercept") <- attrint
	attr(rval, "interceptLabel") <- ints
	attr(rval, "response") <- attr(allterms$mu, "response")
	if(intercept) attr(rval, "interceptIdx") <- seq_along(ints)
	attr(rval, "deps") <- deps
	return(rval)
}
