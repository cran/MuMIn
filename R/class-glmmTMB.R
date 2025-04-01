
`getAllTerms.glmmTMB` <-
function(x, intercept = FALSE, offset = TRUE, ...) {
	
    at <- x$modelInfo$allForm[c("formula", "ziformula", "dispformula")]
    at <- lapply(at, getAllTerms.formula, intercept = FALSE, offset = TRUE)
    names(at) <- c("cond", "zi", "disp")
    deps <- termdepmat_combine(lapply(at, attr, "deps"))
    attrInt <- sapply(at, attr, "intercept")
	
	rval <- unlist(lapply(names(at), function(i)
		if (length(at[[i]])) paste0(i, "(", at[[i]], ")")
		else character(0L)))

	off <- lapply(at, function(tt) if(is.null(off <- attr(tt,"offset"))) integer(0L) else match(off, tt))
	loff <- sapply(off, length)
	off <- if(any(loff > 0L)) {
		unlist(off) + rep(c(0L, cumsum(sapply(at[-length(at)], length))), sapply(off, length))
	} else NULL
	
	if(hasOffset <- !is.null(off)) {
		offsetTerm <- rval[off]
		if(!isTRUE(offset)) {
			depnames <- rval <- rval[-off]
		} else depnames <- rval[-off]
	} else depnames <- rval
	
	dimnames(deps) <- list(depnames, depnames)
    intLabel <- paste0(names(attrInt[attrInt != 0L]), "((Int))")

	sortorder <- lapply(at, attr, "sortorder")
	sortorderl <- vapply(sortorder, length, 0, USE.NAMES = FALSE)
	sortorder <- unlist(sortorder, use.names = FALSE) + 
        rep(c(0, sortorderl[-length(sortorderl)]), sortorderl)
	
	if(intercept) {
		rval <- c(intLabel, rval)
		sortorder <- c(seq.int(along.with = intLabel), sortorder + length(intLabel))
	}
	
	if(hasOffset) attr(rval, "offset") <- offsetTerm
		
	rt <- lapply(at, "attr", "random.terms")
    if(!all(vapply(rt, is.null, FALSE))) {
		rt <- paste0(rep(names(rt), vapply(rt, length, 0L)), "(", unlist(rt), ")")
		random <- reformulate(c(".", rt), response = ".")
		environment(random) <- environment(x$modelInfo$allForm$combForm)
	} else
		rt <- random <- NULL

	attr(rval, "random.terms") <- rt
	attr(rval, "random") <- random
	attr(rval, "response") <- attr(at$cond, "response") 
	attr(rval, "sortorder") <- sortorder
    attr(rval, "intercept") <- attrInt
    attr(rval, "interceptLabel") <- intLabel
    attr(rval, "deps") <- deps
    return(rval)
}


coefTable.glmmTMB <-
function (model, ...) {
    dfs <- df.residual(model)
    cf <- summary(model, ...)$coefficients
    cf1 <- do.call("rbind", cf)
    nm <- paste0(rep(names(cf), sapply(cf, NROW)), "(", rownames(cf1), ")")
    nm <- sub("\\(\\(Intercept\\)\\)$", "((Int))", nm)
    .makeCoefTable(cf1[, 1L], cf1[, 2L], dfs, coefNames = nm)
}


coeffs.glmmTMB <-
function(model) {
    coefTable(model)[, 1L]    
}

`makeArgs.glmmTMB` <- 
function(obj, termNames, opt, ...) {
	
	.addRanTermToFormula <- function(f, r) {
		if(is.null(r)) return(f)
		dot <- as.symbol(".")
		rflhs <- call("+", dot, call("(", r))
		if(is.null(f)) f <- ~ 1
		update.formula(f, as.formula(if(length(f) == 2L) call("~", rflhs) else call("~", dot, rflhs)))
	}
	
	fnm <- c("cond", "zi", "disp")

	randomterms <- lapply(attr(opt$allTerms, "random.terms"), str2lang)
	names(randomterms) <- vapply(lapply(randomterms, "[[", 1L), as.character, "")
	randomterms <- lapply(randomterms, "[[", 2L)
	rval <- umf_terms2formulalist(termNames, opt, replaceInt = "1")[fnm]
	for(i in fnm)
	   while(i %in% names(randomterms)) {
			rval[[i]] <- .addRanTermToFormula(rval[[i]], randomterms[[i]])
			randomterms[[i]] <- NULL
	   }
	
    argnm <- c("formula", "ziformula", "dispformula")
	
	names(rval) <- argnm
    for(i in which(vapply(rval, is.null, FALSE))) {
		rval[[i]] <- ~ 0
		environment(rval[[i]]) <- opt$gmFormulaEnv
	}
	
	# XXX: Why it was `as.symbol(opt$response)` ?
	rval$formula <- as.formula(call("~", opt$response, rval$formula[[2L]]), opt$gmFormulaEnv)
		
	rval
}
