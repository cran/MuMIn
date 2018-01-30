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
		
	ord <- lapply(at, attr, "order")
	ordl <- vapply(ord, length, 0, USE.NAMES = FALSE)
	ord <- unlist(ord, use.names = FALSE) + rep(c(0, ordl[-length(ordl)]), ordl)
	
	if(intercept) {
		rval <- c(intLabel, rval)
		ord <- c(seq.int(along.with = intLabel), ord + length(intLabel))
	}
	
	if(hasOffset) attr(rval, "offset") <- offsetTerm
	attr(rval, "random.terms") <- attr(at$cond, "random.terms")
	attr(rval, "random") <- attr(at$cond, "random") 
	attr(rval, "response") <- attr(at$cond, "response") 
	attr(rval, "order") <- ord
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
    fnm <- c("cond", "zi", "disp")
    argnm <- c("formula", "ziformula", "dispformula")
    rval <- umf_terms2formulalist(termNames, opt, replaceInt = "1")[fnm]
	names(rval) <- argnm
    for(i in which(vapply(rval, is.null, FALSE))) rval[[i]] <- ~ 0
	rval$formula <- as.formula(call("~", as.symbol(opt$response), rval$formula[[2L]]), opt$gmFormulaEnv)
	if(inherits(attr(opt$allTerms, "random"), "formula"))
		rval$formula <- update.formula(rval$formula, attr(opt$allTerms, "random"))
	rval
}
