
`formula.unmarkedFit` <- function (x, ...) x@formula

umf_terms2formulalist <- 
function(termNames, opt, replaceInt = "(1)") {
	i <- termNames %in% opt$interceptLabel
	termNames[i] <- gsub("(Int)", replaceInt, termNames[i], fixed = TRUE)
	fexpr <- lapply(termNames, str2lang)
	
	nm  <- as.character(lapply(fexpr, "[[", 1L))

	fsplt <- split(sapply(fexpr, "[[", 2L), nm)[nm[!duplicated(nm)]]
	farg <- lapply(fsplt, function(z) {
		if(! 1 %in% z) z <- c(0, z)
		rval <- z[[1L]]
		n <- length(z)
		if(n > 1) for(i in 2L:n) rval <- call("+", rval, z[[i]])
		as.formula(call("~", rval),  opt$gmFormulaEnv)
	})
	farg[] <- lapply(farg, `environment<-`, opt$gmFormulaEnv)
	farg
}

umf.formlist <-
function(x, specs = getspecs(x, TRUE), names.item = "model") {

    specs <- specs[!is.na(specs$formula.slot), ]
    
    fl <- if(.hasSlot(x, "formlist")) {
        fl <- x@formlist
        fl[match(specs$formlist.name, names(fl))]
    } else if(.hasSlot(x, "formula") &&
              all(startsWith(specs$formula.slot, "formula"))) {
        fl <- frmsplit(x@formula)
        fl[match(specs$formula.slot, paste0("formula", seq_along(fl)))]
    } else
        lapply(specs$formula.slot, function(name) slot(x, name))
    names(fl) <- specs[[names.item]]
    fl
}

getAllTerms.unmarkedFit <-
function(x, intercept = FALSE, ...) {
    
    formlist <- umf.formlist(x, names.item = "short.name")
    term.prefix <- names(formlist) 
    

	allterms <- lapply(formlist, getAllTerms.formula, intercept = FALSE)
	
	.add.prefix <- function(pfx, s) paste0(pfx, "(", s, ")")
	
	allterms2 <-
	.mapply(function(ato, pfx, addInt) {
		if(length(ato) != 0L)
			ato[] <- .add.prefix(pfx, c(ato))
		for(attrname in c("interceptLabel", "offset", "random.terms"))
			if(!is.null(attr(ato, attrname)))
				attr(ato, attrname) <- .add.prefix(pfx, attr(ato, attrname))
		if(nrow(attr(ato,"deps")) != 0L)
			dimnames(attr(ato,"deps")) <-
				rep(list(setdiff(ato[], attr(ato,"offset"))),  2L)
		if(!is.null(intLab <- attr(ato, "interceptLabel"))) {
			attr(ato, "interceptLabel") <- intLab <- 
				sub("((Intercept))", "(Int)", intLab, fixed = TRUE)
			if(addInt)
				ato[seq_len(length(ato) + length(intLab))] <-
					append(ato, intLab, after = 0L)
		}
		ato
	}, list(ato = allterms, pfx = term.prefix),
		MoreArgs = list(addInt = isTRUE(intercept)))
	
	n <- length(allterms)
	rval <- unlist(allterms2)
	mode(rval) <- "character" # in case of zero-length

	attrint <- vapply(allterms2, attr, integer(1L), "intercept")
	names(attrint) <- term.prefix
	attr(rval, "intercept") <- attrint
	
	for(attrname in c("interceptLabel", "random.terms"))
		attr(rval, attrname) <- unlist(lapply(allterms2, attr, attrname))
	
	if(!is.null(attr(rval, "random.terms"))) {
		ranef.forms <- lapply(allterms2, attr, "random")
		names(ranef.forms) <- term.prefix
		# remove lhs from formulas:
		attr(rval, "random") <-
			lapply(ranef.forms[!vapply(ranef.forms, is.null,
				logical(1L))], "[", -2L)
	}
	
	if(intercept)
		attr(rval, "interceptIdx") <-
			which(unlist(sapply(allterms2, function(ato) {
		   rval <- logical(length(ato))
		   rval[attr(ato, "intercept")] <- TRUE
		   rval
		}, simplify = FALSE)))
	attr(rval, "deps") <- termdepmat_combine(lapply(allterms2, attr, "deps"))
	return(rval)
}


`makeArgs.unmarkedFit` <- 
function(obj, termNames, opt, ...) {
	
	specs <- getspecs(obj)
	formulanames <- as.character(specs$formula.arg)
	single_formula <- all(startsWith(formulanames, "formula"))
	
	# NOTE: elements are named after full short.name
	zarg <- umf_terms2formulalist(termNames, opt)
	
	if(!is.null(opt$random)) {
		for(a in names(opt$random))
			zarg[[a]] <- update.formula(zarg[[a]], opt$random[[a]])
	}

	zarg <- zarg[specs$short.name]
	names(zarg) <- specs$formula.arg
	
	if(single_formula) {
		n <- length(zarg)
		zarg <- zarg[paste0("formula", seq_len(n))]
    	form <- zarg[[1L]]
		if(n > 1L) for(i in seq.int(2L, n)) form <- call("~", form, zarg[[i]][[2L]])
	    form <- as.formula(form, env = environment(zarg[[1L]]))
		environment(form) <- environment(zarg[[1L]])
		list(formula = form)
	} else {
		names(zarg) <- formulanames
		zarg
	}
}
