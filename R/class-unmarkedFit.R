umf_formlist <-
function(x) {
	sn <- slotNames(x)
	if("formlist" %in% sn) {
		x@formlist
	} else if(length(j <- grep("^.+formula$", sn))) {
		sn <- sn[j]
		sapply(sn, slot, object = x, simplify = FALSE)
	} else {
		env <- environment(formula(x))
		f <- formula(x)
		rval <- list()
		while(is.call(f) && f[[1L]] == "~") {
			rval <- c(as.formula(f[c(1L, length(f))], env = env), rval)
			f <- f[[2L]]
		}
		lapply(rval, `environment<-`, env)
	}
}

umf_shortname2estname <-
function(x)
sub("[[:upper:]].+$", "", x)


umf_get_specs <- 
function(x) {
	specs <- umf_specs[[class(x)[1L]]]
	
	if(is.null(specs)) 
	  stop(gettextf("this 'unmarkedFit' subclass or structure is unknown to MuMIn",
	  class(x)[1L]))

	estsn <- sapply(x@estimates@estimates, "slot", "short.name")
	estsn <- estsn[!duplicated(estsn)] #  for "unmarkedFitDS" state="lam",det="p",scale="p"
	estsnpfx <- umf_shortname2estname(estsn)
	
	i <-
		vapply(lapply(specs, "[", 2L, ), function(x) setequal(x[!is.na(x)], estsnpfx), logical(1L)) &
		vapply(lapply(specs, "[", 3L, ), function(x) setequal(x[!is.na(x)], names(estsn)), logical(1L))
		
	if(!any(i)) stop(gettextf("unknown \"%s\" model structure", class(x)[1L]))
	rval <- specs[[which(i)[1L]]]
	attr(rval, "estsnpfx") <- estsnpfx
	attr(rval, "estsn") <- estsn
	rval
}

umf_terms2formulalist <- 
function(termNames, opt, replaceInt = "(1)") {
	i <- termNames %in% opt$interceptLabel
	termNames[i] <- gsub("(Int)", replaceInt, termNames[i], fixed = TRUE)
	fexpr <- lapply(termNames, function(x) parse(text = x)[[1L]])
	
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

umf_ds_fix_names <-
function(x, fm) {
	if(is(fm, "unmarkedFitDS")) {
		x <- sub(paste0("p(", switch(fm@keyfun, uniform = "", halfnorm = "sigma",
				hazard = "shape", exp = "rate", "")), "p(", x, fixed = TRUE)
		x <- sub("\\(\\(([^\\)]+)\\)\\)", "(\\1)", x, perl = TRUE)
		x <- sub("(Intercept)", "(Int)", x, fixed = TRUE)
	}
	x
}

# TODO: optimize
getAllTerms.unmarkedFit <-
function(x, intercept = FALSE, ...) {

	specs <- umf_get_specs(x)
	shortnames <- attr(specs, "estsn")
	formlist <- umf_formlist(x)
	allterms <- lapply(formlist, getAllTerms.formula, intercept = FALSE)
	attrint <- vapply(allterms, attr, 0L, "intercept")
	specs <- specs[, !is.na(specs["formula.arg", ]), drop = FALSE]
	if(!is.null(names(formlist)))
		specs <- specs[,
			match(specs["formulaItemName", ], names(formlist), nomatch = 0L),
				drop = FALSE]
	sn1 <- names(specs["estimate:itemName", ])
	names(sn1) <- specs["estimate:itemName", ]
	j <- match(names(shortnames), names(sn1), nomatch = 0L)
	sn1[j] <- shortnames[j != 0L]
	shortnames <- sn1
	term.prefix <- specs["estimate:short.name", ]
	term.prefix[] <- shortnames
	term.prefix[i] <- names(term.prefix)[i <- is.na(term.prefix)]
	
	if(is(x, "unmarkedFitDS")) {
		if(!anyNA(term.prefix['p'])) {
			z <- allterms[[ i <- which(names(term.prefix) == "p") ]]
			detprefix <- switch(x@keyfun, uniform = "", halfnorm = "sigma",
				hazard = "shape", exp = "rate", "")
			z[] <- paste0(detprefix, z)
			attr(z,"interceptLabel") <- paste0(detprefix, attr(z,"interceptLabel"))
			allterms[[i]] <- z
		}
		if(any(i <- is.na(term.prefix))) term.prefix[i] <- paste0("dummy", 1L:sum(i))
	}
	
	if(any(i <- is.na(term.prefix))) term.prefix[i] <- colnames(specs)[i]
	n <- length(allterms)
	rval <- vector("list", n)
	for(i in which(sapply(allterms, length) != 0L))
		rval[[i]] <- paste0(term.prefix[i], "(", allterms[[i]], ")")
	rval <- unlist(rval)
	
	attrint <- vapply(allterms, attr, 0L, "intercept")
	names(attrint) <- if(is.null(names(attrint)))
	    unname(term.prefix) else
		term.prefix[match(names(attrint), specs["formulaItemName", ])]

	ints <- paste0(names(attrint[attrint != 0L]), "(",
		unlist(lapply(allterms, "attr", "interceptLabel")), ")")
	ints <- sub("((Intercept))", "(Int)", ints, fixed = TRUE)

	deps <- termdepmat_combine(lapply(allterms, attr, "deps"))
	dimnames(deps) <- list(rval, rval)
	
	if(intercept) rval <- c(ints, rval)
	mode(rval) <- "character"		

	attr(rval, "intercept") <- attrint
	attr(rval, "interceptLabel") <- ints
	if(intercept) attr(rval, "interceptIdx") <- seq_along(ints)
	attr(rval, "deps") <- deps
	return(rval)

}

`makeArgs.unmarkedFit` <- 
function(obj, termNames, opt, ...) {
	
	specs <- umf_get_specs(obj)
	formulanames <- specs[1L, !is.na(specs[1L, ])]
	single_formula <- all(formulanames == "formula")
	
	# NOTE: elements are named after full short.name
	zarg <- umf_terms2formulalist(termNames, opt)
	names(zarg) <- 	umf_shortname2estname(names(zarg))
	#zarg <- zarg[specs["estimate:short.name", !is.na(specs["formula.arg", ])]]
	zarg <- zarg[colnames(specs)[!is.na(specs["formula.arg", ])]]
	
	if(single_formula) {
		n <- length(zarg)
    	form <- zarg[[1L]]
		if(n > 1L) for(i in 2L:n) form <- call("~", form, zarg[[i]][[2L]])
	    form <- as.formula(form, env = environment(zarg[[1L]]))
		#XXX ? environment(form) <- environment(zarg[[1L]])
		list(formula = form)
	} else  {
		names(zarg) <- formulanames
		zarg
	}
}

`makeArgs.unmarkedFitDS` <-
function(obj, termNames, opt, ...)  {
	termNames[] <- umf_ds_fix_names(termNames, obj)
	opt$interceptLabel <- umf_ds_fix_names(opt$interceptLabel, obj)
	makeArgs.unmarkedFit(obj, termNames, opt)
}


`formula.unmarkedFit` <- function (x, ...) x@formula
