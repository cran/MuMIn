umf_formlist <-
function(x) {
	if("formlist" %in% slotNames(x)) {
		x@formlist
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

umf_standardize_estnames <-
function(x) {	
	tr <- c(lambda = "lam", gamma = "gam", omega = "om", det = "p",
				col = "gam", ext = "eps", epsilon = "eps")
	i <- match(x, names(tr), nomatch = 0L)
	if(all(i == 0L)) return(x)
	x[i != 0L] <- tr[i]
	x
}

umf_get_specs <- 
function(x) {
	spc <- umf_specs[[class(x)[1L]]]
	if(is.null(spc)) 
	  stop(gettextf("this unmarkedFit subclass or structure is unknown to MuMIn",
	  class(x)[1L]))

	estsn <- sapply(x@estimates@estimates, "slot", "short.name")
	estsn <- estsn[!duplicated(estsn)] #  for "unmarkedFitDS" state="lam",det="p",scale="p"
	estsnpfx <- umf_shortname2estname(estsn)
	
	i <-
		sapply(lapply(spc, "[", 2L, ), function(x) setequal(x[x != ""], estsnpfx)) &
		sapply(lapply(spc, "[", 3L, ), function(x) setequal(x[x != ""], names(estsn)))
		
	if(!any(i)) stop(gettextf("unknown \"%s\" model structure", class(x)[1L]))
	rval <- spc[[which(i)[1L]]]
	attr(rval, "estsnpfx") <- estsnpfx
	attr(rval, "estsn") <- estsn
	rval
}


umf_terms2formulalist <- 
function(termNames, opt) {
	i <- termNames %in% opt$interceptLabel
	termNames[i] <- gsub("(Int)", "(1)", termNames[i], fixed = TRUE)
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

getAllTerms.unmarkedFit <-
function(x, intercept = FALSE, ...) {

	spc <- umf_get_specs(x)
	formlist <- umf_formlist(x)
	formnames <- spc[1L, ]

	allterms <- lapply(formlist, getAllTerms.formula, intercept = FALSE)
	
	fnames <- spc[2L, ]
	i <- match(fnames, attr(spc, "estsnpfx"), nomatch = 0L)
    fnames[i != 0] <- attr(spc, "estsn")[i]
	
	if(is(x, "unmarkedFitDS")) {
		if(!is.na(fnames['p'])) {
			z <- allterms[[ i <- which(names(fnames) == 'p') ]]
			detprefix <- switch(x@keyfun, uniform = "", halfnorm = "sigma",
				hazard = "shape", exp = "rate", "")
			z[] <- paste0(detprefix, z)
			attr(z,"interceptLabel") <- paste0(detprefix, attr(z,"interceptLabel"))
			allterms[[i]] <- z
		}
		if(any(i <- fnames == "")) fnames[i] <- paste0("dummy", 1L:sum(i))
	}
	if(any(i <- fnames == "")) fnames[i] <- colnames(spc)[i]
	
	n <- length(allterms)
	rval <- vector("list", n)
	for(i in seq.int(n))
		rval[[i]] <- if(length(allterms[[i]])) paste0(fnames[i], "(",
			allterms[[i]], ")") else character(0L)
	rval <- unlist(rval)
	
	deps <- termdepmat_combine(lapply(allterms, attr, "deps"))
	dimnames(deps) <- list(rval, rval)
	
	attrint <- sapply(allterms, attr, "intercept")
	names(attrint) <- if(is.null(names(attrint)))
	    unname(fnames[formnames != ""]) else
		fnames[match(names(attrint), formnames)]

	ints <- paste0(names(attrint[attrint != 0L]), "(",
		unlist(lapply(allterms, "attr", "interceptLabel")), ")")
	ints <- sub("((Intercept))", "(Int)", ints, fixed = TRUE)
	
	if(intercept) rval <- c(ints, rval)
	attr(rval, "intercept") <- attrint
	attr(rval, "interceptLabel") <- ints
	if(intercept) attr(rval, "interceptIdx") <- seq_along(ints)
	attr(rval, "deps") <- deps
	return(rval)
}


`makeArgs.unmarkedFit` <- 
function(obj, termNames, opt, ...) {
	spc <- umf_get_specs(obj)
	
	 # TODO: set attr(, "argsOrder") <- k
	if(isTRUE(attr(spc, "revArgs"))) {
		k <- which(spc[1L, ] != "")
		spc[, k] <- spc[, rev(k), drop = FALSE]
		colnames(spc)[k] <- colnames(spc)[rev(k)]
	}
		
	formulanames <- spc[1L, spc[1L, ] != ""]
	single_formula <- all(formulanames == "formula")
	
	# NOTE: elements are named after full short.name
	zarg <- umf_terms2formulalist(termNames, opt)
	names(zarg) <- umf_standardize_estnames(umf_shortname2estname(names(zarg)))

	#stopifnot(all(names(zarg) %in% colnames(spc)[spc[1, ] != ""])) # DEBUG
	
    zarg <- zarg[match(colnames(spc)[spc[1, ] != ""], names(zarg))]
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


###

.umf_compute_specs <- 
function(x) {
	fnsym <- get_call(x)[[1L]]
	argnm <- names(formals(eval(fnsym)))
	formargnm <- argnm[grep("[a-z]*formula$", argnm)]
	estshnm <- sapply(x@estimates@estimates, "slot", "short.name")
	estshpfx <- umf_shortname2estname(estshnm)
	nform <- sum(all.names(formula(x)) == "~")
	names(formargnm) <- sub("formula$", "", formargnm, perl = TRUE)
	
	translnm <- umf_standardize_estnames(c(names(formargnm), estshpfx))
	
	u <- unique(translnm)
	u <- u[u != ""]
	
	m <- matrix("", ncol = length(u), nrow = 3L)
	dimnames(m) <- list(c("argument", "est:short.name", "est:label"), u)
	
	if(all(formargnm == "formula")) {
		m[1L, 1L:min(nform, ncol(m))] <- "formula"
	} else
		m[1L, match(translnm[seq.int(length(formargnm))], u)] <- formargnm
	
	m[2, j <- match(umf_standardize_estnames(estshpfx), u)] <- estshpfx
	m[3, j] <- names(estshnm)
	
	attr(m, "umf_class") <- class(x)
	attr(m, "n_formulas") <- nform
	return(m)
}
