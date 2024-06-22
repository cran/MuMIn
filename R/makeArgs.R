# combination of term names (a character vector), additional list of arbitrary 
# options is accepted. This is much a reverse action to getAllTerms
makeArgs <- function(obj, termNames, opt, ...) UseMethod("makeArgs", obj)

# opt == argsOptions

#argsOptions <- list(
#	response = attr(allTerms0, "response"),
#	intercept = nInts,        ### ONLY .default
#	interceptLabel = interceptLabel,
#	random = attr(allTerms0, "random"),
#	gmCall = gmCall,          ### ONLY .default
#	gmEnv = gmEnv,
#	allTerms = allTerms0,
#	gmCoefNames = gmCoefNames,
#	gmDataHead = if(!is.null(gmCall$data)) {
#		if(eval(call("is.data.frame", gmCall$data), gmEnv))
#			eval(call("head", gmCall$data, 1L), gmEnv) else gmCall$data
#		} else NULL,
#	gmFormulaEnv = gmFormulaEnv
#	)


.getCoefNames <- 
function(formula, data, contrasts, envir = parent.frame()) {
	colnames(eval(call("model.matrix.default",
		object = formula, data = data, contrasts.arg = contrasts), envir = envir))
}

makeArgs.default <- 
function(obj, termNames, opt, ...) {
	reportProblems <- character(0L)
	termNames[termNames %in% opt$interceptLabel] <- "1"
	## XXX: what if length(opt$intercept) > 1 ???
	f <- reformulate(c(if(!opt$intercept) "0" else if (!length(termNames)) "1", termNames), response = opt$response)

	environment(f) <- opt$gmFormulaEnv
	ret <- list(formula = f)
	if(!is.null(opt$gmCall$start)) {
		coefNames <- fixCoefNames(.getCoefNames(f, opt$gmDataHead,
			opt$gmCall$contrasts, envir = opt$gmEnv))
		idx <- match(coefNames, opt$gmCoefNames)
		if(anyNA(idx)) reportProblems <-
			append(reportProblems, "cannot subset 'start' argument. Coefficients in the model do not exist in 'global.model'")
		else ret$start <- substitute(start[idx], list(start = opt$gmCall$start,
			idx = idx))
	}
	#attr(ret, "formulaList") <- list(f)
	attr(ret, "problems") <- reportProblems
	ret
}

makeArgs.gls <- 
makeArgs.wgee <- 
function(obj, termNames, opt, ...) {
	ret <- makeArgs.default(obj, termNames, opt)
	names(ret)[1L] <- "model"
	ret
}

`makeArgs.asreml` <- 
makeArgs.MCMCglmm <-
makeArgs.lme <- 
function(obj, termNames, opt, ...) {
	ret <- makeArgs.default(obj, termNames, opt)
	names(ret)[1L] <- "fixed"
	ret
}


`makeArgs.glmmadmb` <- 
`makeArgs.merMod` <-    ## since lme4-0.99999911-0
`makeArgs.mer` <- 
function(obj, termNames, opt, ...) {
	ret <- makeArgs.default(obj, termNames, opt)
	if(!is.null(opt$random)) ret[['formula']] <-
		update.formula(ret[['formula']], opt$random)
	ret
}

# clmm needs explicit "1" if no other FX terms
`makeArgs.clmm` <- 		## Class 'clmm'  from package 'ordinal':
function(obj, termNames, opt, ...) {
	ret <- makeArgs.merMod(obj, termNames, opt, ...)
	if(length(termNames) == 1L && identical(termNames[1L], opt$interceptLabel))
		ret$formula[[3L]] <- call("+", 1, ret$formula[[3L]])
	ret
}

`makeArgs.coxph` <- 
function(obj, termNames, opt, ...) {
	ret <- makeArgs.default(obj, termNames, opt)
	ret$formula <- update.formula(ret$formula, . ~ . + 1)
	ret
}

`makeArgs.betareg` <- 
function(obj, termNames, opt, ...) {
	i <- termNames %in% opt$interceptLabel
	termNames[i] <- gsub("(Intercept)", "1", termNames[i], fixed = TRUE)
	j <- grepl("^\\(phi\\)_", termNames)
	# TODO: zero-length terms in reformulate
	zarg <- list(beta = formula(terms.formula(reformulate(termNames[!j]), simplify = TRUE)))
	if(any(j))
		zarg$phi <- formula(terms.formula(reformulate(substring(termNames[j], 7L)), simplify = TRUE))
	
	zarg <- lapply(zarg, `environment<-`, opt$gmFormulaEnv)
	fexpl <- zarg$beta[[2L]]
	if(!is.null(zarg$phi)) fexpl <- call("|", fexpl, zarg$phi[[2L]])
		else zarg$phi <- NULL
	ret <- list(formula = call("~", opt$response, fexpl))
	#attr(ret, "formulaList") <- zarg
	ret
}

`makeArgs.hurdle` <- 
`makeArgs.zeroinfl` <-
function(obj, termNames, opt, ...) {

	intType <- substring(opt$interceptLabel, 0,
		regexpr("_", opt$interceptLabel, fixed = TRUE) - 1L)

	i <- termNames %in% opt$interceptLabel
	termNames[i] <- gsub("(Intercept)", "1", termNames[i], fixed = TRUE)
	pos <- regexpr("_", termNames, fixed = TRUE)

	fnames <- c("count", "zero")
	
	# TODO: zero-length terms in reformulate
	zarg <- split(substring(termNames, pos + 1L, 256L),
		substring(termNames, 1L, pos - 1L))
	for(j in fnames) zarg[[j]] <-
		if(is.null(zarg[[j]])) {
			if(j %in% intType) ~1 else ~0
		} else formula(terms.formula(reformulate(as.character(zarg[[j]]),
			intercept = j %in% intType),
			simplify = TRUE))

	zarg <- lapply(zarg, `environment<-`, opt$gmFormulaEnv)
	zarg <- zarg[fnames]
	fexpl <- zarg$count[[2L]]
	if(!is.null(zarg$zero)) fexpl <-
		call("|", fexpl, zarg$zero[[2L]]) else zarg$zero <- NULL
	ret <- list(formula = call("~", opt$response, fexpl))
	#attr(ret, "formulaList") <- zarg
	ret
}

`makeArgs.coxme` <-
`makeArgs.lmekin` <-
function(obj, termNames, opt, ...) {
	ret <- makeArgs.default(obj, termNames, opt)
	ret$formula <- update.formula(update.formula(ret$formula, . ~ . + 1),
		opt$random)
	ret
}

`makeArgs.mark` <- 
function(obj, termNames, opt, ...) {
	
	interceptLabel <- "(Intercept)"
	termNames <- sub(interceptLabel, "1", termNames, fixed = TRUE)
	rxres <- regexpr("^([a-zA-Z]+)\\((.*)\\)$", termNames, perl = TRUE)
	cs <- attr(rxres, "capture.start")
	cl <- attr(rxres, "capture.length")
	parname <- substring(termNames, cs[, 1L], cs[, 1L] + cl[,1L] - 1L)
	parval <- substring(termNames, cs[, 2L], cs[, 2L] + cl[,2L] - 1L)

	formulaList <- lapply(split(parval, parname), function(x) {
		int <- x == "1"
		x <- x[!int]
		res <- if(!length(x))
				if(int) ~ 1 else ~ 0 else 
			reformulate(x, intercept = any(int))
		environment(res) <- opt$gmFormulaEnv
		res
	})
	
	mpar <- if(is.null(obj$model.parameters))
		eval(opt$gmCall$model$parameters) else
		obj$model.parameters
	for(i in names(mpar)) mpar[[i]]$formula <- formulaList[[i]]
	#ret <- list(model.parameters = mpar)
	
	if(opt$gmCall[[1L]] == "run.mark.model") {
		arg.model <- opt$gmCall$model
		arg.model$parameters <- mpar
		ret <- list(model = arg.model)
	} else {
		ret <- list(model.parameters = mpar)
	}
	#attr(ret, "formulaList") <- formulaList
	ret
}


`makeArgs.aodml` <-
function(obj, termNames, opt, ...) {
	if(sys.nframe() > 2L && (parent.call <- sys.call(-2L))[[1L]] == "dredge" &&
	   !is.null(get_call(obj)$fixpar))
		stop(simpleError("'aodml' models with constant parameters cannot be handled by 'dredge'",
						 call = parent.call))
	makeArgs.default(obj, termNames, opt, ...)
}
