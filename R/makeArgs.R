# 'makeArgs' internal generic to generate arguments from combination of
# term names (a character vector), additional list of arbitrary options is accepted.
# This is much a reverse action to getAllTerms
makeArgs <- function(obj, termNames, comb, opt, ...) UseMethod("makeArgs", obj)

.getCoefNames <- function(formula, data, contrasts, envir = parent.frame()) {
	colnames(eval(call("model.matrix.default",
		object = formula, data = data, contrasts.arg = contrasts), envir = envir))
}

makeArgs.default <- function(obj, termNames, comb, opt, ...) {
	reportProblems <- c()
	termNames0 <- termNames

	termNames[termNames %in% opt$interceptLabel] <- "1"
	f <- reformulate(c(if(!opt$intercept) "0", termNames), response=opt$response)
	environment(f) <- opt$gmFormulaEnv
	ret <- list(formula = f)
	if(!is.null(opt$gmCall$start)) {
		coefNames <- fixCoefNames(.getCoefNames(f, opt$gmDataHead, opt$gmCall$contrasts, envir = opt$gmEnv))

		idx <- match(coefNames, opt$gmCoefNames)

		if(any(is.na(idx))) {
			reportProblems <- append(reportProblems, "cannot subset 'start' argument. Coefficients in generated model do not exist in the global model")
		} else {
			ret$start <- substitute(start[idx], list(start = opt$gmCall$start, idx = idx))
		}
	}
	attr(ret, "formulaList") <- list(f)
	attr(ret, "problems") <- reportProblems
	ret
}

makeArgs.gls <- function(obj, termNames, comb, opt, ...) {
	ret <- makeArgs.default(obj, termNames, comb, opt)
	names(ret)[1] <- "model"
	ret
}
makeArgs.lme <- function(obj, termNames, comb, opt, ...) {
	ret <- makeArgs.default(obj, termNames, comb, opt)
	names(ret)[1] <- "fixed"
	ret
}
makeArgs.mer <- function(obj, termNames, comb, opt, ...) {
	ret <- makeArgs.default(obj, termNames, comb, opt)
	ret$formula <- update.formula(ret$formula, opt$random)
	ret
}

# used by makeArgs.unmarkedFit*
`.makeUnmarkedFitFunnyFormulas` <- function(termNames, opt, fnames) {
	i <- termNames %in% opt$interceptLabel
	termNames[i] <- gsub("(Int)", "(1)", termNames[i], fixed = TRUE)
	zexpr <- lapply(termNames, function(x) parse(text = x)[[1L]])
	zsplt <- split(sapply(zexpr, "[[", 2L), as.character(sapply(zexpr, "[[", 1L)))
	zarg <- lapply(zsplt, function(x) reformulate(as.character(x)))
	zarg <- zarg[fnames]
	names(zarg) <- fnames
	i <- sapply(zarg, is.null)
	zarg[i] <- rep(list(~ 1), sum(i))
	zarg <- lapply(zarg, `environment<-`, opt$gmFormulaEnv)
	#print(zarg)
	zarg
}

# This works for all the child classes provided that the arguments containing
# formulas are named *formula, and they are at the beginnning of the call, and
# if there is no more than 6 of them, and if there is a single 'formula'
# argument, it is double-one-sided.
# Guessing the argument names is only negligibly slower than when they are
# provided to the function, so the methods for child classes are commented out.

`makeArgs.unmarkedFit` <- function(obj, termNames, comb, opt,
	fNames = sapply(obj@estimates@estimates, slot, "short.name"),
	formula.arg = "formula",
	...) {
	zarg <- .makeUnmarkedFitFunnyFormulas(termNames, opt, fNames)
	if(missing(formula.arg)) { # fallback
		nm <- names(opt$gmCall)[-1]
		formula.arg <- nm[grep(".*formula$", nm[1:7])]
		#print(formula.arg)
	}
	if(length(formula.arg) == 1L && formula.arg == "formula") {
		i <- sapply(zarg, is.null)
		zarg[i] <- rep(list(~ 1), sum(i))
		ret <- list(formula = call("~", zarg[[2]], zarg[[1]][[2]]))
		attr(ret, "formulaList") <- zarg
		ret
	} else {
		names(zarg) <- formula.arg
		zarg
	}
}

`makeArgs.unmarkedFitDS` <-
function(obj, termNames, comb, opt, ...)  {
	termNames <- sub("^p\\(sigma(.+)\\)", "p(\\1)", termNames, perl = TRUE)
	termNames[termNames == "p((Intercept))"] <- "p(1)"
	makeArgs.unmarkedFit(obj, termNames, comb, opt, c("lam", "p"),
		"formula")
}


#`makeArgs.unmarkedFitColExt` <- function(obj, termNames, comb, opt, ...)
#	makeArgs.unmarkedFit(obj, termNames, comb, opt, c("psi", "col", "ext", "p"),
#		c("psiformula", "gammaformula", "epsilonformula", "pformula"))
#
#
#`makeArgs.unmarkedFitGMM` <- function(obj, termNames, comb, opt, ...)
#	makeArgs.unmarkedFit(obj, termNames, comb, opt,
#		c("lambda", "phi", "p"),
#		c("lambdaformula", "phiformula", "pformula") )
#
#`makeArgs.unmarkedFitPCO` <- function(obj, termNames, comb, opt, ...)
#	makeArgs.unmarkedFit(obj, termNames, comb, opt,
#		c("lam", "gamConst", "omega", "p"),
#		c("lambdaformula", "gammaformula", "omegaformula", "pformula"))
#
#`makeArgs.unmarkedFitOccu` <- function(obj, termNames, comb, opt, ...)
#	makeArgs.unmarkedFit(obj, termNames, comb, opt, c("psi", "p"),
#		"formula")
