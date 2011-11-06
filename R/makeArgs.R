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
	#DebugPrint(ret)
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
.makeUnmarkedFitFunnyFormulas <- function(termNames, opt, fnames) {
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
	zarg
}

makeArgs.unmarkedFitColExt <- function(obj, termNames, comb, opt, ...)  {
	zarg <- .makeUnmarkedFitFunnyFormulas(termNames, opt, c("psi", "col", "ext", "p"))
	#names(formals(colext)[seq_along(zarg)])
	names(zarg) <- c("psiformula", "gammaformula", "epsilonformula", "pformula")
	zarg
}

makeArgs.unmarkedFit  <- function(obj, termNames, comb, opt, ...)  {
	fNames <- sapply(obj@estimates@estimates, slot, "short.name")
	zarg <- .makeUnmarkedFitFunnyFormulas(termNames, opt, fNames)
	i <- sapply(zarg, is.null)
	zarg[i] <- rep(list(~ 1), sum(i))
	ret <- list(formula = call("~", zarg[[2]], zarg[[1]][[2]]))
	attr(ret, "formulaList") <- zarg
	ret
}

makeArgs.unmarkedFitOccu <- function(obj, termNames, comb, opt, ...)  {
	zarg <- .makeUnmarkedFitFunnyFormulas(termNames, opt, c("psi", "p"))
	i <- sapply(zarg, is.null)
	zarg[i] <- rep(list(~ 1), sum(i))
	ret <- list(formula = call("~", zarg$p, zarg$psi[[2]]))
	attr(ret, "formulaList") <- zarg
	ret
}
