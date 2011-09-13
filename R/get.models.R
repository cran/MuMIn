`get.models` <-
function(dd, subset = delta <= 4, ...) {

	#subset <- if (missing(subset)) quote() else substitute(subset)
	subset <- eval(substitute(subset), envir = dd, enclos = parent.frame())
	gmod <- attr(dd, "global")
	calls <- attr(dd, "calls")[subset]

	arg <- list(substitute(gmod), NA, ...)
	env <- attr(tryCatch(terms(gmod), error=function(...) terms(formula(gmod))),".Environment")

	models <- lapply(calls, eval, envir=env)

	attr(models, "rank.call") <- attr(dd, "rank.call")
	attr(models, "rank") <- attr(dd, "rank")

	return(models)
}
