`get.models` <-
function(object, subset = delta <= 4, ...) {

	#subset <- if (missing(subset)) quote() else substitute(subset)
	subset <- eval(substitute(subset), envir = object, enclos = parent.frame())
	gmod <- attr(object, "global")
	calls <- attr(object, "calls")[subset]

	arg <- list(substitute(gmod), NA, ...)
	env <- attr(tryCatch(terms(gmod), error=function(...) terms(formula(gmod))),
		".Environment")

	models <- lapply(calls, eval, envir=env)

	attr(models, "rank.call") <- attr(object, "rank.call")
	attr(models, "rank") <- attr(object, "rank")

	return(models)
}
