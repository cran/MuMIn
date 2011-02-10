`get.models` <-
function(dd, subset = delta <= 4, ...) {

	#subset <- if (missing(subset)) quote() else substitute(subset)
	subset <- eval(substitute(subset), envir = dd, enclos = parent.frame())
	gmod <- attr(dd, "global")
	frm <- attr(dd, "formulas")[subset]

	arg <- list(substitute(gmod), NA, ...)
	env <- attr(tryCatch(terms(gmod), error=function(...) terms(formula(gmod))),".Environment")

	models <- lapply(frm, function(.x) {
		arg[[2]] <- .x
		do.call("update", arg, quote = TRUE, envir=env)
	})

	if (!is.null(attr(dd, "rank.call"))) {
  		attr(models, "rank.call") <- attr(dd, "rank.call")
	}

	return(models)
}
