`get.models` <-
function(dd, subset = delta <= 4, ...) {

	subset <- eval(substitute(subset), envir = dd, enclos = parent.frame())
	gmod <- attr(dd, "global")
	frm <- attr(dd, "formulas")[subset]

	sgmod <- substitute(gmod)
	models <- lapply(frm, function(.x) eval(call("update", sgmod, .x), sys.parent(3)))

	if (!is.null(attr(dd, "rank.call"))) {
  		attr(models, "rank.call") <- attr(dd, "rank.call")
	}

	return(models)
}

