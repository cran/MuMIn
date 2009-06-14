`print.averaging` <-
function(x, ...) {
	cat("\nModel summary:\n")
	print(signif(x$summary,3))

	cat("\nVariables:\n")
	print(x$variable.codes, quote= F)

	cat("\nAveraged model parameters:\n")
	print(signif(x$avg.model, 3))

	cat("\nRelative variable importance:\n")
	print(round(x$relative.importance, 2))
}

