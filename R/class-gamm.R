# gamm/gamm4 support


`update.gamm` <-
function(object, ...) {
# or, if call is as attribute: object$call <- attr(object, "call")
	update.default(object, ...)
}

`print.gamm` <-
function(x, ...) {
	cat("\nCall:\n", paste(asChar(x$call, nlines = -1L), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
	cat("--- \n")
	print(x[[if(inherits(x, "gamm4")) "mer" else "lme"]])
	cat("--- \n")
	print(x$gam)
	invisible(x)
}

`formula.gamm` <-
function (x, ...) formula(x$gam, ...)