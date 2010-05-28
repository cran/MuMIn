`print.model.selection` <-
function(x, ...) {
	x$weight <- round(x$weight, 3)

	nn <- attr(x, "terms")

	names(x)[seq(along=nn)] <- sapply( strsplit(nn, ":"), function(xx) paste(sapply(xx, abbreviate, 6 )  , collapse=":") )

	cat ("Model selection table", "\n")
	print(signif(as.matrix(x), digits=4), na.print="")
	if (!is.null(attr(x, "random.terms"))) {
		cat("Random terms:", paste(attr(x, "random.terms"), collapse=", "), "\n")
 	}
}

