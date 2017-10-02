`.DebugPrint` <-
function (x) {
    if (isTRUE(getOption("debug.MuMIn"))) {
		fun <- asChar(sys.call(sys.parent())[[1L]])
		name <- substitute(x)
		cat(sprintf("<%s> ~ ", fun))
		if(is.language(name)) cat(asChar(name), "= \n")
        print(x)
    }
}
	
`.Debug` <- 
function(expr) {
	if(isTRUE(getOption("debug.MuMIn"))) {
		eval.parent(substitute(expr))
	}
}

`srcc` <- function() {
	ret <- eval(expression(source("clipboard", local = TRUE)), .GlobalEnv)
	return(if(ret$visible) ret$value else invisible(ret$value))
}
