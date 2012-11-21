.onLoad <- function(libname, pkgName) {
	# cat(sprintf("onLoad(%s, %s) \n", libname, pkgName))

# workaround for original logLik from unmarked with no nobs/df attributes
	if("package:unmarked" %in% search()) {
		do.call("setMethod", list("logLik", "unmarkedFit", logLik.unmarkedFit),
			envir=.GlobalEnv)
	} else {
		setHook(packageEvent("unmarked", "attach"), function(...) {
			do.call("setMethod", list("logLik", "unmarkedFit",
				logLik.unmarkedFit), envir=.GlobalEnv)
		})
	}
}
