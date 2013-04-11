.onLoad <- function(libname, pkgName) {
	# cat(sprintf("onLoad(%s, %s) \n", libname, pkgName))

	asNeeded <- function(pkgName, Fun) {
		if(paste("package", pkgName, sep = ":") %in% search()) Fun() else 
			setHook(packageEvent(pkgName, "attach"), Fun)
	}
	
	asNeeded("unmarked", function(...)
		do.call("setMethod", list("logLik", "unmarkedFit",  logLik.unmarkedFit),
				envir = .GlobalEnv))

}
