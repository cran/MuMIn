.onLoad <- function(libname, pkgName) {
	# Ugly tricks to put own replacement methods on top (don't try this at home):
	asNeeded <- 
    function(pkgName, fun) {
		if(paste("package", pkgName, sep = ":") %in% search()) fun() else 
			setHook(packageEvent(pkgName, "attach"), fun)
	}
	
	regmethod <- 
    function(funname, classname, s4 = FALSE,
		fun = get(paste0(funname, ".", classname), getNamespace(.packageName)),
		envir = .GlobalEnv)
		do.call(if(s4) methods::setMethod else "registerS3method", 
            list(funname, classname, fun),
			envir = envir)
          
    regMethodsOnPkgAttach <-   
    function(pkgName, funname, classname, s4 = FALSE)
        asNeeded(pkgName, function(...) for(a in classname) regmethod(funname, a, s4))
	    
    regMethodsOnPkgAttach("unmarked", "logLik", "unmarkedFit", TRUE)
    regMethodsOnPkgAttach("nlme", "predict", c("lme", "gls"))
    regMethodsOnPkgAttach("xtable", "xtable", 
        c("summary.averaging", "averaging", "model.selection"))
 
	
}
