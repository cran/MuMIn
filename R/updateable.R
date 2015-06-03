`updateable` <-
function (FUN, eval.args = NULL, Class) {
	Fname <- substitute(FUN)
	FUN <- match.fun(FUN)
	FUNV <- function() {
		ocl <- cl <- match.call()
		parentframe <- parent.frame()
		for(i in eval.args) ocl[[i]] <- eval(ocl[[i]], parentframe) # [[4]]
		cl[[1L]] <- Fname
		res <- eval(cl, parentframe)
		if(!isS4(res) && is.list(res))
			res$call <- ocl else
			attr(res, "call") <- ocl
		class(res) <- Class  # [[8]]
		res
	}
	
	body(FUNV)[c(if(missing(eval.args)) 4L, if(missing(Class)) 8L)] <- NULL
    formals(FUNV) <- formals(FUN)
	rm(FUN, inherits = FALSE)
    FUNV
}

`updateable2` <-
function (FUN, Class) .Defunct("updateable")


`get_call` <- function(x) {
	rval <-
	if(isS4(x)) {
		if(any(i <- (sln <- c("call", "CALL", "Call")) %in% slotNames(x)))
			slot(x, sln[i][1L]) else
			if(!is.null(attr(x, "call")))
				attr(x, "call") else NULL
	} else {
		if(!is.atomic(x) && (i <- match("call", names(x), nomatch = 0L)) != 0L) {
			x[[i]]
		} else if(!is.null(attr(x, "call"))) {
			attr(x, "call")
		} else
			NULL
	}
	if(is.null(rval)) stats::getCall(x) else rval
}

##==============================================================================

`uGamm` <-
 function(formula, random = NULL, ..., lme4 = inherits(random, "formula")) {
	pkg <- if(lme4) "gamm4" else "mgcv"
	if (!require(pkg, character.only = TRUE)) stop("cannot load package ", sQuote(pkg))
	funcl <- call("get", if(lme4) "gamm4" else "gamm", ns <- asNamespace(pkg))
	clx <- cl <- match.call()
	clx$lme4 <- NULL	
	clx <- match.call(clx, definition = eval(funcl, envir = ns))
	clx[[1L]] <- funcl
	res <- eval.parent(clx)
	res$call <- cl
	class(res) <- c(if(lme4) "gamm4", "gamm")
	res
}


`gamm` <- 
function(...) .Deprecated("uGamm", old = "MuMIn::gamm")
