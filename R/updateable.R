`updateable` <-
function (FUN) {
	Fname <- substitute(FUN)
	FUN <- match.fun(FUN)
	FUNV <- function() {
		ocl <- cl <- match.call()
		cl[[1L]] <- Fname
		res <- eval(cl, parent.frame())
		if(!isS4(res) && is.list(res))
			res$call <- ocl else
			attr(res, "call") <- ocl
		res
	}
    formals(FUNV) <- formals(FUN)
    FUNV
}

`updateable2` <-
function (FUN, Class) {
	Fname <- substitute(FUN)
	FUN <- match.fun(FUN)
	FUNV <- function() {
		ocl <- cl <- match.call()
		cl[[1L]] <- Fname
		res <- eval(cl, parent.frame())
		if(!isS4(res) && is.list(res))
			res$call <- ocl else
			attr(res, "call") <- ocl
		class(res) <- Class
		res
	}
	if(missing(Class)) body(FUNV)[[6L]] <- NULL
    formals(FUNV) <- formals(FUN)
	rm(FUN)
    FUNV
}


`.getCall` <- function(x) {
	if(isS4(x)) {
		if(any(i <- (sln <- c("call", "CALL", "Call")) %in% slotNames(x)))
			slot(x, sln[i][1L]) else
			if(!is.null(attr(x, "call")))
				attr(x, "call") else NULL
	} else {
		if(!is.atomic(x) && !is.null(x$call)) {
			x$call
		} else if(!is.null(attr(x, "call"))) {
			attr(x, "call")
		} else
			NULL
	}
}





`getCall.default` <-
function (x, ...)
.getCall(x)



##==============================================================================

`uGamm` <-
 function(formula, random = NULL, ..., lme4 = inherits(random, "formula")) {
	pkg <- if(lme4) "gamm4" else "mgcv"
	if (!require(pkg, character.only = TRUE)) stop("'gamm' requires package '",
												pkg, "' to be installed")
	funcl <- call("get", if(lme4) "gamm4" else "gamm", ns <- asNamespace(pkg))
	clx <- cl <- match.call()
	clx$lme4 <- NULL	
	clx <- match.call(clx, definition = eval(funcl, envir = ns))
	clx[[1L]] <- funcl
	#list(clx, cl)
	res <- eval(clx, parent.frame())
	res$call <- cl
	class(res) <- c(if(lme4) "gamm4", "gamm", "list")
	res
 }


`gamm` <- 
function(...) .Deprecated("uGamm", old = "MuMIn::gamm")


