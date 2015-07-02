# evaluate 'expr' in 'env' after adding variables passed as '...'
evalExprInEnv <- function(expr, env, enclos, ...) {
	list2env(list(...), env)
	eval(expr, envir = env, enclos = enclos)
}

# substitute names for varName[1], varName[2], ... in expression
`.subst4Vec` <- function(expr, names, varName, n = length(names), fun = "[") {
	varName <- as.name(varName)
	eval(call("substitute", expr,
		env = structure(lapply(seq_len(n), function(i) call(fun, varName, i)),
						names = names)),
		envir = NULL)
}

# like substitute, but does evaluate 'expr'.
subst <-
function(expr, envir = NULL, ...) {
	eval.parent(call("substitute", expr, c(envir, list(...))))
}

asChar <- function(x, control = NULL, nlines = 1L, ...) 
	if(is.character(x)) x[1L:nlines] else
	deparse(x, control = control, nlines = nlines, ...)

## .sub_* functions used with '.exprapply' as 'func'
.sub_Term <- function(x) {
	if(length(x) < 2L) cry(x, "'Term' needs one argument")
	as.name(asChar(x[[2L]]))
}

.sub_dot <- function(x, fac, at, vName) {
	if(length(x) != 2L) cry(x, "exactly one argument needed, %d given.", length(x) - 1L)
	if(length(x[[2L]]) == 2L && x[[2L]][[1L]] == "+") {
		fun <- "all"
		sx <- as.character(x[[2L]][[2L]])
	} else {
		fun <- "any"
		sx <- as.character(x[[2L]])
	}
	dn <- dimnames(fac)
	if(!(sx %in% dn[[2L]])) cry(x, "unknown variable name '%s'", sx)
	as.call(c(as.name(fun), call("[", vName, as.call(c(as.name("c"), 
		match(dn[[1L]][fac[, sx]], at))))))
}

.sub_dc_has <- function(e, fname) {
		#e[[1L]] <- call(":::", as.name(.packageName), fname)
		e[[1L]] <- fname
		for(i in 2L:length(e)) e[[i]] <- call("has", e[[i]])
		e
}

.sub_args_as_vars <- function(e) {
	for(i in 2L:length(e))
		if(!is.name(e[[i]]))
			e[[i]] <- as.name(asChar(e[[i]]))
	e
}

.sub_has <- function(e) {
	n <- length(e)
	for(i in seq.int(2L, n)) {
		ex <- if(length(e[[i]]) == 2L && e[[i]][[1L]] == "!")
			call("is.na", e[[i]][[2L]]) else
			call("!", call("is.na", if(is.name(e[[i]])) e[[i]] else
						   as.name(asChar(e[[i]]))))
		res <- if(i == 2L) ex else call("&", res, ex)
	}
	call("(", res)
}

.sub_V <- function(x, cVar, fn) {
	if(length(x) > 2L) cry(x, "discarding extra arguments", warn = TRUE)
	i <- which(fn == x[[2L]])[1L]
	if(is.na(i)) cry(x, "'%s' is not a valid name of 'varying' element",
					  as.character(x[[2L]]), warn = TRUE)
	call("[[", cVar, i)
}

# substitute function calls in 'e'. 'func' must take care of the substitution job.
`.exprapply` <- function(e, name, func, ...) 
exprApply(e, name, func, ..., symbols = FALSE)


`exprApply` <-
function (expr, what, FUN = print, ..., symbols = FALSE) {
    FUN <- match.fun(FUN)
	self <- sys.function()
	if((ispairlist <- is.pairlist(expr)) || is.expression(expr)) {
		for (i in seq_along(expr))	expr[i] <-
			list(self(expr[[i]], what, FUN, ..., symbols = symbols))
		return(if(ispairlist) as.pairlist(expr) else expr)
	}
    n <- length(expr)
    if (n == 0L)
		return(expr) else
	if (n == 1L) {
		if (!is.call(expr)) {
            if (symbols && (is.na(what) || any(expr == what))) 
                expr <- FUN(expr, ...)
            return(expr)
        }
    } else {
		if(expr[[1L]] == "function") {
			if(n == 4L) {
				n <- 3L
				expr[[4L]] <- NULL ## remove srcref
			}
		} 
        for (i in seq.int(2L, n)) {
			y <- self(expr[[i]], what, FUN, symbols = symbols, ...)
			expr[i] <- list(y)
		}
    }
    if (is.na(what) || any(expr[[1L]] == what)) expr <- FUN(expr, ...)
    return(expr)
}
