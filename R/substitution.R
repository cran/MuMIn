# substitute has(a, !b, ...) for !is.na(a) & is.na(b) ..., in expression
`.substHas` <- function(e) {
	if(is.expression(e)) e <- e[[1L]]
	n <- length(e)
	if(n == 1L) return(e)
	if(e[[1L]] != "has") {
		for(i in 1L:n) e[[i]] <- .substHas(e[[i]])
		return(e)
	}
	res <- NULL
	for(i in seq.int(2L, n)) {
		ex <- if(length(e[[i]]) == 2L && e[[i]][[1L]] == "!")
			call("is.na", e[[i]][[2L]]) else
			call("!", call("is.na", if(is.name(e[[i]])) e[[i]] else
						   as.name(deparse(e[[i]], control = NULL))))

		res <- if(i == 2L) ex else call("&", res, ex)
	}
	res <- call("(", res)
	return(res)
}

# substitute function calls in 'e'. 'func' must take care of the substitution job.
`.substFunc` <- function(e, name, func = identity, ...) {
	if(is.expression(e)) e <- e[[1L]]
	n <- length(e)
	if(n == 0L) return(e) else if (n == 1L) {
		if (!is.call(e)) return(e)
	} else for(i in 2L:n) e[i] <- list(.substFunc(e[[i]], name, func, ...))
	if(any(e[[1L]] == name)) e <- func(e, ...)
	return(e)
}

# evaluate 'expr' in 'env' after adding variables passed as '...'
.evalExprIn <- function(expr, env, enclos, ...) {
	list2env(list(...), env)
	eval(expr, envir = env, enclos = enclos)
}

# substitute names for varName[1], varName[2], ... in expression
`.subst4Vec` <- function(expr, names, varName, n = length(names), fun = "[") {
	eval(call("substitute", expr,
		env = structure(lapply(seq_len(n), function(i) call(fun, varName, i)),
						names = names)),
		envir = NULL)
}


## .sub_* functions used with '.substFunc' as 'func'
.sub_Term <- function(x) {
	if(length(x) < 2L) .cry(x, "'Term' needs one argument")
	as.name(deparse(x[[2L]], control = NULL))
}


.sub_dot <- function(x, fac, at, vName) {
	if(length(x) != 2L) .cry(x, "exactly one argument needed, %d given.", length(x) - 1L)
	if(length(x[[2L]]) == 2L && x[[2L]][[1L]] == "+") {
		fun <- "all"
		sx <- as.character(x[[2L]][[2L]])
	} else {
		fun <- "any"
		sx <- as.character(x[[2L]])
	}
	dn <- dimnames(fac)
	if(!(sx %in% dn[[2L]])) .cry(x, "unknown variable name '%s'", sx)
	as.call(c(as.name(fun), call("[", vName, as.call(c(as.name("c"), 
		match(dn[[1L]][fac[, sx]], at))))))
}

.sub_dc_has <- function(e, fname) {
		e[[1L]] <- call(":::", as.name(.packageName), fname)
		for(i in 2L:length(e)) e[[i]] <- call("has", e[[i]])
		e
}

.sub_has <- function(e) {
	n <- length(e)
	for(i in seq.int(2L, n)) {
		ex <- if(length(e[[i]]) == 2L && e[[i]][[1L]] == "!")
			call("is.na", e[[i]][[2L]]) else
			call("!", call("is.na", if(is.name(e[[i]])) e[[i]] else
						   as.name(deparse(e[[i]], control = NULL))))
		res <- if(i == 2L) ex else call("&", res, ex)
	}
	call("(", res)
}

.sub_V <- function(x, cVar, fn) {
	if(length(x) > 2L) .cry(x, "discarding extra arguments", warn = TRUE)
	i <- which(fn == x[[2L]])[1L]
	if(is.na(i)) .cry(x, "'%s' is not a valid name of 'varying' element",
					  as.character(x[[2L]]), warn = TRUE)
	call("[[", cVar, i)
}

`exprApply` <-
function(expr, what, FUN = identity, ..., symbols = FALSE) {
	if(asExpr <- is.expression(expr)) expr <- expr[[1L]]
	n <- length(expr)
	if(n == 0L)
		return(expr)
	else if (n == 1L) {
		if(!is.call(expr)) {
			if (symbols && any(expr == what)) expr <- FUN(expr, ...)
			return(expr)
		}
	} else {
		for(i in seq.int(2L, n)) expr[i] <- 
			list(exprApply(expr[[i]], what, FUN, symbols =  symbols, ...))
	}
	if(any(expr[[1L]] == what)) expr <- FUN(expr, ...)
	if(asExpr) return(as.expression(expr)) 
	return(expr)
}
