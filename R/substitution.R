# evaluate 'expr' in 'env' after adding variables passed as '...'
evalExprInEnv <- function(expr, env, enclos, ...) {
	list2env(list(...), env)
	eval(expr, envir = env, enclos = enclos)
}

# change `names[]` for varName[1], varName[2], ... in expression
# Not using `substitute` anymore to omit function calls.
# Ignore also expressions within I(), elements extracted with $ or @
`.subst.names.for.items` <-
function(expr, names, varName, n = length(names), fun = "[") {
	exprApply(expr, names, symbols = TRUE,
		function(x, v, fun, varName, parent) {
			if(is.call(parent) && any(parent[[1L]] == c("I", "$", "@")))
				return(x)
			if(length(x) == 1L)
				return(call(fun, varName, match(asChar(x), v)))
			x
		}, v = names, fun = fun, varName = as.name(varName))
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
.subst.term <- function(x) {
	if(length(x) < 2L) cry(x, "'Term' needs one argument")
	as.name(asChar(x[[2L]]))
}

.subst.with <- 
function (x, fac, allTerms, vName, envir = parent.frame()) {
    if (length(x) > 4L) cry(x, "too many arguments [%d]", length(x) - 1L)
    if (length(x[[2L]]) == 2L && x[[2L]][[1L]] == "+") {
        fun <- "all"
        sx <- asChar(x[[2L]][[2L]], backtick = FALSE)
    } else {
        fun <- "any"
        sx <- asChar(x[[2L]], backtick = FALSE)
    }
    dn <- dimnames(fac)
    if (!(sx %in% dn[[2L]])) cry(x, "unknown variable name '%s'", sx)
    xorder <- if(length(x) >= 3L) as.integer(eval(x[[3L]], envir))
		else unique(rowSums(fac))
    i <- which(fac[, sx])
    j <- which(is.element(rowSums(fac[i, , drop = FALSE]), xorder))
    if(length(j) == 0L) cry(x, "no terms match the criteria")    
    as.call(c(as.name(fun), call("[", vName, as.call(c(as.name("c"), match(dn[[1L]][i[j]], allTerms))))))
}

.subst.vars.for.args <- function(e) {
	for(i in 2L:length(e))
		if(!is.name(e[[i]]))
			e[[i]] <- as.name(asChar(e[[i]]))
	e
}

.subst.has <- 
function(e) {
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

.subst.has.dc <- 
function(e) {
		for(i in 2L:length(e)) e[[i]] <- call("has", e[[i]])
		e
}

.subst.v <- 
function(x, cVar, fn) {
	if(length(x) > 2L) cry(x, "discarding extra arguments", warn = TRUE)
	i <- which(fn == x[[2L]])[1L]
	if(is.na(i)) cry(x, "'%s' is not a valid name of 'varying' element",
					  as.character(x[[2L]]), warn = TRUE)
	call("[[", cVar, i)
}

# substitute function calls in 'e'. 'func' must take care of the substitution job.
`exprapply0` <- 
function(e, name, func, ...)
exprApply(e, name, func, ..., symbols = FALSE)

`exprApply` <-
function (expr, what, FUN, ..., symbols = FALSE) {
    FUN <- match.fun(FUN)
	if(all(names(formals(FUN)) != "parent"))
		formals(FUN)[["parent"]] <- NA
	.exprapply(expr, what, FUN, ..., symbols = symbols)
}

`.exprapply` <-
function (expr, what, FUN, ..., symbols = FALSE, parent = NULL) {
	self <- sys.function()
	if((ispairlist <- is.pairlist(expr)) || is.expression(expr)) {
		for (i in seq_along(expr))	expr[i] <-
			list(self(expr[[i]], what, FUN, ..., symbols = symbols, parent = expr))
		return(if(ispairlist) as.pairlist(expr) else expr)
	}
    n <- length(expr)
    if (n == 0L)
		return(expr) else
	if (n == 1L) {
		if (!is.call(expr)) {
            if (symbols && (anyNA(what) || any(expr == what)))
                expr <- FUN(expr, ..., parent = parent)
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
			y <- self(expr[[i]], what, FUN, ..., symbols = symbols, parent = expr)
			if(!missing(y)) expr[i] <- list(y)
		}
    }
    if (anyNA(what) || (length(expr[[1L]]) == 1L && any(expr[[1L]] == what)))
		expr <- FUN(expr, ..., parent = parent)
    return(expr)
}

