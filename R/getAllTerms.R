`getAllTerms.default` <-
#function(x, ...) getAllTerms.formula(as.formula(formula(x)), ...)
function(x, ...) getAllTerms.terms(terms(as.formula(formula(x))), ...)


`getAllTerms.lm` <-
function(x, intercept = FALSE, ...)
	if(inherits(x, "gam"))
		getAllTerms.terms(terms(formula(x), ...), intercept = intercept) else
		getAllTerms.terms(terms(x, ...), intercept = intercept)

`getAllTerms.terms` <-
function(x, offset = TRUE, intercept = FALSE, ...) {

	interceptLabel <- "(Intercept)"
	variables <- attr(x, "variables")[-1L]

	if (!is.null(attr(x, "offset"))){
		offs <- sapply(variables[attr(x, "offset")], deparse)
	} else offs <- NULL

	ret <- attr(x, "term.labels")

	# Get term names, with higher order term components arranged alphabetically
	if (length(ret) > 0) {
		factors <- attr(x, "factors")
		factors <- factors[order(rownames(factors)), , drop = FALSE]
		v <- rownames(factors)
		ret <- apply(factors != 0L, 2L, function(x) paste(v[x], collapse = ":"))
	}

	# Leave out random terms (lmer type)
	#ran <- attr(x, "variables")[-1][-c(attr(x, "offset"), attr(x, "response"))]
	ran <- variables
	ran <- as.character(ran[sapply(ran,
		function(x) length(x) == 3L && x[[1L]] == "|")])
	ifx <- !(ret %in% ran)

	ret <- ret[ifx] # ifx - indexes of fixed terms
	#retUnsorted <- ret

	# finally, sort by order and then alphabetically
	#ret <- unname(ret[order(attr(x, "order")[ifx], ret)])
	ord <- order(attr(x, "order")[ifx], gsub("I\\((.*)\\)", "\\1", ret))
	ret <- unname(ret[ord])

	if(intercept && attr(x, "intercept")) {
		ret <- c(interceptLabel, ret)
		ord <- c(1, ord + 1L)
	}

	if (!is.null(offs[1L])) {
		if (offset) {
			ret <- c(ret, offs)
			ord <- c(ord, length(ord) + 1L)
		}
		attr(ret, "offset") <- offs
	}
	attr(ret, "intercept") <- attr(x, "intercept")
	attr(ret, "interceptLabel") <- interceptLabel

	if (length(ran) > 0L) {
		attr(ret, "random.terms") <- ran
		attr(ret, "random") <- reformulate(c(".", paste("(", ran, ")",
			sep = "")), response = ".")
	}

	response <- attr(x, "response")
	response <- if(response == 0L) NULL else variables[[response]]
	attr(ret, "response") <- response
	attr(ret, "order") <- order(ord)

	return(ret)
}

`getAllTerms.formula` <-
function(x, ...) getAllTerms.terms(terms.formula(x), ...)

`getAllTerms.lme` <-
function(x, ...) {
	ret <- getAllTerms.terms(terms(x), ...)
	attr(ret, "random") <- . ~ .

	# Code from nlme:::print.reStruct, modified slightly
	reStruct <- x$modelStruct$reStruct
	nobj <- length(reStruct)
	if (is.null(namx <- names(reStruct)))
		names(reStruct) <- nobj:1L
	aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
	aux[lower.tri(aux)] <- ""
	reStruct[] <- rev(reStruct)
	aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
	aux[lower.tri(aux)] <- ""
	attr(ret, "random.terms") <- paste(lapply(lapply(reStruct, attr, "formula"),
		"[[", 2L), "|",
		rev(apply(aux, 1L, function(z) paste(z[z != ""], collapse = " %in% "))))

	return(ret)
}

# `getAllTerms.glmer` <- # For backwards compatibility
# `getAllTerms.lmer` <-  # with older versions of lme4
`getAllTerms.mer` <-
function(x, ...) getAllTerms(lme4::formula(x), ...)

# Apparently there is no (explicit) intercept in coxph, but 'terms' gives
# attr(,"intercept") == 1.
`getAllTerms.coxph` <- function (x, ...) {
	ret <- getAllTerms.default(x, ...)
	attr(ret, "intercept") <- 0L
	attr(ret, "interceptLabel") <- NULL
	return(ret)
}

`getAllTerms.glmmML` <- function (x, ...) {
	ret <- getAllTerms.terms(terms(x), ...)
	#attr(ret, "random.terms") <-  deparse(call("|", 1,  x$call$cluster))
	attr(ret, "random.terms") <-  paste("1 |",  x$call$cluster)
	ret
}

`getAllTerms` <-
function(x, ...) UseMethod("getAllTerms")
