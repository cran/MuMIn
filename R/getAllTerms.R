`getAllTerms.default` <-
function(x, ...) getAllTerms(as.formula(formula(x)), ...)

`getAllTerms.terms` <-
function(x, offset = TRUE, ...) {
	if (!is.null(attr(x, "offset"))){
		offs <- sapply((attr(x, "variables")[-1])[attr(x, "offset")], deparse)
	} else {
		offs <- NULL
	}
	ret <- attr(x, "term.labels")

	# Get term names, with higher order term components arranged alphabetically
	if (length(ret) > 0) {
		factors <- attr(x, "factors")
		factors1 <- rownames(factors)
		ret <- apply(factors > 0, 2, function(i) paste(sort(factors1[i]), collapse=":"))
	}

	# Leave out random terms (lmer type)
	#ran <- attr(x, "variables")[-1][-c(attr(x, "offset"), attr(x, "response"))]
	ran <- attr(x, "variables")[-1]
	ran <- as.character(ran[sapply(ran,
		function(x) length(x) == 3 && x[[1]] == as.name("|"))])
	ifx <- !(ret %in% ran)

	ret <- ret[ifx]

	# finally, sort by order and then alphabetically
	ret <- unname(ret[order(attr(x, "order")[ifx], ret)])


	if (!is.null(offs[1])) {
		if (offset)
			ret <- c(ret, offs)
		attr(ret, "offset") <- offs
	}
	attr(ret, "intercept") <- attr(x, "intercept")
	if (length(ran) > 0) {
		attr(ret, "random.terms") <- ran
		attr(ret, "random") <- reformulate(c(".", paste("(", ran, ")",
			sep = "")), response = ".")
	}

	return(ret)
}

`getAllTerms.formula` <-
function(x, ...) getAllTerms.terms(terms(x), ...)

`getAllTerms.lme` <-
function(x, ...) {
	ret <- getAllTerms(terms(x))
	attr(ret, "random") <- . ~ .

	# Code from nlme:::print.reStruct, modified slightly
	reStruct <- x$modelStruct$reStruct
	nobj <- length(reStruct)
	if (is.null(namx <- names(reStruct)))
		names(reStruct) <- nobj:1
	aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
	aux[lower.tri(aux)] <- ""
	reStruct[] <- rev(reStruct)
	aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
	aux[lower.tri(aux)] <- ""
	attr(ret, "random.terms") <- paste(lapply(lapply(reStruct, attr, "formula"),
		"[[", 2), "|",
		rev(apply(aux, 1, function(z) paste(z[z != ""], collapse = " %in% "))))

	return(ret)
}

`getAllTerms.glmer` <- # For backwards compatibility
`getAllTerms.lmer` <-  # with older versions of lme4
`getAllTerms.mer` <-
function(x, ...) getAllTerms(formula(x), ...)

# Apparently there is no (explicit) intercept in coxph, but 'terms' gives
# attr(,"intercept") == 1.
`getAllTerms.coxph` <- function (x, ...) {
	ret <- getAllTerms.default(x)
	attr(ret,"intercept") <- 0
	return(ret)
}


`getAllTerms` <-
function(x, ...) UseMethod("getAllTerms")
