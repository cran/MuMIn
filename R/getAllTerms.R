`getAllTerms.default` <-
function(x, ...) {
	return(getAllTerms(as.formula(formula(x))))
}

`getAllTerms.formula` <-
function(x, ...) {
	mTerms <- terms(x)
	intercept <- attr(mTerms, "intercept")
	if (!is.null(attr(mTerms, "offset"))){
		offs <-
		sapply(as.list(attr(mTerms,"variables")[-1])[attr(mTerms,"offset")],
		deparse)
	} else {
		offs <- NULL
	}
	ret <- attr(mTerms, "term.labels")
	if (length(ret) > 0) {
		ret <- ret[order(ret)]
		i <- grep(" ", ret)
		ret[i] <- paste("(", ret[i] , ")")

		mTerms <- terms(reformulate(ret))
		ret <- attr(mTerms, "term.labels")
	}

	attr(ret, "intercept") <- intercept

	if (!is.null(offs[1]))
		attr(ret, "offset") <- offs
	return(ret)
}

`getAllTerms.glmer` <-
function(x, ...) {
     ret <- getAllTerms(as.formula(formula(x)))
     i <- grep(" \\| ", ret)

     intercept <- attr(ret, "intercept")

	rnd <- ret[i]
     ret <- ret[-i]
     attr(ret, "random.terms") <- rnd
     rnd.formula <- paste("(", rnd, ")", sep="", collapse=" + ")
     rnd.formula <- as.formula(paste(". ~ .", rnd.formula, sep="+"))

	attr(ret, "random") <- rnd.formula
	attr(ret, "intercept") <- intercept

     ret
}

`getAllTerms.lme` <-
function(x, ...) {
	getAllTerms(as.formula(formula(x)))
}

`getAllTerms.lmer` <-
function(x, ...) {
     ret <- getAllTerms(as.formula(formula(x)))
     i <- grep(" \\| ", ret)

     intercept <- attr(ret, "intercept")

	rnd <- ret[i]
     ret <- ret[-i]
     attr(ret, "random.terms") <- rnd

     rnd.formula <- paste("(", rnd, ")", sep="", collapse=" + ")
     rnd.formula <- as.formula(paste(". ~ .", rnd.formula, sep="+"))

	attr(ret, "random") <- rnd.formula
	attr(ret, "intercept") <- intercept

     ret
}

`getAllTerms.mer` <-
function(x, ...) {
     ret <- getAllTerms(as.formula(formula(x)))
     i <- grep(" \\| ", ret)

     intercept <- attr(ret, "intercept")

	rnd <- ret[i]
     ret <- ret[-i]
     attr(ret, "random.terms") <- rnd

     rnd.formula <- paste("(", rnd, ")", sep="", collapse=" + ")
     rnd.formula <- as.formula(paste(". ~ .", rnd.formula, sep="+"))

	attr(ret, "random") <- rnd.formula
	attr(ret, "intercept") <- intercept

     ret
}

`getAllTerms` <-
function(x, ...) UseMethod("getAllTerms")
