`getAllTerms.default` <-
function(x, ...) {
	return(getAllTerms(as.formula(formula(x))))
}

`getAllTerms.formula` <-
function(x, ...) {
	mTerms <- terms(x)
	ret <- attr(terms(x),"term.labels")
	if (length(ret) > 0) {
		ret <- ret[order(ret)]
		i <- grep(" ", ret)
		ret[i] <- paste("(", ret[i] , ")")

		mTerms <- terms(as.formula(paste(". ~", paste(ret, sep=" ", collapse=" + "))))
		ret <- attr(mTerms, "term.labels")
	}

	attr(ret, "intercept") <- attr(mTerms, "intercept")
	ret
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
function (x, ...) UseMethod("getAllTerms")