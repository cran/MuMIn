`null.fit` <- function(x, evaluate = FALSE, envir = environment(as.formula(formula(x)))) {
	cl <- .getCall(x)
	mClasses <- c("glmmML", "lm", "lme", "gls", "mer", "unmarkedFit")
	mClass <- mClasses[inherits(x, mClasses, which = TRUE) != 0L][1]
	if(is.na(mClass)) mClass <- "default"
	formulaArgName <- "formula"
	Fun <- "glm"
	switch(mClass,
		glmmML = {
			REML <- FALSE
			if(is.null(cl$family)) cl$family <- as.name("binomial")
		}, gls = {
			arg <- formals(match.fun(cl[[1]]))
			REML <- match.arg(cl$method, eval(arg$method)) == "REML"
			formulaArgName <- "model"
			cl$weights <- NULL
		}, lme = {
			arg <- formals(match.fun(cl[[1L]]))
			REML <- match.arg(cl$method, eval(arg$method)) == "REML"
			formulaArgName <- "fixed"
			cl$weights <- NULL
		}, mer = {
			arg <- formals(match.fun(cl[[1L]]))
			REML <- if(is.null(cl$REML)) arg$REML else cl$REML
		}, unmarkedFit = {
			nm <- names(cl)[-1]
			if("formula" %in% nm) {
				cl$formula <- ~1~1
			} else {
				formula.arg <- nm[grep(".+formula$", nm[1:7])]
				for (i in formula.arg) cl[[i]] <- ~1
			}
			cl$starts <- NULL
			Fun <- NA
		}, {
			REML <- FALSE
		}
	)

	if(!is.na(Fun)) cl[[1L]] <- as.name(Fun)
	if(identical(Fun, "glm")) {
		if(formulaArgName != "formula")
			names(cl)[names(cl) == formulaArgName] <- "formula"
		cl$formula <- update(as.formula(cl$formula), . ~ 1)
		cl$method <- cl$start <- cl$offset <- contrasts <- NULL
		#cl <- cl[c(TRUE, names(cl)[-1] %in% names(formals(match.fun(cl[[1L]]))))]
		cl <- cl[c(TRUE, names(cl)[-1] %in% names(formals("glm")))]
	}
	if(evaluate) eval(cl, envir = envir) else cl
}

#
##lm1 <- lm(weight ~ group)
##plot(weight ~ group)
#glm1 <- glm(exp(weight) ~ group, Gamma("log"))
#
#glm2 <- local ({
#glm(exp(weight) ~ group, Gamma("log"))
#})
#environment(formula(glm1))
#environment(formula(glm2))
#r.squaredLR(glm1)
#r.squaredLR(lm1)
#nullCall(lm1)
#summary(lm1)
# Multiple R-squared: 0.07308,	Adjusted R-squared: 0.02158
#str(summary(lm(weight ~ group)))
# r.squaredLR(lm1)
#r.squaredLR(x, null = nullCall(x, TRUE, parent.frame()))
#nullCall(x, evaluate = FALSE, envir = parent.frame())

#args(null.fit)


`r.squaredLR` <- function(x, null = null.fit(x, TRUE)) {
	L0 <- c(logLik(null))
	L1 <- c(logLik(x))
	n <- nobs(x)
	ret <- 1 - exp(-2 / n * (L1 - L0))
	max.r2 <- 1 - exp(2 / n * L0)
	attr(ret, "adj.r.squared") <- ret / max.r2
	ret
}
