`null.fit` <- function(x, evaluate = FALSE, envir = environment(as.formula(formula(x)))) {
	cl <- .getCall(x)
	mClasses <- c("glmmML", "lm", "lme", "gls", "mer", "unmarkedFit", "coxph",
		"coxme")
	mClass <- mClasses[inherits(x, mClasses, which = TRUE) != 0L][1]
	if(is.na(mClass)) mClass <- "default"
	formulaArgName <- "formula"
	Fun <- "glm"
	call2arg <- function(x) formals(match.fun(x[[1L]]))
	switch(mClass,
		glmmML = {
			REML <- FALSE
			if(is.null(cl$family)) cl$family <- as.name("binomial")
		}, gls = {
			REML <- match.arg(cl$method, eval(call2arg(cl)$method)) == "REML"
			formulaArgName <- "model"
			cl$weights <- NULL
		}, lme = {
			REML <- match.arg(cl$method, eval(call2arg(cl)$method)) == "REML"
			formulaArgName <- "fixed"
			cl$weights <- NULL
		}, mer = {
			arg <- formals(match.fun(cl[[1L]]))
			REML <- if(is.null(cl$REML)) call2arg(cl)$REML else cl$REML
		}, unmarkedFit = {
			nm <- names(cl)[-1L]
			if("formula" %in% nm) {
				cl$formula <- ~1~1
			} else {
				formula.arg <- nm[grep(".+formula$", nm[1:7])]
				for (i in formula.arg) cl[[i]] <- ~1
			}
			cl$starts <- NULL
			Fun <- NA
		}, coxph = , coxme = {
			Fun <- "coxph"
			cl$formula <- update.formula(eval(cl$formula), . ~ 1)
		}, lm = {
			REML <- FALSE
		}, {
			stop("do not know (yet) how to construct a null model for class ",
				paste(sQuote(class(x))))
		}
	)

	if(!is.na(Fun)) cl[[1L]] <- as.name(Fun)
	if(identical(Fun, "glm")) {
		if(formulaArgName != "formula")
			names(cl)[names(cl) == formulaArgName] <- "formula"
		cl$formula <- update(as.formula(cl$formula), . ~ 1)
		cl$method <- cl$start <- cl$offset <- contrasts <- NULL
	}
	cl <- cl[c(TRUE, names(cl)[-1L] %in% names(call2arg(cl)))]
	if(evaluate) eval(cl, envir = envir) else cl
}



`r.squaredLR` <- function(x, null = null.fit(x, TRUE, parent.frame())) {
	L0 <- c(logLik(null))
	L1 <- c(logLik(x))
	n <- nobs(x)
	ret <- 1 - exp(-2 / n * (L1 - L0))
	max.r2 <- 1 - exp(2 / n * L0)
	attr(ret, "adj.r.squared") <- ret / max.r2
	ret
}
