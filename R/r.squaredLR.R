`null.fit` <-
function(object, evaluate = FALSE, RE.keep = FALSE, envir = NULL, ...) {
	# backward compatibility:
	if("x" %in% ...names()) {
		object <- ...elt(match("x", ...names()))
		warning("the argument ", sQuote("x"), " has been removed. Use ",
            sQuote("object"), " instead")
	}
    
	# TODO: detect if RE.keep is TRUE and object is not a mixed model

	cl <- get_call(object)
	if(!is.environment(envir)) envir <- environment(as.formula(formula(object)))
	
	if(RE.keep) {
		if(inherits(object, c("mer", "merMod", "coxme", "lmekin"))) {
			cl$formula <- .nullREForm(as.formula(cl$formula))
			environment(cl$formula) <- envir
		} else if(inherits(object, "gamm")) {
			mefm <- object[[if("lme" %in% names(object)) "lme" else "mer"]]
			
			if(inherits(mefm, "merMod")) {
				Fun <- if(inherits(mefm, "glmerMod"))
					"glmer" else if(inherits(mefm, "lmerMod")) {
					cl$family <- NULL
					"lmer"
				}
				cl$REML <- as.logical(object$mer@devcomp$dims[['REML']])
				frm <- cl$formula
				frm[[3L]] <- call("+", 1, as.formula(cl$random)[[2L]])
				cl$random <- NULL
				environment(cl$formula) <- envir
			} else if (inherits(mefm, "lme")) {
				Fun <- "lme"
				cl$fixed <- update.formula(as.formula(cl$formula), . ~ 1)
				cl$formula <- cl$family <- NULL
				cl$method <- object$lme$method
				environment(cl$fixed) <- envir
			}
			cl[[1L]] <- as.symbol(Fun)
		} else if(inherits(object, c("glmmML", "glimML"))) {
			cl$formula <- update.formula(as.formula(cl$formula), . ~ 1)
			environment(cl$formula) <- envir
		} else if(inherits(object, "lme")) {
			cl$fixed <- update.formula(as.formula(cl$fixed), . ~ 1)
			environment(cl$fixed) <- envir
		} else {
			stop("do not know (yet) how to construct a null model with RE for class ",
				 prettyEnumStr(class(object), sep.last = ", "))					
		}
		return(if(evaluate) eval(cl, envir = envir) else cl)
	}	
	
	mClasses <- c("glmmML", "lm", "lme", "gls", "mer", "merMod", "lmekin",
				  "unmarkedFit", "coxph", "coxme", "zeroinfl", "gamm",
                  "survreg")
	mClass <- mClasses[inherits(object, mClasses, which = TRUE) != 0L][1L]

	if(is.na(mClass)) mClass <- "default"
	formulaArgName <- "formula"
	Fun <- "glm"
	call2arg <- function(x) formals(match.fun(x[[1L]]))
	switch(mClass,
		glmmML = {
			if(is.null(cl$family)) cl$family <- as.name("binomial")
		}, gls = {
			formulaArgName <- "model"
			cl$weights <- NULL
		}, lme = {
			formulaArgName <- "fixed"
			cl$weights <- NULL
		}, lmekin =, merMod =, mer = {
			arg <- formals(match.fun(cl[[1L]]))
		}, unmarkedFit = {
			nm <- names(cl)[-1L]
			if("formula" %in% nm) {
				cl$formula <- ~1~1
			} else {
				formula.arg <- nm[grep(".+formula$", nm[1L:7L])]
				for (i in formula.arg) cl[[i]] <- ~1
			}
			cl$starts <- NULL
			Fun <- NA
		}, coxph =, coxme = {
			Fun <- "coxph"
			cl$formula <- update.formula(eval(cl$formula), . ~ 1)
		}, survreg = , zeroinfl =, lm = {
			Fun <- NA
			cl$formula <- update.formula(as.formula(cl$formula), . ~ 1)
		}, gamm = {
			Fun <- "gam"
			cl$formula <- update.formula(as.formula(cl$formula), . ~ 1)
			cl$random <- NULL
		}, {
			stop("do not know (yet) how to construct a null model for class ",
				sQuote(class(object)))
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
		

# from lme4:::findbars:
.findbars <- function (term) {
    if (is.name(term) || !is.language(term))
        return(NULL)
    if (term[[1L]] == as.name("("))
        return(.findbars(term[[2L]]))
    if (!is.call(term))
        stop("term must be of class call")
    if (term[[1L]] == as.name("|"))
        return(term)
    if (length(term) == 2L)
        return(.findbars(term[[2L]]))
    c(.findbars(term[[2L]]), .findbars(term[[3L]]))
}

`.nullREForm` <-
function(formula) {
	re <- lapply(.findbars(formula), function(x) call("(", x))
	f <- 1
	for(i in seq_along(re)) f <- call("+", f, re[[i]])
	formula[[length(formula)]] <- f
	formula
}

.getLLML <- 
function(x) {
    cls <- class(x)
    llfun <- if(isS4(logLik)) selectMethod("logLik", cls) else logLik
    if(isS3stdGeneric(llfun))
        for(cl in cls)
            if(is.function(llfun <- getS3method("logLik", cl, optional = TRUE)))
                break
    if(is.null(llfun))
        stop("no 'logLik' method found for object of class ",
            prettyEnumStr(cls, sep.last = ", "))
    arg <- list(object = x, REML = FALSE)
    do.call(llfun, arg[names(arg) %in%  names(formals(llfun))])
}


`r.squaredLR` <-
function(object, null = NULL, null.RE = FALSE, ...) {
	if("x" %in% ...names()) {
		object <- ...elt(match("x", ...names()))
		warning("the argument ", sQuote("x"), " has been removed. Use ",
            sQuote("object"), " instead")
	}

	if(!missing(null) && !missing(null.RE))
		warning("argument 'null.RE' ignored if 'null' is provided")
	if(is.null(null))
		null <- null.fit(object, TRUE, null.RE, parent.frame())

	L0 <- as.vector(.getLLML(null))
	L1 <- .getLLML(object)
	n <- if(is.null(attr(L1, "nobs"))) nobs(object) else attr(L1, "nobs")
	#n <- sum(weights(object))
	ret <- 1 - exp(-2 / n * (as.vector(L1) - L0))
	max.r2 <- 1 - exp(2 / n * L0)
	attr(ret, "adj.r.squared") <- ret / max.r2
	ret
}

