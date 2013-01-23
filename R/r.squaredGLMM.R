
`r.squaredGLMM` <-
function(x, nullfx = NULL)
	UseMethod("r.squaredGLMM")

`r.squaredGLMM.default` <-
function(x, nullfx = NULL) 
	.NotYetImplemented()
	
`r.squaredGLMM.lme` <-
function(x, ...) {
	#VarFx  <- var(as.vector(fixef(x) %*% t(model.matrix(x, data = x$data))))
	VarFx <- var(fitted(x, level = 0L))
	vc <- vapply(lapply(x$modelStruct$reStruct, VarCorr, sigma = x$sigma),
		   "[[", double(1L), 1L)
	varAll <- sum(VarFx, vc)
	res <- c(VarFx, varAll) / (varAll + x$sigma^2)
	names(res) <- c("R2m", "R2c")
	res
}

`r.squaredGLMM.merMod` <-
`r.squaredGLMM.mer` <-
function(x, nullfx = NULL) {
	vc <- VarCorr(x)
	.rsqGLMM(x, fam = family(x),
		varFx = var(as.vector(fixef(x) %*% t(model.matrix(x)))),
		varRan = sapply(vc, "[", 1L),
		resVar = attr(vc, "sc")^2,
		fxNullCoef = fixef(if(is.null(nullfx)) 
				null.fit(x, RE.keep = TRUE, evaluate = TRUE) else nullfx
				)
		)
}

`r.squaredGLMM.glmmML` <-
function(x, nullfx = NULL) {
	if(is.null(x$x)) {
		#stop("to return model.matrix, glmmML must be fit with 'x = TRUE'")
		cl <- x$call
		cl[[1L]] <- as.name("model.matrix")
		cl$object <- cl$formula
		cl <- cl[names(cl) %in% c("", "object", "data",  "weights", "subset", "na.action")]
		X <- eval(cl)
	} else X <- x$x
	.rsqGLMM(x, family(x),
			 varFx = var(as.vector(coef(x) %*% t(X))),
			 varRan = x$sigma^2, resVar = NULL,
			 fxNullCoef = coef(if(is.null(nullfx)) 
				null.fit(x, RE.keep = TRUE, evaluate = TRUE) else nullfx
				)
			)
}


`r.squaredGLMM.lm` <-
function(x, ...)
	.rsqGLMM(x, family(x),
		 varFx = var(as.vector(coef(x) %*% t(model.matrix(x)))),
		 #varFx = var(fitted(x)),
		 varRan = 0,
		 resVar = sum(if(is.null(x$weights)) resid(x)^2 else
					   resid(x)^2 * x$weights) / df.residual(x),
		 fxNullCoef = coef(null.fit(x, evaluate = TRUE)))




`.rsqGLMM` <-
function(x, fam, varFx, varRan, resVar, fxNullCoef) {
	v <- switch(paste(fam$family, fam$link, sep = "."), 
		gaussian.identity = resVar,
		binomial.logit = 3.28986813369645, # pi^2 / 3
		binomial.probit = 1,
		poisson.log = log(1 / exp(fxNullCoef) + 1),
		poisson.sqrt = 0.25,
		stop("do not know (yet) how to calculate variance for this family/link combination")
	)
	varAll <- sum(varFx, varRan)
	res <- c(varFx, varAll) / (varAll + v)
	names(res) <- c("R2m", "R2c")
	res
}

#`model.matrix.glmmML` <- function (object, ...)  {
#	if(is.null(object$x)) stop("to return model.matrix, glmmML must be fit with 'x = TRUE'")
#	object$x
#}
