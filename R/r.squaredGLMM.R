
`r.squaredGLMM` <-
function(x)
	UseMethod("r.squaredGLMM")

`r.squaredGLMM.default` <-
function(x) 
	.NotYetImplemented()
	
`r.squaredGLMM.lme` <-
function(x) {
	VarFx <- var(fitted(x, level = 0L))

	mmRE <- model.matrix(x$modelStruct$reStruct,
						 data = x$data[rownames(x$fitted), ])
	n <- nrow(mmRE)

	sigma2 <- x$sigma^2
	varRe <- sum(sapply(x$modelStruct$reStruct, function(z) {
		sig <- pdMatrix(z) * sigma2
		mM1 <-  mmRE[, rownames(sig)]
		sum(diag(mM1 %*% sig %*% t(mM1))) / n
	}))
	varTot <- sum(VarFx, varRe)
	res <- c(VarFx, varTot) / (varTot + sigma2)
	names(res) <- c("R2m", "R2c")
	res
}


`r.squaredGLMM.merMod` <-
`r.squaredGLMM.mer` <-
function(x) {
	fam <- family(x)
	
	## for poisson(log), update 'x' to include individual-level variance (1 | 1:nobs(x)):
	if((pois.log <- fam$family == "poisson" && fam$link == "log") &&
	   !any(sapply(x@flist, nlevels) == nobs(x))) {
			cl <- getCall(x)
			frm <- formula(x)
			cl$formula <- update.formula(frm, substitute( . ~ . + (1 | gl(N,
				1)), list(N = nobs(x))))
			x <- eval(cl, envir = environment(frm), enclos = parent.frame())
		#print(getCall(x))
	}

	## can also use lme4:::subbars for that, but like this should be more
	## efficient, because the resulting matrix has only random fx variables:
	ranform <- (lapply(findbars(formula(x)), "[[", 2L))
	frm <- ranform[[1L]]
	for(a in ranform[-1L]) frm <- call("+", a, frm)
	frm <- as.formula(call("~", frm))
	#frm <- .substFun4Fun(formula(x), "|", function(e, ...) e[[2L]])
	cl <- getCall(x)
	
	mmAll <- model.matrix(frm, data = model.frame(x))
		##Note: Argument 'contrasts' can only be specified for fixed effects
		##contrasts.arg = eval(cl$contrasts, envir = environment(formula(x))))	
	
	vc <- VarCorr(x)

	n <- nrow(mmAll)
	fx <- fixef(x) # fixed effect estimates
	fxpred <- as.vector(model.matrix(x) %*% fx)
	
	if(pois.log) {
		vname <- names(x@flist)[sapply(x@flist, nlevels) == n][1L]
		varResid <-  vc[[vname]][1L]
		beta0 <- mean(fxpred)
		vc <- vc[names(vc) != vname]
	} else {
		varResid <- attr(vc, "sc")^2
		beta0 <- NULL
	}
	
	
	if(!all(c(unlist(sapply(vc, rownames))) %in% colnames(mmAll)))
		stop("random term names do not match those in model matrix. \n",
			 "Have you changed 'options(contrasts)' since the model was fitted?")
	
	varRe <- if(length(vc) == 0L) 0L else
		sum(sapply(vc, function(sig) {
			mm1 <-  mmAll[, rownames(sig)]
			sum(diag(mm1 %*% sig %*% t(mm1))) / n
		}))
	
	#varRe0 <- if(length(vc) == 0L) 0L else
	#          sum(sapply(vc, function(sig) sig[[1]]))
		
	.rsqGLMM(fam = family(x), varFx = var(fxpred), varRe = varRe,
			 varResid = varResid, beta0 = beta0)
}

`r.squaredGLMM.glmmML` <-
function(x) {
	if(is.null(x$x))
		stop("glmmML must be fitted with 'x = TRUE'")
		
	fam <- family(x)
	if(pois.log <- fam$family == "poisson" && fam$link == "log") {
		#if(length(x$posterior.modes) != nobs(x))
			.cry(NA, "cannot calculate unit-variance for poisson(log) family glmmML")
		#}
	} 
	fxpred <- as.vector(x$x %*% coef(x))
	.rsqGLMM(family(x), varFx = var(fxpred), varRe = x$sigma^2, varResid = NULL,
			 beta0 = mean(fxpred))
}




`r.squaredGLMM.lm` <-
function(x) {
	fam <- family(x)
	.rsqGLMM(fam,
		 varFx = var(as.vector(model.matrix(x) %*% coef(x))),
		 #varFx = var(fitted(x)),
		 varRe = 0,
		 varResid = sum(if(is.null(x$weights)) resid(x)^2 else
					   resid(x)^2 * x$weights) / df.residual(x),
		 beta0 = if(fam$family == "poisson" && fam$link == "log") 
			log(mean(model.response(model.frame(x)))) else 
			NULL
		 )
}



`.rsqGLMM` <-
function(fam, varFx, varRe, varResid, beta0) {
	varDistr <- switch(paste(fam$family, fam$link, sep = "."), 
		gaussian.identity = varResid,
		binomial.logit = 3.28986813369645, #  = pi^2 / 3
		binomial.probit = 1,
		poisson.log = {
			expBeta0 <- exp(beta0)
			if(expBeta0 < 6.0)
				.cry(sys.call(-1L), "exp(beta0) of %0.1f is close to zero, estimate may be unreliable \n", 
					expBeta0, warn = TRUE)
			varResid + log(1 / expBeta0 + 1)
		},
		poisson.sqrt = 0.25,
		.cry(sys.call(-1L), "do not know how to calculate variance for this family/link combination")
	) ## == Se + Sd
	varTot <- sum(varFx, varRe)
	res <- c(varFx, varTot) / (varTot + varDistr)
	names(res) <- c("R2m", "R2c")
	res
}

