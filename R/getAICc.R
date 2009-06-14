`getAICc` <-
function(model, k = 2) {
	if (any(inherits(model,  c("lmer", "glmer")))) {
		mLogLik <- logLik(model, model@status["REML"])
		N <- NROW(model@frame)
	} else {
		mLogLik <- logLik(model)
		N <- length(resid(model))
	}

	mK <- attr(mLogLik, "df")
	mAIC <- -2 * c(mLogLik) + k * mK
	ret <- mAIC + 2 * mK * (mK + 1)/(N - mK - 1)
	attr(ret, "df") <- mK
	attr(ret, "AIC") <- mAIC
	return (ret)
}

