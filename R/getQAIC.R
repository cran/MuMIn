`getQAIC` <-
function(model, chat) {
	if (any(inherits(model,  c("lmer", "glmer")))) {
		mLogLik <- logLik(model, model@status["REML"])
		N <- NROW(model@frame)
	} else {
		mLogLik <- logLik(model)
		N <- length(resid(model))
	}

	k <- attr(mLogLik, "df") + 1
 	ret <- (deviance(model) / chat) + 2 * k
	return (ret)
}

