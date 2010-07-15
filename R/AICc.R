`AICc` <-
function(object, ..., k = 2) {
	if(length(list(...))) {
		object <- list(object, ...)
		val <- as.data.frame(t(sapply(object,
 							function(el) {
					     		z <- getAICc(el, k = k)
					     		c(attr(z, "df"), attr(z, "AIC"), z)
						  	}
						  )))
		names(val) <- c("df", "AIC", "AICc")
		   Call <- match.call()
		   Call$k <- NULL
		row.names(val) <- as.character(Call[-1])
		return(val)
	} else {
		return(getAICc(object, k = k))
	}
}

`getAICc` <-
function(object, k = 2) {
	# No longer needed?
	#if (any(inherits(object, c("lmer", "glmer")))) {
	#	mLogLik <- logLik(object, object@status["REML"])
	#	N <- NROW(object@frame)
	#} else {	}

	mLogLik <- logLik(object)
	N <- length(residuals(object))

	mK <- attr(mLogLik, "df")
	mAIC <- -2 * c(mLogLik) + k * mK
	ret <- mAIC + 2 * mK * (mK + 1)/(N - mK - 1)
	attr(ret, "df") <- mK
	attr(ret, "AIC") <- mAIC
	return (ret)
}
