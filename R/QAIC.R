`QAIC` <-
function(object, ..., chat) {
	#chat <- summary(gm)$dispersion
	ll <- .getLogLik()

	.getQAIC <- function(object, chat) {
		mLogLik <- ll(object)
		N <- attr(mLogLik, "nobs")
		k <- attr(mLogLik, "df") + 1
		ret <- (deviance(object) / chat) + 2 * k
		return (ret)
	}

	if(length(list(...))) {
		object <- list(object, ...)
		val <- data.frame(QAIC=sapply(object, .getQAIC, chat = chat))
		   Call <- match.call()
		   Call$chat <- NULL
		row.names(val) <- as.character(Call[-1])
		return(val)
	} else {
		return(.getQAIC(object, chat = chat))
	}
}


`QAICc` <-
function(object, ..., chat) {
	#chat <- summary(gm)$dispersion
	ll <- .getLogLik()

	.getQAICc <- function(object, chat) {
		mLogLik <- ll(object)
		N <- attr(mLogLik, "nobs")
		k <- attr(mLogLik, "df") + 1
		ret <- (deviance(object) / chat) + (2 * k) * (1 + ((k + 1) / (N - k - 1)))
		return (ret)
	}

	if(length(list(...))) {
		object <- list(object, ...)
		val <- data.frame(QAIC=sapply(object, .getQAICc, chat = chat))
		   Call <- match.call()
		   Call$chat <- NULL
		row.names(val) <- as.character(Call[-1])
		return(val)
	} else {
		return(.getQAICc(object, chat = chat))
	}
}

#QAIC = -2log Likelihood/c-hat + 2K
#QAICc = -2log Likelihood/c-hat + 2K + 2K(K + 1)/(n - K - 1)
