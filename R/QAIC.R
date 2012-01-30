# AIC = -2 * LL + 2 * n
# LL = -1/2 * log(dev) * n + C
# AIC = log(dev) * n + C + 2 * n

`QAIC` <-
function(object, ..., chat, k = 2) {
	#chat <- summary(gm)$dispersion
	# chat <- deviance(object) / df.residual(object)
	logLik <- .getLogLik()

	if(chat < 1) {
		warning("'chat' given is < 1, increased to 1")
		chat <- 1
	}
	npar.adj <- if(chat == 1) 0 else 1

	qaic <- function(object, chat, k) {
		ll <- logLik(object)
		n <- attr(ll, "nobs")
		# df is the number of parameters plus 1 for estimating c-hat
		df <- attr(ll, "df") + npar.adj
		#neg2ll <- log(deviance(object)) * n # + Constant...
		neg2ll <- -2 * c(ll)
		#ret <- (neg2ll * n / chat) + k * df
		#ret <- -2 * ll / chat + k * df
		ret <- neg2ll / chat + k * df
		return (ret)
	}

	if(length(list(...))) {
		object <- list(object, ...)
		val <- data.frame(QAIC = sapply(object, qaic, chat = chat, k = k))
		Call <- match.call()
		Call$chat <- NULL
		row.names(val) <- as.character(Call[-1L])
		return(val)
	} else {
		return(qaic(object, chat, k))
	}
}


`QAICc` <-
function(object, ..., chat, k = 2) {
	logLik <- .getLogLik()
	if(chat < 1) {
		warning("'chat' given is < 1, increased to 1")
		chat <- 1
	}
	npar.adj <- if(chat == 1) 0 else 1

	qaicc <- function(object, chat, k) {
		ll <- logLik(object)
		n <- attr(ll, "nobs")
		if (is.null(n)) n <- nobs(object, use.fallback = FALSE)
		# df is the number of parameters plus 1 for estimating c-hat
		df <- attr(ll, "df") + npar.adj
		#neg2ll <- log(deviance(object)) * n # + Constant...
		neg2ll <- -2 * c(ll) # + Constant...
		ret <- (neg2ll / chat) + (k * df) * (1 + ((df + 1) / (n - df - 1)))
		return (ret)
	}

	if(length(list(...))) {
		object <- list(object, ...)
		val <- data.frame(QAIC = sapply(object, qaicc, chat = chat, k = k))
		   Call <- match.call()
		   Call$chat <- NULL
		row.names(val) <- as.character(Call[-1L])
		return(val)
	} else {
		return(qaicc(object, chat, k))
	}
}

#QAIC = -2log Likelihood/c-hat + 2K
#QAICc = -2log Likelihood/c-hat + 2K + 2K(K + 1)/(n - K - 1)
