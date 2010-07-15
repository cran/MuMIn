`QAIC` <-
function(object, ..., chat) {
	#chat <- summary(gm)$dispersion

	`getQAIC` <- function(object, chat) {
		mLogLik <- logLik(object)
		N <- length(residuals(object))
		k <- attr(mLogLik, "df") + 1
		ret <- (deviance(object) / chat) + 2 * k
		return (ret)
	}

	if(length(list(...))) {
		object <- list(object, ...)
		val <- data.frame(QAIC=sapply(object, getQAIC, chat = chat))
		   Call <- match.call()
		   Call$chat <- NULL
		row.names(val) <- as.character(Call[-1])
		return(val)
	} else {
		return(getQAIC(object, chat = chat))
	}
}
