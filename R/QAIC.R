`QAIC` <-
function(object, ..., chat) {
	#chat <- summary(gm)$dispersion
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

