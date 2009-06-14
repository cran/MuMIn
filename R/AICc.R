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

