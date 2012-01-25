`beta.weights` <-
function(model) {
	response.sd <- sd(model.response(model.frame(model)))
	msd <- apply(model.matrix(model), 2L, sd)
	bx <- msd / response.sd

	coefmat <- coefTable(model)[, 1L:2L]
	ret <- cbind(coefmat, coefmat * bx)
	dimnames(ret) <- list(names(model$coefficients),
		c("Estimate", "Std. Error", "Beta", "Std. Err. Beta"))
	return (ret)
}
