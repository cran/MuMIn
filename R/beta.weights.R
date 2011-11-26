`beta.weights` <-
function(model) {
	response.sd <- sd(model.response(model.frame(model)))
	msd <- apply(model.matrix(model), 2L, sd)
	bx <- msd / response.sd

	coefmat <- coefTable(model)
	ret <- cbind(coefmat[, 1L:2L], coefmat[, 1L:2L] * bx)
	dimnames(ret) <- list(names(model$coefficients),
		c("Estimate", "Std. Err.", "Beta", "Std. Err. Beta"))
	return (ret)
}
