`beta.weights` <-
function(model) {

	ttbl <- tTable(model)
	m.coef <- ttbl[,1]
	m.se <- ttbl[,2]

	response.sd <- sd(eval(attr(model$terms, "variables"),
		envir=model$model)[[attr(model$terms, "response")]])
	m.terms.sd <- sd(model.matrix(model))
	bx <- m.terms.sd / response.sd

	m.b.coef <- m.coef * bx
	m.b.se <- m.se * bx

	ret <- data.frame(m.coef, m.se, m.b.coef, m.b.se)

	colnames(ret) <- c("Estimate", "Std. Err.", "Beta", "Std. Err. Beta")
	rownames(ret) <- names(model$coefficients)

	ret <- as.matrix(ret)
	return (ret)
}
