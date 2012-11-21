# Calculate Akaike weights
`Weights` <-
function(x)  UseMethod("Weights")

`Weights.model.selection` <-
function(x) x[, "weight"] / sum(x[, "weight"])

`Weights.averaging` <-
function(x) x$summary$Weight

`Weights.data.frame` <-
function(x) {
	if(ncol(x) == 2L && colnames(x)[2L] %in% c("AIC", "AICc", "BIC", "QAIC", "QAICc")
	&& is.numeric(x[, 2L]))
		Weights.default(x[, 2L])
	else NA
}

`Weights.default` <-
function(x) {
	delta <- x - min(x)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))
	return (weight)
}