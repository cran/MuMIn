
get.response <-
function(x, data = NULL, ...)
UseMethod("get.response")

get.response.formula <-
function(x, data = NULL, ...) {
	x <- terms(x)
	if(!inherits(attr(data, "terms"), "terms")) 
		data <- model.frame(x, data = data, ...)
	data[, asChar(attr(x, "variables")[[attr(x, "response") + 1L]])]
}

get.response.lm <-
function(x, data = NULL, ...) {
	if(missing(data) && (family(x)$family != "binomial") && !is.null(x$y))
		x$y else
		#get.response.default(x, data = data, ...)
		NextMethod()
}
# NOTE: for 'binomial' 'y' is a vector not nmatrix2

get.response.averaging <-
function(x, data = NULL, ...) {
	if(is.null(attr(x, "modelList")))
		stop("'x' has no model list")
	get.response(attr(x, "modelList")[[1L]], data = data, ...)
}

get.response.default <-
function(x, data = NULL, ...) {
	if(is.null(data)) {
		# model frame:
		if(is.data.frame(x) && !is.null(tf <- attr(x, "terms"))) {
			tf <- terms(x)
			return(x[, asChar(attr(tf, "variables")[[attr(tf, "response") + 1L]])])
		} else data <- model.frame(x)
	}
	#model.frame(x)[, asChar(getResponseFormula(x))]
	if(is.null(data)) data <- model.frame(x)
	get.response(terms(x), data = data, ...)
}
