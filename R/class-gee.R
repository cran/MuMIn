##=============================================================================
## Classes: gee & geeglm
##=============================================================================

`coefTable.gee` <-
`coefTable.geeglm` <-
function(model, ..., type = c("naive", "robust")) {
	cf <- summary(model, ...)$coefficients
	type <- match.arg(type)
	j <- if(type == "naive") 2L else 4L
	.makeCoefTable(cf[, 1L], cf[, j], coefNames = rownames(cf))
}


`coefTable.geese` <-
function(model, ..., type = c("naive", "robust")) {
	cf <- summary(model, ...)$mean
	type <- match.arg(type)
	j <- if(type == "naive") 2L else 4L
	.makeCoefTable(cf[, 1L], cf[, j], coefNames = rownames(cf))
}

`coef.geese` <- 
function (object, ...) object$beta


## What if 'data' changed in the meantime?
# model.matrix.gee <-
# function (object, ...) {
	# cl <- .getCall(fgee)
	# cl[[1L]] <- as.name("model.matrix")
	# cl$object <- cl$formula
	# cl$id <- cl$corstr <- cl$formula <- NULL
	# eval(cl, parent.frame())
# }


##=============================================================================
## Class: yags
##=============================================================================

`coef.yagsResult` <-
function (object, ...)
structure(object@coefficients, names = object@varnames)

`coefTable.yagsResult` <-
function(model, ..., type = c("naive", "robust")) {
	type <- match.arg(type)
	vcv <- slot(model, if(type == "naive") "naive.parmvar" else "robust.parmvar")
	.makeCoefTable(model@coefficients, sqrt(diag(vcv)), coefNames = model@varnames)
}

`getCall.yagsResult` <-
	function(x, ...) x@Call
	
	
`nobs.yagsResult` <-
function (object, ...) length(object@residuals)
	
`formula.yagsResult` <-
function (x, ...) 
eval(x@Call$formula, parent.frame())


