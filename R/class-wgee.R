coef.wgee <- function (object, ...) object$beta

coefTable.wgee <- 
function (model, ..., type = c("robust", "naive")) {
    if (match.arg(type) == "naive")
		stop("naive std. errors are not available for 'wgee' models.")
    .makeCoefTable(model$beta, summary(model)$se.robust, coefNames = names(model$beta))
}

nobs.wgee <-
function (object, ...) 
length(object$mu_fit)

formula.wgee <-
function (x, ...)
x$model