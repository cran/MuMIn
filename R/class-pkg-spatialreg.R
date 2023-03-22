nobs.Spautolm <-
function (object, ...)
    object$fit$N

nobs.Sarlm <-
function (object, ...)
   length(resid(object))

coefTable.Sarlm <-
coefTable.Spautolm <-
function (model, ...) {
    cf <- summary(model)$Coef
    .makeCoefTable(cf[, 1L], cf[, 2L])
}

