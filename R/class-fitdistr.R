
coefTable.fitdistr <-
function(model, ...) 
.makeCoefTable(coef(model), sqrt(diag(model$vcov)))

family.fitdistr <-
function (object, ...)
structure(list(family = "none", link = "NA"), class = "family")

formula.fitdistr <-
function(x, ...) 
reformulate(names(x$estimate), env = parent.frame())

nobs.fitdistr <-
function (object, ...) object$n

fitdistr2 <-
function (x, densfun, start, ...) {
    rval <- MASS::fitdistr(x, densfun, start, ...)
    densfunc <- c("beta", "cauchy", "chi-squared", "exponential", "gamma",
        "geometric", "log-normal", "lognormal", "logistic", "negative binomial",
        "normal", "Poisson", "t", "weibull")
    if(is.character(densfun)) {
        densfun <- charmatch(densfun, densfunc)
    } else {
        densfun <- deparse1(substitute(densfun))
    }
    rval$densfun <- densfun
    rval$call <- match.call()
    rval
}

