family.glmmadmb <-
function(object, ...) {
    famstr <- switch(object$family,
        gamma = "Gamma",
        binom = "binomial",
        nbinom = { # nbinom1 or nbinom2 ? Impossible to discern...
            s <- object$call$family
            if(!is.character(s)) {
                warning("\"nbinom\" family may be incorrectly determined if ", sQuote(deparse(s)),
                    " has changed since model fitting")
                s <- eval(s, environment(formula(object)))
                if(!is.character(s) || length(s) != 1L || !tolower(s) %in% c("nbinom", "nbinom2", "nbinom1"))
                    stop("cannot determine \"nbinom\" family type. Use a character string literal in the model call.")
                if(s == "nbinom") s <- "nbinom2"
            }
            s
        }, 
        object$family)
    
    if(famstr %in% c("poisson", "binomial", "gaussian", "Gamma"))
        return(do.call(match.fun(famstr), list(object$link)))
    
    rval <- list(family = famstr, link = object$link, linkfun = object$linkfun, linkinv = object$ilinkfun)
    class(rval) <- "family"
    rval
}

# residual std. dev. or dispersion parameter
sigma.glmmadmb <-
function (object, ...) {
    switch(family(object)$family,
        Gamma = 1 / sqrt(object$alpha),
        object$alpha
       )
    # XXX: not checked with other non-standard families
}
