
`nobs.cplm` <-
function (object, ...) 
sum(!is.na(resid(object)))

family.cplm <- 
function (object, ...) {
    getFrom("statmod", "tweedie")(var.power = object@p, link.power = object@link.power) 
    #lambda <- object@link.power
    #link <- c("log", "identity", "sqrt", "inverse")[match(lambda, c(0, 1, 0.5, -1), nomatch = 0L)]
    #if(length(link) == 0L) link <- paste0("mu^", as.character(lambda))
    #rval <- list(family = "CpPoisson", link = link)
    #class(rval) <- "family"
    #rval
}
