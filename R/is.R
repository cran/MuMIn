.isREMLFit <-
isREML <- 
function (x) {
    if (inherits(x, "merMod")) 
        return(lme4::isREML(x))
    if (inherits(x, c("lme", "gls", "gam")) && is.character(x$method)) 
        return(x$method[1L] %in% c("lme.REML", "REML"))
    return(FALSE)
}


isGEE <- 
function(object) 
inherits(object, c("geeglm", "geese", "gee", "geem", "wgee", "yagsResult"))
