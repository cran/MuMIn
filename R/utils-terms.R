
.multipleformula2list <- 
function(x) {
    f <- formula(x)
    env <- environment(f)
    rval <- list()
    while(is.call(f) && f[[1L]] == "~") {
        rval[[length(rval) + 1L]] <- as.formula(f[c(1L, length(f))], env = env)
        f <- f[[2L]]
    }
    rval <- rev(rval)
    if(!identical(f, rval[[1L]][[2L]])) {
        # has lhs - attach the response to each right-sided formula:
        rval <- lapply(rval, function(g) {
            g <- g[c(1L, NA, 2L)]
            g[[2L]] <- f
            g
        })
    }
    lapply(rval, `environment<-`, env)
}

