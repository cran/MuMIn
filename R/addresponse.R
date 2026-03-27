
addresponse <-
function(x, f, y = getresponse(x), env = parent.frame()) {
    stopifnot(is.language(f) || (is.numeric(f) && (f == 1 || f == -1 || f == 0))) # 0 or -1 is a valid x
    stopifnot(is.null(y) || is.language(y))
	if(length(f) == 1L || f[[1L]] != "~") {
       lhs <- f
       if(is.null(y)) {
           f <- ~ .
           f[[2L]] <- lhs
       } else {
            f <- . ~ .
            f[[2L]] <- y
            f[[3L]] <- lhs
       }
    } else {
        if(!inherits(f, "formula"))
            oldClass(f) <- "formula"
        f <- if(is.null(y)) f else {
            if(length(f) == 2L) f <- f[c(1L, NA, 2L)]
            f[[2L]] <- y
            f
        }
    }
    environment(f) <- env
    f
}

getresponse <-
function(x, f) {
    tt <- terms(formula(x))
    if(attr(tt, "response") == 0L) NULL else 
        attr(tt, "variables")[-1L][[attr(tt, "response")]]
}
