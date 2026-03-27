
formula.formlist <-
function(x, env = parent.frame(),
         style = c("~", "|", "+", "tilde", "bar", "combined"),
         response = NULL,
         ...) {
    
    style <- if(inherits(x, "formlist") && !is.null(attr(x, "mform.style")))
        attr(x, "mform.style") else match.arg(style)
		
    rhslist <- lapply(x, "[[", 2L)
    
    fexpr <- switch(style, "|" =, bar = {
        #"bar separated"
        ff <- Reduce(\(a,b) call("|", a, b), rhslist)
		addresponse(f = ff, y = response, env = env)
    }, "~" = , tilde =, {
        #"tilde separated"
        rhslist[[1L]] <- addresponse(f = rhslist[[1L]], y = response)
        Reduce(\(a, b) call("~", a, b), rhslist)
    }, "+" =, combined = {
        #"all model formulas combined togeter"
        ff <- str2lang(paste(vapply(rhslist, deparse1, ""), collapse = "+"))
        f <- addresponse(f = ff, y = response, env = env)
    })
    formula(fexpr, env)
}

# IF tilde-separated formulas, response is only for the first model
# IF bar-separated formulas, response is with all component models
formlist <- 
function(x, env = environment(f)) {
    f <- formula(x)
    
    .op <- \(expr) if(is.call(expr)) expr[[1L]] else ""

    # decompose multi-sided (n > 2), i.e. [y] ~ a ~ b ~ c
    f1 <- f
    rval <- list()
    while(.op(f1) == "~") {
        rval[[length(rval) + 1L]] <- f1[[length(f1)]]
        f1 <- f1[[2L]]
    }
    response <- if(!identical(f1, rval[[length(rval)]])) f1 else NULL
    
    mform.tilde <- length(rval) != 1L
    
    # one or two-sided [y] ~ a + b
    if(length(rval) == 1L) {
        rhs <- rval[[1L]] # remove ~
        rval <- list()
        # bar-split formula y ~ a | b
        while(.op(rhs) == "|") {
            rval[[length(rval) + 1L]] <- rhs[[3L]]
            rhs <- rhs[[2L]]
        }
        rval[[length(rval) + 1L]] <-  rhs
    }
    
    rval <- rev(rval)
    if(is.null(response)) {
        noresp <- seq.int(1L, length(rval))
        resp <- 0L
    } else {
        if(mform.tilde) {
            # has lhs and is ~-separated - attach the response to FIRST right-sided formula:
            resp <- 1L
            noresp <- seq.int(2L, length(rval)) # because length(rval) > 1L
        } else {
            # has lhs - attach the response to each right-sided formula:
            resp <- seq.int(1L, length(rval))
            noresp <- 0L
        }
    }
    rval[resp] <- lapply(rval[resp], \(g) formula(call("~", response, g), env = env))
    rval[noresp] <- lapply(rval[noresp], \(g) formula(call("~", g), env = env))
    
    if(length(rval) != 1L)
        attr(rval, "mform.style") <- if(mform.tilde) "~" else "|"
    
    oldClass(rval) <- "formlist"
    rval
}
