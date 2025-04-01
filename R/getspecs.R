
getspecs <-
function(x, with.model = FALSE, fitType = x@fitType) {
	if(inherits(x, "unmarkedFit")) {
		#.modelspecs[[c("unmarked", class(x)[1L], x@fitType)]]
        cls <- class(x)[1L]
        i <- match("unmarked", names(.modelspecs))
        j <- match(cls, names(.modelspecs[[i]]), nomatch = 0L)
        if(j == 0L) stop(gettextf("unknown \"unmarkedFit\" subclass: \"%s\"", cls))
        k <- match(x@fitType[1L], names(.modelspecs[[c(i, j)]]), nomatch = 0L)
        if(k == 0L) stop(gettextf("unknown \"%s\" fit type: \"%s\"", cls, x@fitType[1L]))
        specs <- .modelspecs[[c(i, j, k)]]
        estid <- paste0(names(x@estimates@estimates), ":",
                sapply(x@estimates@estimates, "slot", "short.name"))
        specs <- specs[match(estid, paste0(specs$item.name, ":", specs$short.name)), , drop = FALSE]
        stopifnot(!anyDuplicated(specs$item.name))
		specs
    } else {
       i <- match(class(x), names(.modelspecs), nomatch = 0L)
       if(all(i == 0L)) stop("unknown model class")
       .modelspecs[[i[i != 0L][1L]]]
    }
}

print.model.specs <-
function (x, ...) {
    info <- attr(x, "object.info")
    cat(sprintf("Model specification table for class '%3$s' [%1$s::%2$s()]\n",
           info["package"], info["func"], info["className"]))
    print.data.frame(x, row.names = FALSE)
}
