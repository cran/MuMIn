`importance` <-
function (x) UseMethod("importance")

`importance.averaging` <-
function(x) return(x$importance)

`importance.model.selection` <-
function(x) {
	tt <- attr(x, "terms")
	z <- x[, tt, drop = FALSE]
	z <- z[, !apply(apply(z, 2L, is.na), 2, all) & !(tt %in% attr(tt, "interceptLabel")),
		drop = FALSE]
	wt <- x[, "weight"]
	return(sort(apply(z, 2L, function(y) sum(wt[!is.na(y)])), decreasing = TRUE))
}




#function(x) return(apply(x[, attr(x, "terms")], 2L,
#	function(z) sum(x[, "weight"][!is.na(z)])))

`importance.default` <- function(x)
	model.avg(x)$importance