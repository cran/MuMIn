`importance` <-
function (x) UseMethod("importance")

`importance.averaging` <-
function(x) return(x$importance)

`importance.model.selection` <-
function(x) return(apply(x[, attr(x, "terms")], 2L,
	function(z) sum(x[, "weight"][!is.na(z)])))

`importance.default` <- function(x)
	model.avg(x)$importance