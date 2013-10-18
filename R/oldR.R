# compatibility with older versions of R
if(!("intercept" %in% names(formals(stats::reformulate)))) {
	`reformulate` <- function (termlabels, response = NULL, intercept = TRUE) {
		ret <- stats::reformulate(termlabels, response = response)
		if (!intercept) ret <- update.formula(ret, .~. -1)
		attr(ret, ".Environment") <- parent.frame()
		ret
	}
}

if (!exists("nobs", mode = "function", where = "package:stats", inherits = FALSE)) {
`nobs` <- function(object, ...) UseMethod("nobs")
`nobs.default` <- function(object, ...) NROW(resid(object, ...))
#`nobs.glm` <- function (object, ...) sum(!is.na(object$residuals))
}
