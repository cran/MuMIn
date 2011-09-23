`AICc` <- function (object, ..., k = 2, REML = NULL)
UseMethod("AICc")

`AICc.logLik` <- function (object, ..., k = 2) {
	df <- attr(object, "df")
	no <- attr(object, "nobs")
	if (is.null(no)) stop("'logLik' object must have \"nobs\" attribute")
	aic <- -2 * as.numeric(object) + k * df
	aic + 2 * df * (df + 1) / (no - df - 1)
}

`AICc.default` <- function (object, ..., k = 2, REML = NULL) {
	loglik <- .getLogLik()
	
	.aicc <- function(ll, df, no)
		(-2 * ll + k * df) + (2 * df * (df + 1) / (no - df - 1))

	if (length(list(...))) {
		lls <- sapply(list(object, ...), function(x) {
			ll <- if (!is.null(REML) && inherits(x, c("mer", "lme", "gls", "lm")))
				loglik(x, REML = REML) else loglik(x)
			no <- attr(ll, "nobs")
			if (is.null(no)) no <- nobs(x, use.fallback = FALSE)
			c(as.numeric(ll), attr(ll, "df"), if (is.null(no)) NA_integer_ else no)
		})
		val <- data.frame(df = lls[2L, ], AICc= .aicc(lls[1L, ], lls[2L, ], lls[3L, ]))

		Call <- match.call()
		Call$k <- Call$REML <- NULL
		row.names(val) <- as.character(Call[-1L])
		val
	} else {
		ll <- if (!is.null(REML) && inherits(object, c("mer", "lme", "gls", "lm")))
				loglik(object, REML = REML) else loglik(object)
		no <- attr(ll, "nobs")
		if (is.null(no)) no <- nobs(object, use.fallback = FALSE)
		.aicc(as.numeric(ll), attr(ll, "df"), no)
	}
}

