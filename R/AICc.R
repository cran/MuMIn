`AICc` <- 
function (object, ..., k = 2, REML = NULL)
UseMethod("AICc")

.get.nobs <- function(llik, error = FALSE) {
	no <- attr(llik, "nall")
	if(is.null(no)) no <- attr(llik, "nobs")
	if (error && is.null(no)) stop("'logLik' object must have a \"nobs\" attribute")
	no
}


`AICc.logLik` <- 
function (object, ..., k = 2) {
	if (length(list(...)) > 0L) warning("additional arguments ignored")
	df <- attr(object, "df")
	no <- .get.nobs(object, TRUE)
	aic <- -2 * as.numeric(object) + k * df
	aic + 2 * df * (df + 1) / (no - df - 1)
}



`AICc.default` <- 
function (object, ..., k = 2, REML = NULL) {
	loglik <- .getLogLik()
	
	.aicc <- function(ll, df, no)
		(-2 * ll + k * df) + (2 * df * (df + 1) / (no - df - 1))
		
	need.reml.arg <- c("mer", "merMod", "lme", "gls", "lm")

	if (length(list(...))) {
		lls <- sapply(list(object, ...), function(x) {
			ll <- if (!is.null(REML) && inherits(x, need.reml.arg))
				loglik(x, REML = REML) else loglik(x)
			no <- .get.nobs(ll, FALSE)
			if (is.null(no)) no <- nobs(x, use.fallback = FALSE)
			c(as.numeric(ll), attr(ll, "df"), if (is.null(no)) NA_integer_ else no)
		})
		val <- data.frame(df = lls[2L, ], AICc = .aicc(lls[1L, ], lls[2L, ], lls[3L, ]))

		Call <- match.call()
		Call$k <- Call$REML <- NULL
		row.names(val) <- as.character(Call[-1L])
		val
	} else {
		ll <- if (!is.null(REML) && inherits(object, need.reml.arg))
				loglik(object, REML = REML) else loglik(object)
		no <- .get.nobs(ll, FALSE)
		if (is.null(no)) no <- nobs(object, use.fallback = FALSE)
		.aicc(as.numeric(ll), attr(ll, "df"), no)
	}
}

