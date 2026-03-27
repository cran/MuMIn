`coefTable` <-
function (model, ...) UseMethod("coefTable")

.makeCoefTable <- 
function(x, se, df = NA_real_, coefNames = names(x)) {
	if(n <- length(x)) {
		xdefined <- !is.na(x)
		ndef <- sum(xdefined)
		if(ndef < n) {
			if(length(se) == ndef) {
				y <- rep(NA_real_, n); y[xdefined] <- se; se <- y
			}
			if(length(df) == ndef) {
				y <- rep(NA_real_, n); y[xdefined] <- df; df <- y
			}
		}
	}
	if(n && n != length(se)) stop("length(x) is not equal to length(se)")
	ret <- matrix(NA_real_, ncol = 3L, nrow = length(x),
		dimnames = list(coefNames, c("Estimate", "Std. Error", "df")))
	if(n) ret[, ] <- cbind(x, se, rep(if(is.null(df)) NA_real_ else df,
		length.out = n), deparse.level = 0L)
	class(ret) <- c("coefTable", "matrix")
	ret
}

`print.coefTable` <-
function (x, ...)
stats::printCoefmat(x[, if(all(is.na(x[, 3L]))) -3L else TRUE, drop = FALSE],
	has.Pvalue = FALSE)

summary.coefTable <-
function (object, ...) {
	tvalue <- object[, 1L] / object[, 2L]
	if (all(is.na(object[, 3L]))) {
		pvalue <- 2 * pnorm(-abs(tvalue))
		rval <- cbind(object, tvalue, pvalue)
		cn <- c("z value", "Pr(>|z|)")
	} else if (any(is.finite(tvalue))) {
		pvalue <- 2 * pt(-abs(tvalue), object[, 3L])
		cn <- c("t value", "Pr(>|t|)")
    } else {
        pvalue <- tvalue <- NaN
		cn <- c("t value", "Pr(>|t|)")
    }
	rval <- cbind(object, tvalue, pvalue)
	colnames(rval)[4L:5L] <- cn
	class(rval) <- c("summary.coefTable", class(object))
	rval
}

`print.summary.coefTable` <-
function (x, signif.stars = getOption("show.signif.stars"), ...) {
	j <- if(all(is.na(x[, 3L]))) -3L else TRUE
	stats::printCoefmat(x[, j],
		has.Pvalue = any(is.finite(x[, 5L])),
		signif.stars = signif.stars,
		...)
}

plot.mcoefTable <-
function (x, y, labAsExpr = FALSE, n = 101, w = 5, ...) {
	lab_as_expr <- function(x) {
		x <- gsub(":", "%*%", x, perl = TRUE)
		x <- gsub("\\B_?(\\d+)(?![\\w\\._])", "[\\1]", x, perl = TRUE)
		str2lang(x)
	}
	
	xd <- function(z, n, w, x = NULL) {
		rval <- matrix(NA_real_, ncol = 3L, nrow = n)
		rval[, 1L] <- if(is.null(x)) seq(z[1L] - (w * z[2L]), z[1L] + (w * z[2L]),
			length.out = n) else x
		if(!is.na(z[3L]))
			rval[, 2L] <- dt((rval[, 1L] - z[1L]) / z[2L], z[3L]) / z[2L]
		rval[, 3L] <- dnorm(rval[, 1L], z[1L], z[2L])
		rval
	}
	
	m <- nrow(x)
	lab <- if(labAsExpr) lab_as_expr(rownames(x)) else rownames(x)
	nmodels <- dim(x)[3L]
	
	col <- 1L:nmodels

	par(mfrow = n2mfrow(m), ...)
	for(i in 1L:m) {
		#cat("--", dimnames(x)[[1]][i], "--\n")

		mat <- matrix(0, n, 2L * nmodels + 1L)
		for(k in 1L:nmodels) {
			j <- seq.int(length.out = 2 + (k == 1), from = 1 + (k - 1)* 2 + (k != 1))
			xlim <- c(min(x[i, 1L, ]) - w * min(x[i, 2L, ]), max(x[i, 1L,]) + w * max(x[i,2L,]))
			vx <- seq(xlim[1L], xlim[2L], length.out = n)
			v <- xd(x[i, , k], n = n, w = w, x = vx)
			mat[, j] <- v[, (2 - (k == 1)):3]
		}
		
		plot.new()
		plot.window(xlim = range(mat[, 1L]), ylim = range(mat[,-1L], na.rm = TRUE))
		j <- seq.int(2, length.out = nmodels, by = 2)
		if(any(!is.na(mat[, j])))
			matplot(mat[, 1], mat[, j], type = "l", lty = 2, add = TRUE, col = col)
		matplot(mat[, 1L], mat[, j + 1L], type = "l", lty = 1, add = TRUE, col = col)
		abline(v = x[i, 1L, ], lty = 1L, col = col)
		abline(v = 0, lty =3L, col = 8)
		axis(1L)
		axis(2L)
		box()
		
		title(lab[[i]])
		
	}
	invisible()
}

plot.coefTable <-
function (x, y, labAsExpr = FALSE, n = 101, w = 5,
		  include.zero = TRUE,
		  col = 2, lwd = par("lwd"), lty = 1, lend = par("lend"),
		  ...) {

	xd <- function(z, n, w) {
		#rval <- matrix(NA_real_, ncol = 3L, nrow = n)
		rval <- matrix(NA_real_, ncol = 2L, nrow = n)
     	x0 <- z[1L] - w * z[2L]
	    x1 <- z[1L] + w * z[2L]
		if(include.zero) {
			x0 <- min(0, x0)
			x1 <- max(0, x1)
		}
		rval[, 1L] <- seq(x0, x1, length.out = n)
		rval[, 2L] <- if(!anyNA(z[3L])) {
				dt((rval[, 1L] - z[1L]) / z[2L], z[3L]) / z[2L]
			} else 
				dnorm(rval[, 1L], z[1L], z[2L])
		rval
	}
	
	m <- nrow(x)
	lab <- if(labAsExpr) .lab2expr(rownames(x)) else rownames(x)
	
	par(mfrow = n2mfrow(m))
	for(i in 1L:m) {	
		v <- xd(x[i, ], n = n, w = w)
		plot.new()
		plot.window(xlim = range(v[, 1L]),
					ylim = range(v[,-1L], na.rm = TRUE))
		#lines(v[, 1L], v[, 2L])
		lines(v[, 1L], v[, 2L], col = col, lwd = lwd, lty = lty, lend = lend)
		abline(v = c(x[i, 1L], 0), lty = c(1L, 3L))
		axis(1L)
		axis(2L)
		box()
		title(substitute(A == B %+-% C, list(A = lab[[i]],
			B = round(x[i, 1L], 1L), C = round(x[i, 2L], 2L))))
	}
	invisible()
}
