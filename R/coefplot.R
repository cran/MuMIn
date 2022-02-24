
coefplot <-
function (x, lci, uci, labels = NULL,	
    width = 0.15, shift = 0,
	horizontal = TRUE,
	main = NULL, xlab = NULL, ylab = NULL,
	xlim = NULL, ylim = NULL,
	labAsExpr = TRUE, mar.adj = TRUE, lab.line = 0.5,
	lty = par("lty"), lwd = par("lwd"), pch = 21,
	col = par("col"), bg = par("bg"),
	dotcex = par("cex"), dotcol = col,
	staplelty = lty, staplelwd = lwd, staplecol = col,
    zerolty = "dotted", zerolwd = lwd, zerocol = "gray",
    las = 2,
    ann = TRUE, axes = TRUE, add = FALSE,
	type = "p",
    ...) {
	
	do.plot <- !identical(type, "n")
	horizontal <- as.logical(horizontal)[1L]
    
	y <- x
	x <- seq_len(n <- length(y))
	
    if(is.matrix(lci) && ncol(lci) == 2L) {
		uci <- lci[, 2L]
		lci <- lci[, 1L]
	}
	if(horizontal) x <- rev(x)
	
    add <- isTRUE(add)
	if(!add) {
		if(is.null(labels)) labels <- names(y)
		if(labAsExpr) labels <- .lab2expr(labels)
		mai <- par("mai")
	}

	if(...length() != 0L) {
		op2 <- par(...)
		on.exit(par(op2), add = TRUE)
	}
	
	inset <- max(.33, (width + abs(shift)) * 1.05)
	if(isTRUE(horizontal)) {

		if(!add) {
			if(mar.adj && las %in% c(1, 2)) { 
				mai[2L] <- max(mai[2L],
					grconvertX(.5 + lab.line, from = "lines", to = "inches") +
					strwidth(labels, cex = par("cex.axis"),  units = "inches"))
			}
			if(is.null(xlim)) xlim <- range(c(lci, uci))
			if(is.null(ylim)) ylim <- c(1 - inset, length(x) + inset)
		}
		zl <- c(NA, 0)
		ly0 <- ly1 <- x + shift
		lx0 <- lci
		lx1 <- uci
		wy0 <- x + shift - width
		wy1 <- x + shift + width
		wx0 <- wx1 <- c(lci, uci)
		py <- x + shift
		px <- y
		axn <- c(1L, 2L)
	} else { # vertical
		if(!add) {
			if(mar.adj && las == 2) { 
				mai[1L] <- max(mai[1L],
					grconvertY(.5 + lab.line, from = "lines", to = "inches") +
					strwidth(labels, cex = par("cex.axis"),  units = "inches"))
			}
			if(is.null(ylim)) ylim <- range(c(lci, uci))
			if(is.null(xlim)) xlim <- c(1 - inset, length(x) + inset)
		}
		zl <- c(0, NA)
		lx0 <- lx1 <- x + shift
		ly0 <- lci
		ly1 <- uci
		wx0 <- x + shift - width
		wx1 <- x + shift + width
		wy0 <- wy1 <- c(lci, uci)
		px <- x + shift
		py <- y
		axn <- c(2L, 1L)

	} # horizontal
	
	#on.exit(par(op), add = TRUE)
	if(!add) {
		#op <-
		par(mai = mai)
		plot.new()
		plot.window(xlim, ylim)
		if(!is.na(zerolty[1L]))
			abline(h = zl[1L], v = zl[2L], lty = zerolty, col = zerocol,
				lwd = zerolwd)
		title(main, xlab = xlab, ylab = ylab)
	}

	if(do.plot) {
		segments(lx0, ly0, lx1, ly1, col = col, lty = lty, lwd = lwd,
				 lend = 1L)
		if(width > 0)
			segments(wx0, wy0, wx1, wy1, col = staplecol, lty = staplelty,
				lwd = staplelwd)
		points(px, py, pch = pch, cex = dotcex, col = dotcol, bg = bg)
	}
	if(!add && isTRUE(axes)) {
		axis(axn[1L])
		axis(axn[2L], at = x, labels = labels, tick = 0, las = las, 
			 mgp = c(3, lab.line, 0))
		box()
	}
  	invisible(cbind(px, py, lx0, ly0, lx1, ly1))
}

plot.averaging <-
function(x, full = TRUE, level = 0.95,
    intercept = TRUE, parm = NULL, labels = NULL,	
    width = 0.1, shift = max(.2, width * 2.1 + .05),
	horizontal = TRUE,
	xlim = NULL, ylim = NULL,
	main = "Model-averaged coefficients",
	xlab = NULL, ylab = NULL,
	add = FALSE,
    ...) {
    
    full <- as.logical(full[1L])
	if(is.na(full)) full <- c(TRUE, FALSE)
	horizontal <- as.logical(horizontal)[1L]

	coeftypes <- c("full", "subset")[2L - full]
	if(!isTRUE(as.logical(add))) {
		lab <- substitute(lab %+-% level * `%` ~ CI,
			list(lab = sprintf("%s average ", prettyEnumStr(coeftypes)),
				level = round(level * 100, 1)))
		
		if(horizontal && missing(xlab)) xlab <- lab
		if(!horizontal && missing(ylab)) ylab <- lab
	}

	beta <- x$coefficients
	if(is.null(parm)) {
		parm <- if(isFALSE(intercept)) {
			interceptLabel <- if(!is.null(attr(x, "modelList"))) {
				unique(unlist(lapply(attr(x, "modelList"),
					function(m) attr(getAllTerms(m), "interceptLabel"))))
			} else if(!is.null(attr(x, "interceptLabel")))
				attr(x, "interceptLabel")
			which(! colnames(beta) %in% interceptLabel)
		} else seq.int(ncol(beta))
	} else if(isFALSE(intercept))
		warning("argument 'intercept' ignored, since 'parm' is given")

	beta <- beta[, parm, drop = FALSE]

	lim <- numeric(0L)
	val <- vector("list", n <- length(full))
	for(i in seq.int(n)) {
		ci <- confint(x, full = full[i], level = level, parm = parm)
		lim <- range(lim, ci, finite = TRUE)
		val[[i]] <- list(x = beta[2L - full[i], ], lci = ci)
	}

	xlim <- if(horizontal && is.null(xlim)) lim else xlim
	ylim <- if(!horizontal && is.null(ylim)) lim else ylim
	
    np <- sum(pn <- sapply(val, function(a) length(a$x)))
	fi <- rep(2 - full, pn)
	pi <- split(1L:np, rep(seq(length(pn)), pn))
    gparnm <- c("lty", "lwd", "pch", "col", "bg", "dotcex",
	  "dotcol", "staplelty", "staplelwd", "staplecol", "zerolty", "zerolwd", "zerocol")
	dots <- list(...)
	gpar <- dots[gpi <- names(dots) %in% gparnm]
	dots <- dots[!gpi]
    if(is.null(gpar$bg)) gpar$bg <- c("black", "white")[fi]
    if(is.null(gpar$pch)) gpar$pch <- 22L
	gpar <- as.data.frame(lapply(gpar, rep, length.out = np), stringsAsFactors = FALSE)
	
	off <- ((1:n) - ((n + 1) / 2)) * shift
    
	rval <- vector("list", n)
	names(rval) <- coeftypes
	
	for(i in seq.int(n)) {
		rval[[i]] <- do.call("coefplot", c(val[[i]],
			list(labels = labels,
			horizontal = horizontal,
			width = width, shift = off[i],
			add = add || i != 1L,
			xlab = xlab, ylab = ylab, main = main,
			xlim = xlim, ylim = ylim), gpar[pi[[i]], , drop = FALSE], dots))
	}
	
  	invisible(if(n == 1L) rval[[1L]] else rval)
}
