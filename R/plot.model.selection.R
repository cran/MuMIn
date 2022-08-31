
`plot.model.selection` <-
function(x,
	ylab = NULL, xlab = NULL, main = "Model selection table",
	labels = NULL, terms = NULL, 
    labAsExpr = TRUE,
	vlabels = rownames(x),
	mar.adj = TRUE,
	col = NULL,
	col.mode = 2,
    bg = "white",
	border = par("col"),
	par.lab = NULL, par.vlab = NULL,
	axes = TRUE, ann = TRUE, ...) {

    if(is.null(xlab)) xlab <- NA  
	if(is.null(ylab)) ylab <- expression("Cumulative model weight" ~~ (omega))

	vlabels.wts.cutoff <- 0.01
	
	if(is.numeric(col.mode)) {
		if(col.interpolate <- col.mode[1L] > 0) {
			col.interp.bias <- col.mode[1L]
			col.by.value <- FALSE
		} else if(col.mode[1L] < 0) {
			col.by.value <- TRUE
		} else
			col.by.value <- FALSE
	} else if (is.character(col.mode) && startsWith("value", tolower(col.mode[1L]))) {
		col.by.value <- TRUE
		col.interpolate <- FALSE
	} else {
		col.by.value <- col.interpolate <- FALSE
	}

	if(...length() != 0L) {
        dots <- list(...)
        if(!is.null(dots$col2)) {
            warning("argument 'col2' is now defunct")
            dots$col2 <- NULL
        }
		op <- do.call("par", c(dots, no.readonly = TRUE))
		on.exit(par(op))
	}
    
    wts <- Weights(x)
    ok <- wts > 1e-5
	
	wts <- wts[ok]
    cwts <- cumsum(wts)     
    #xp <- !is.na(itemByType(x, "terms", drop = FALSE))
    xp <- !is.na(itemByType(x, c("terms", "varying"), drop = FALSE))
	
	if(is.null(terms)) terms <- TRUE
	
    xp <- xp[ok, terms, drop = FALSE]
    m <- ncol(xp)
    n <- nrow(xp)

	if (isTRUE(col.by.value)) {
        if(is.null(col)) {
			col1 <- hcl.colors(10L, palette = "Blue-Red 3")
			col2 <- "gray50"
		} else {
			if(length(col) < 2L)
				stop("colours by value need 'col' with at least two elements")
			col2 <- col[1L]
			col1 <- col[-1L]
		}
    
        ncq <- length(col1)
        ncqscale <- (ncq - 1) / 2
   
        cft <- do.call("cbind", lapply(x[, terms], function(x) {
            if(all(is.na(x))) {
                as.integer(x)
            } else if(is.character(x)) { # in case stringsAsFactors was not used in dredge
                as.numeric(factor(x)) + ncq
            } else if(is.factor(x)) {
                as.numeric(x) + ncq
            } else {
                m <- max(abs(x), na.rm = TRUE)
                floor((x / m * ncqscale) + ncqscale + 1L)
            }
        }))

        mode(cft) <- "integer"
        col <- array(c(col1, col2[1L])[cft], dim = dim(x))
        
        col.interp.bias <- 0L
        col.interpolate <- FALSE
    } else if(is.null(col)) {
		cola <- grDevices::hcl.colors(25L, palette = "Blues 3", rev = TRUE)
		colb <- grDevices::rgb(desaturate(t(col2rgb(cola)), .66), maxColorValue = 255)
		col <- cbind(cola, colb, deparse.level = 0L)
	}

		
	if(isTRUE(col.interpolate)) {
		if(!is.matrix(col)) {
			col <- matrix(col)
		} else if(m < ncol(col))
			col <- col[, seq.int(min(m, ncol(col))), drop = FALSE]
	
		colwts <- wts / max(wts)
		if(nrow(col) < 3L) col <-
			apply(col, 2L, function(a) grDevices::colorRampPalette(a)(3L))
		col1 <- array("", dim = c(length(colwts), ncol(col)))
		for(i in 1L:ncol(col)) col1[, i] <-
			rgb(grDevices::colorRamp((col[, i]),
				bias = col.interp.bias)(colwts), maxColorValue = 255)
		col <- col1
	} else { #recycle colors row- and column-wise
		if(is.matrix(col))
			col <- col[rep(seq.int(nrow(col)), length.out = n),
				   rep(seq.int(ncol(col)), length.out = m)]
	}

	if(is.matrix(col)) {
		a <- rep(seq(0, by = n, length.out = ncol(col)),
			 each = n, length.out = length(xp))
		x2 <- array(rep(seq.int(n), m), dim = dim(xp)) * xp
		x2[!xp] <- NA
		#x2 <- unname(x2)
	} else {
		a <- 0L
		x2 <- array(NA_integer_, dim = dim(xp))
		x2[xp] <- rep(seq.int(length(col)), length.out = sum(xp))		
		#x2[] <- (x2 - 1L) %% length(col) + 1L
	}
	
    plot.new()	# need it here for reading 'par'
    
	if(isTRUE(ann)) {
		commonpar <- list(col = par("col.axis"), font = par("font.axis"),
			cex = par("cex.axis") * par("cex"))
		
		if(missing(labels)) labels <- colnames(xp)
		
		if(!is.null(labels)) {
			if(length(labels) != m)
				stop("length of 'labels' is not equal to number of terms")
			if(labAsExpr && is.character(labels))
				labels <- .lab2expr(labels)
				
			arglab <- c(list(las = 2L, line = 0.33, padj = 0.5), commonpar)
			for(i in names(par.lab)) arglab[i] <- par.lab[i]
				
			vlabels <- vlabels[ok][vli <- wts > vlabels.wts.cutoff]
		
			argvlab <- c(list(las = 2L, mgp = c(1, .5, 0), hadj = 0), commonpar)
			for(i in names(par.vlab)) argvlab[i] <- par.vlab[i]
			if(!is.numeric(argvlab$mgp)) argvlab$mgp <- par("mgp")
			if(is.numeric(argvlab$line)) {
				argvlab$mgp[2L] <- argvlab$line
				argvlab$line <- NULL
			}
			# for rhs-axis, replace 'axpar' with 'axpar.axis'
			axp <- c("col", "font", "cex")
			j <- match(axp, names(argvlab), nomatch = 0L)
			names(argvlab)[j] <- paste0(names(argvlab)[j], ".axis")
			
			if(isTRUE(mar.adj)) {
				# top labels:
				sw <- max((if(arglab$las == 1) strheight else strwidth)(
						labels, font = arglab$font,
						cex = arglab$cex / par("cex"), units = "in"))
				mai <- par("mai")
				mai[3L] <- max(mai[3L], sw + grconvertY(arglab$line + .33 +
					(if(is.null(main)) 0 else 2),
					"lines", "inches"))
				
				# right-hand side labels:
				ml <- argvlab$mgp[2L]
				if(argvlab$las == 3) {
					sw <- 0
					ml <- ml + 1
				} else {
					sw <- (1 - argvlab$hadj) * 
						max(strwidth(vlabels, font = argvlab$font,
						cex = argvlab$cex, units = "in"))
					ml <- ml + .25
				}
				mai[4L] <- max(mai[4L], sw + grconvertX(ml, "lines", "inches"))
				op2 <- par(mai = mai)
				on.exit(par(op2), add = TRUE)
			} # mar.adj
		} # labels
	} # ann
    
	plot.window(xlim = c(0.5, m + .5), ylim = c(1, 0), xaxs = "i", yaxs = "i")
    rect(0.5, 0, m + .5, 1, col = bg, border = 0)

    #image(seq.int(ncol(x2)), c(0, cwts), t(a + x2), col = col, add = TRUE)
	cx <- seq(0.5, ncol(x2) + .5)
	cy <- c(0, cwts)
	
	ixy <- expand.grid(y = seq.int(length(cwts)), x = seq.int(ncol(x2)))
	#j <- !is.na(x2) ### XXX === xp!
	#koBrowseHere()
	ixy <- ixy[xp, , drop = FALSE]
	plot.window(xlim = c(0.5, m + .5), ylim = c(1, 0), xaxs = "i", yaxs = "i")
	rect(cx[ixy[, 2L]], cy[ixy[, 1L]], cx[ixy[, 2L] + 1L], cy[ixy[, 1L] + 1L],
		 col = col[c(x2 + a)][xp])
    abline(h = cwts, v = seq(.5, length.out = m))
    box()

   	if(isTRUE(ann)) {
		do.call(mtext, c(list(labels, at = seq.int(m), side = 3L), arglab))
		do.call(axis, c(list(side = 4L, at = cwts[vli] - 0.5 * wts[vli],
			tick = FALSE, labels = vlabels), argvlab))
		title(main = main, line = par("mar")[3L] - 1.5)
		title(ylab = ylab)
		title(xlab = xlab, line = 2)
	}
	if(axes) {
		axis(2L, col = border, col.ticks = border)
		box(col = border)
	}
	invisible(NULL)
}
