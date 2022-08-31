`print.model.selection` <-
function(x, abbrev.names = TRUE, warnings = getOption("warn") != -1L, ...) {
	
    origx <- x
	class(x) <- "data.frame"
	xterms <- attr(origx, "terms")  # TERMS
	
    if(is.null(xterms) || !all(xterms %in% colnames(x)[seq_along(xterms)])) {
		print.data.frame(x, ...)
	} else {
		if(abbrev.names) xterms <- abbreviateTerms(xterms, 6L, 3L, deflate = TRUE)
		colnames(x)[seq_along(xterms)] <- xterms
		globcl <- attr(origx, "global.call")
		if(!is.null(globcl)) {
			cat("Global model call: ")
			print(globcl)
			cat("---\n")
			random.terms <- attr(getAllTerms(attr(origx, "global")), "random.terms")
			if(!is.null(random.terms)) random.terms <- list(random.terms)
		} else random.terms <- attr(origx, "random.terms")
		
		dig <- c(terms = NA, varying = NA, extra = NA, df = 0L, loglik = 3L, ic = 1L, delta = 2L,
				 weight = 3L)
		column.types <- attr(origx, "column.types")
		#stopifnot(names(dig) == levels(column.types)) ## DEBUG
		
		decprint <- dig[column.types[colnames(x)]]

		i <- vapply(x, is.numeric, FALSE) & is.na(decprint)
		x[, i] <- signif(x[, i], 4L)
		k <- which(!is.na(decprint))
		for(i in k) x[, i] <- round(x[, i], digits = decprint[i])
			
		vLegend <- NULL
		if(abbrev.names) {
			vCols <- type2colname(column.types, "varying")
			vCols <- vCols[(vCols %in% colnames(x)) & !(vCols %in% c("class"))]
			if(!is.null(vCols) && length(vCols) != 0L) {
                vlen <- nchar(vCols)
                vLegend <- vector(length(vCols), mode = "list")
                names(vLegend) <- vCols
                x[, vCols] <- droplevels(x[, vCols, drop = FALSE])
				for(i in vCols) {
					if(!is.factor(x[, i])) next
					lev <- levels(x[, i])
					lev <- lev[!(lev %in% c("", "NULL"))]
					shlev <- abbreviateTerms(lev, nchar(i), deflate = TRUE)
					x[, i] <- factor(x[, i], levels = lev, labels = shlev)
					if(any(j <- shlev != lev)) vLegend[[i]] <-
						paste(shlev[j], "=", sQuote(lev[j]))
				}
				vLegend <- vLegend[!vapply(vLegend, is.null, TRUE)]
			}
		}

		uqran <- unique(unlist(random.terms, use.names = FALSE))
		abbran <- abbreviateTerms(gsub("1 | ", "", uqran, fixed = TRUE), 1L,
			deflate = TRUE)
		colran <- vapply(random.terms, function(s) paste(abbran[match(s, uqran)],
			collapse = "+"), "")

		if(addrandcol <- length(unique(colran)) > 1L) {
			k <- which(colnames(x) == "df")[1L]
			x <- cbind(x[, 1L:(k - 1L), drop = FALSE], random = colran,
                x[, k:ncol(x), drop = FALSE], deparse.level = 0L)
		}

		cat("Model selection table \n")
		if(nrow(x) == 0L) {
			print.default(colnames(x), quote = FALSE)
			cat("<0 rows>", "\n")
		} else
			print.default(as.matrix(x)[, !vapply(x, function(y) all(is.na(y)), FALSE),
			drop = FALSE], na.print = "", quote = FALSE, right = TRUE)
			
		indent <- " "

		if(abbrev.names && length(vLegend) != 0L) {
			cat("Abbreviations:", sep = "\n")
			lab <- format(paste0(indent, names(vLegend), ":"))
			for(i in seq_along(vLegend)) {
				cat(vLegend[[i]], sep = ", ", fill = TRUE, labels =
					c(lab[i], rep(paste0(rep(" ", nchar(lab[i])),
					collapse = ""), length(vLegend[[i]]) - 1L)))
			}
		}
		
		cat("Models ranked by", asChar(attr(attr(origx, 'rank'), "call")), "\n")
		if(!is.null(random.terms)) {
			if(addrandcol) {
				cat("Random terms: \n")
				cat(paste0(indent, format(abbran), ": ", uqran), sep = "\n")
			} else {
				cat("Random terms (all models): \n")
				cat(paste(uqran), sep = ", ", fill = TRUE, labels = indent)
				cat("\n")
			}
		}
		
		if (warnings && !is.null(attr(origx, "warnings"))) {
			cat("\n")
			print.warnings(attr(origx, "warnings"))
		}
	}
	invisible(origx)
}
