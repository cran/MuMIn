`rbind.model.selection` <-
function (..., deparse.level = 1, make.row.names = TRUE) {
	
  	n <- ...length()
	if(n == 1L) return(...elt(1L))
    
    items <- list(...)

	if(!all(vapply(items, inherits, FALSE, "model.selection")))
		stop("not all objects are \"model.selection\"")

	### XXX: This modifies original objects!!!
	items <- lapply(items, "class<-", "data.frame")
	## ... reverting to original (?) class on exit:
	on.exit(lapply(items, "class<-", c("model.selection", "data.frame")))

	.allitemsidentical <- 
    function(x) all(vapply(x[-1L], identical, FALSE, x[[1L]]))

	if(!.allitemsidentical(lapply(lapply(items, attr, "rank"), attr, "call")))
		stop("tables are not ranked by the same IC")
	if(!.allitemsidentical(lapply(items, "attr", "nobs")))
		stop("models are fitted to different number of observations")

	.combine <-
	function(x, y, pos, len = length(y)) {
		if(is.factor(x) || is.factor(y)) {
			if(is.factor(x)) {
				if(!is.factor(y)) y <- factor(y)
			} else if(is.factor(y)) x <- factor(x)
			alllev <- unique(c(levels(x), levels(y)))
			x <- factor(x, levels = alllev, labels = alllev)
		}
		x[pos:(pos + len - 1L)] <- y
		x
	}

	ct <- unname(lapply(items, attr, "column.types"))
	vct <- unlist(ct, recursive = FALSE)
	vct <- vct[order(as.integer(unlist(ct)))]

	#vct <- vct[order(as.integer(unlist(ct)), unlist(lapply(ct, seq_along)))]
	vct <- vct[!duplicated(names(vct))]
	# TODO: check mismatch in column.types
	nm <- names(vct)

	rval <- as.data.frame(array(NA, dim = c(sum(sapply(items, nrow)), length(nm)),
								dimnames = list(NULL, nm)))
	row1 <- 1L
	for(z in items) {
		n <- nrow(z)
		nmz <- nm[nm %in% names(z)]
		for(j in nmz) 
            rval[, j] <- .combine(rval[, j], z[, j], row1, n)
		row1 <- row1 + n
	}

    combineAttrs <- c("model.calls", "coefTables")

    hasModelList <- 
        vapply(items, function(x) is.list(attr(x, "modelList")), logical(1L))
    if(any(hasModelList)) {
        if(!all(hasModelList)) {
            warning("not all combined tables include model lists. The missing items will be NULL.")
            for(i in which(!hasModelList))
                attr(items[[i]], "modelList") <-
                    structure(vector("list", nrow(items[[i]])),
                        names = rownames(items[[i]]))
        }
        combineAttrs[length(combineAttrs) + 1L] <- "modelList"
    }

	newattr <- list(column.types = vct)
	for(i in combineAttrs)
		newattr[[i]] <- unlist(lapply(items, attr, i), recursive = FALSE, use.names = FALSE)
	k <- c("rank", "nobs")
	newattr[k] <- attributes(items[[1L]])[k]

	tmp <- lapply(items, attr, "terms")
	newattr[["terms"]] <- 
        structure(unique(unlist(tmp, recursive = FALSE, use.names = FALSE)),
		  interceptLabel = unique(unlist(lapply(tmp, attr, "interceptLabel"))))


	for(i in names(newattr)) 
        attr(rval, i) <- newattr[[i]]
	class(rval) <- c("model.selection", "data.frame")
	if(make.row.names) {
		rn1 <- rep(names(items), sapply(items, nrow))
		rn1[i] <- paste0(rn1[i <- rn1 != ""], ".")
		rlabs <- paste0(rn1, unlist(lapply(items, rownames)))
		if(anyDuplicated(rlabs))
			rlabs <- make.unique(as.character(rlabs), sep = "")
	} else {
		rlabs <- as.character(1L:nrow(rval))
	}
	rownames(rval) <- rlabs

	o <- order(rval[, names(vct)[vct == "ic"]])
	rval <- rval[o, recalc.delta = TRUE]
	attr(rval, "merged-order") <- o
	rval
}

`merge.model.selection` <-
function (x, y, suffixes = c(".x", ".y"), ...) {
	rval <- rbind(x, y, make.row.names = FALSE)
	if (!is.null(suffixes)) row.names(rval) <-
		c(paste0(row.names(x), suffixes[1L]),
            paste0(row.names(y), suffixes[2L]))[attr(rval, "merged-order")]
	attr(rval, "merged-order") <- NULL
	rval
}
