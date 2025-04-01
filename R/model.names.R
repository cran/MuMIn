model.names <-
function(object, ..., labels = NULL, use.letters = FALSE) {
	if (missing(object) && length(models <- list(...)) > 0L) {
		object <- models[[1L]]
	} else if (inherits(object, "list")) {
		if(length(object) ==  0L) stop("at least one model must be given")
		models <- object
		object <- models[[1L]]
	} else models <- list(object, ...)
	if(length(models) == 0L) stop("at least one model must be given")
	.modelNames(models = models, uqTerms = labels, use.letters = use.letters)
}

.modelNames <-
function(models = NULL, allTerms, uqTerms, use.letters = FALSE, ...) {
	if(missing(allTerms)) allTerms <- lapply(models, getAllTerms)
	if(missing(uqTerms) || is.null(uqTerms))
		uqTerms <- unique(unlist(allTerms, use.names = FALSE))

	n <- length(uqTerms)

	if(use.letters && n > length(LETTERS)) 
        stop("more terms than there are letters")
	sep <- if(!use.letters && n > 9L) "+" else ""

	labels <- if (use.letters) LETTERS[seq_len(n)] else as.character(seq_len(n))
	ret <- sapply(allTerms, function(x) paste(labels[sort(match(x, uqTerms))],
		collapse = sep))

	dup <- table(ret)
	dup <- dup[dup > 1L]

	if(length(dup) > 0L) {
		idup <- which(ret %in% names(dup))
		ret[idup] <- sapply(idup, function(i) paste0(ret[i],
			letters[sum(ret[seq.int(i)] == ret[i])]))
	}
	ret[!nzchar(ret)] <- "(Null)"
	attr(ret, "variables") <- structure(seq_along(uqTerms), names = uqTerms)
	ret
}
