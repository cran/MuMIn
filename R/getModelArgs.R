getModelArgs <-
function() {
	.makemnames <- function(cl) {
		if(!is.null(names(cl))) {
			argnames <- names(formals(sys.function(-2L)))
			cl <- cl[! names(cl) %in% argnames[-c(1L, 2L)]]
			names(cl) <- names(cl)[names(cl) != argnames[1L]]
		}
		unlist(.makeListNames(cl[-1L]))
	}
	
	cl <- sys.call(-1L)
	
	pf <- parent.frame()
	missingObject <- evalq(missing(object), pf)
	if(missingObject) NULL else object <- get("object", pf, inherits = FALSE)
	models <- eval(call("list", as.name("...")), pf)

	
	if (missing(object) && length(models) > 0L) {
		object <- models[[1L]]
		names(models) <- .makemnames(cl)
	} else if (is.list(object) && !is.object(object)) {
		if(length(object) ==  0L) stop("at least one model must be given")
		models <- object
		object <- models[[1L]]
		names(models) <- unlist(.makeListNames(models))
	} else if (inherits(object, "averaging")) {
		modelList <- attr(object, "modelList")
		if(length(models)) cry(-1, "extra arguments ignored", warn = TRUE)
		if(is.list(modelList) && length(modelList)) {
			models <- modelList
		} else  cry(-1, "'object' is an \"averaging\" object but has no model list",
			warn = FALSE)
		
	} else {
		models <- c(list(object), models)
		if(length(models) > 1L) {
			names(models) <- .makemnames(cl)
		} else {
			names(models)[1L] <- unlist(.makeListNames(list(substitute(object))))
		}
	}
	if(length(models) == 0L) stop("at least one model must be given")
	invisible(models)
}

# TODO:
#is.listOfCalls <-
#function(x)  is.list(x) && all(vapply(x, is.call, FALSE))

