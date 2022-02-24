## Helper functions:

# RE model matrix
.remodmat <- function(object) 
UseMethod(".remodmat")

#.remodmat.default <- function(object) model.matrix(.ranform(formula(object)), data = model.frame(object))
.remodmat.default <-
function(object) {
    env <- environment(formula(object))
    rval <- lapply(.findbars(formula(object)),
        function(f) model.matrix(as.formula(call("~", f[[2L]]), env = env),
                     data = model.frame(object)))
    rval <- do.call("cbind", rval)
    rval[, !duplicated(colnames(rval)), drop = FALSE]
}
    
.remodmat.merMod <-
function(object) {
    rval <- do.call("cbind", model.matrix(object, type = "randomListRaw"))
	rval[, !duplicated(colnames(rval)), drop = FALSE]
}

.remodmat.lme <-
function(object)
    model.matrix(object$modelStruct$reStruct, data = object$data[rownames(object$fitted), 
			, drop = FALSE])


.nullUpdateWarning <- 
function(message = 
"the null model is correct only if all variables used by the original model remain unchanged.",
Call = NULL) {
	if(!isTRUE(getOption("MuMIn.noUpdateWarning")))
		cry(Call, message, warn = TRUE)
}

# .nullFitRE: update `object` to intercept only model, keeping original RE terms.
# TODO: reOnlyModelCall or reOnlyFormula
.nullFitRE <- function(object, envir) UseMethod(".nullFitRE")
.nullFitRE.default <- 
function(object, envir = parent.frame()) {
    cl <- get_call(object)
	if(! "formula" %in% names(cl)) 
		stop("cannot create a null model for an object without named \"formula\" argument in its call. ")
        
    if(any(grepl("^..\\d$", all.names(cl))))
        stop("object's call contains dotted names: ", sQuote(deparse(cl, control = NULL)),
             "and cannot be evaluated. See '?updateable' for a workaround.")
 
    cl$formula <- .nullREForm(formula(object))
	.nullUpdateWarning()
    eval(cl, envir)
}

.nullFitRE.lme <-
function(object, envir = parent.frame()) {
	cl <- getCall(object)
	cl$fixed <- update.formula(cl$fixed, . ~ 1)
    if(inherits(object, "glmmPQL")) cl$verbose <- FALSE
	.nullUpdateWarning()
	eval(cl, envir)
}

# sum up RE variance using VarCorr list
# For RE-intercept identical to a sum of diagonals of VC matrices.
# sum(sapply(lapply(vc, diag), sum))
.varRESum <-
function(vc, X) {
    if(is.null(vc)) return(0)
	n <- nrow(X)
	sum(sapply(vc, function(sig) {
		mm1 <-  X[, rownames(sig), drop = FALSE]
		sum(matmultdiag(mm1 %*% sig, ty = mm1)) / n
	}))
}

## extracts random effect formula. e.g:
.ranform <-
function (form) {
	### XXX: would give an error: values must be length 1 ...
	###      for very long RE formulas
	ans <- reformulate(vapply(lapply(.findbars(form),
		"[[", 2L), deparse, "", width.cutoff = 500L))
    # XXX: Why?
    #update.formula( , ~ . + 1)
	environment(ans) <- environment(form)
	ans
}

# update model adding an observation level RE term
.OLREFit <- function(object) UseMethod(".OLREFit")
.OLREFit.default <- function(object) .NotYetImplemented()
.OLREFit.merMod <- function(object) {
    if (!any(sapply(object@flist, nlevels) == nobs(object))) {
		cl <- get_call(object)
        frm <- formula(object)
		nRowData <- eval(call("eval", as.expression(call("NROW", cl$formula[[2L]])),
            envir = cl$data), envir = environment(frm),
            enclos = parent.frame())
		fl <- length(frm)
		frx <- . ~ . + 1
		frx[[3L]][[3L]] <- call("(", call("|", 1, call("gl", nRowData, 1)))
		cl$formula <- update.formula(frm, frx)		
		object <- tryCatch(eval(cl, envir = environment(frm), enclos = parent.frame()),
			error = function(e) {
				cry(conditionCall(e), conditionMessage(e), warn = TRUE)
				cry(cl, "fitting model with the observation-level random effect term failed. Add the term manually")
		})
		.nullUpdateWarning("the result is correct only if all variables used by the model remain unchanged.")
    }
    object
}



.binomial.sample.size <-
function(object) {
    tt <- terms(formula(object))
    y <- model.frame(object)[, rownames(attr(tt, "factors"))[attr(tt, "response")]]
    if(is.null(dim(y))) 1 else mean(rowSums(y))        
}
