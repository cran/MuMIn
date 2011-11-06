# This is merely to get rid of the annoying behaviour in summary.glmML.
# it does not do anything except for printing the model output.
`summary.glmmML` <- function(object, ...) object



# family
`family.default` <- function (object, ...)  {
	cl <- getElement(object, "call")
	if(is.null(cl)) return(NULL)
	fam <- cl$family
	if(is.null(fam)) fam <- formals(match.fun(cl[[1L]]))$family
	if(is.null(fam)) return(gaussian())
	switch(mode(fam), call=eval(fam), name =, character = match.fun(fam)())
}


`family.gls` <-
`family.lme` <-
stats:::family.lm



`nobs.rq` <-
function (object, ...) length(object$y)
