\name{exprApply}
\alias{exprApply}
\encoding{utf-8}
\title{Apply a function to calls inside an expression}
\description{
Apply function \code{FUN} to each occurence of a call to \code{what()} (or 
a symbol \code{what}) in an unevaluated expression. It can be used for advanced 
manipulation of expressions.
Intended primarily for internal use.
}

\usage{
exprApply(expr, what, FUN, ..., symbols = FALSE)
}

\arguments{
  \item{expr}{an unevaluated expression. }
  \item{what}{character string giving the name of a function. Each call to
      \code{what} inside \code{expr} will be passed to \code{FUN}. \code{what} 
      can be also a character representation of an operator or parenthesis
      (including \link[=Paren]{curly} and \link[=Extract]{square} brackets) as 
      these are primitive functions in \R.
	  Set \code{what} to \code{NA} to match all names.
 }
  \item{FUN}{a function to be applied. }
  \item{symbols}{logical value controlling whether \code{FUN} should be applied
      to symbols as well as calls.}
  \item{\dots}{optional arguments to \code{FUN}.}
}

\value{
    A (modified) expression.
}

\details{
\code{FUN} is found by a call to \code{\link{match.fun}} and can be either
a function or a symbol (e.g., a backquoted name) or a character string
specifying a function to be searched for from the environment of the call to
\code{exprApply}. 
}

\note{
If \code{expr} has a \link[=srcfile]{source reference} information
(\code{"srcref"} attribute), modifications done by \code{exprApply} will not be
visible when printed unless \code{srcref} is removed. However, \code{exprApply}
does remove source reference from any \code{function} expression inside
\code{expr}.
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{

Expression-related functions: \code{\link{substitute}}, 
\code{\link{expression}}, \code{\link{quote}} and \code{\link{bquote}}.

Similar function \code{walkCode} exists in package \pkg{codetools}. 

Functions useful inside \code{FUN}: \code{\link{as.name}}, \code{\link{as.call}},
\code{\link{call}}, \code{\link{match.call}} etc.

}

\examples{
### simple usage:
# print all Y(...) terms in a formula (note that symbol "Y" is omitted):
# Note: if `print` is defined as S4 "standardGeneric", we need to use 
# 'print.default' rather than 'print', or put the call to 'print' inside a 
# function, i.e. as `function(x) print(x)`:
exprApply(~ X(1) + Y(2 + Y(4)) + N(Y + Y(3)), "Y", print.default)


# replace X() with log(X, base = n)
exprApply(expression(A() + B() + C()), c("A", "B", "C"), function(expr, base) {
    expr[[2]] <- expr[[1]]
    expr[[1]] <- as.name("log")
    expr$base <- base
    expr
}, base = 10)

###
# TASK: fit lm with two poly terms, varying the degree from 1 to 3 in each.
# lm(y ~ poly(X1, degree = a) + poly(X2, degree = b), data = Cement)
# for a = {1,2,3} and b = {1,2,3}

# First we create a wrapper function for lm. Within it, use "exprApply" to add
# "degree" argument to all occurences of "poly()" having "X1" or "X2" as the
# first argument. Values for "degree" are taken from arguments "d1" and "d2"

lmpolywrap <- function(formula, d1 = NA, d2 = NA, ...) { 
    cl <- origCall <- match.call()
    cl[[1]] <- as.name("lm")
    cl$formula <- exprApply(formula, "poly", function(e, degree, x) {
        i <- which(e[[2]] == x)[1]
        if(!is.na(i) && !is.na(degree[i])) e$degree <- degree[i]
        e
    }, degree = c(d1, d2), x = c("X1", "X2"))
    cl$d1 <- cl$d2 <- NULL
    fit <- eval(cl, parent.frame())
    fit$call <- origCall # replace the stored call
    fit
}

# global model:
fm <- lmpolywrap(y ~ poly(X1) + poly(X2), data = Cement)

# Use "dredge" with argument "varying" to generate calls of all combinations of
# degrees for poly(X1) and poly(X2). Use "fixed = TRUE" to keep all global model
# terms in all models.
# Since "dredge" expects that global model has all the coefficients the 
# submodels can have, which is not the case here, we first generate model calls,
# evaluate them and feed to "model.sel"

modCalls <- dredge(fm, 
    varying = list(d1 = 1:3, d2 = 1:3), 
    fixed = TRUE,
    evaluate = FALSE
)

model.sel(models <- lapply(modCalls, eval))

# Note: to fit *all* submodels replace "fixed = TRUE" with: 
# "subset = (d1==1 || {poly(X1)}) && (d2==1 || {poly(X2)})"
# This is to avoid fitting 3 identical models when the matching "poly()" term is
# absent.
}

\keyword{manip}
