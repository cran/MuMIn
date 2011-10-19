\name{model.sel}
\alias{mod.sel}
\alias{model.sel}
\alias{mod.sel.default}
\alias{mod.sel.model.selection}
\alias{model.sel}
\encoding{utf-8}

\title{model selection table}
\description{
Builds a model selection table
}

\usage{

mod.sel(object, ...)
model.sel(object, ...)

\method{mod.sel}{model.selection}(object, rank = NULL, rank.args = NULL, ...)
\method{mod.sel}{default}(object, ..., rank = NULL, rank.args = NULL)

}

\arguments{
	\item{object}{A fitted model object, a list of such objects, or a
	\code{"model.selection"} object.}

    \item{\dots}{ More fitted model objects }

    \item{rank}{Optional, custom rank function (information criterion) to use
        instead of \code{AICc}, e.g. \code{QAIC} or \code{BIC}, may be omitted
        if \code{object} is a model list returned by \code{get.models}. }

    \item{rank.args}{Optional \code{list} of arguments for the \code{rank}
        function. If one is an expression, an \code{x} within it is substituted
        with a current model.}

}

\value{

An object of class \code{"model.selection"} with columns containing useful information
about each model: the coefficients, value of the information criterion used, Delta(IC) and weight.

}

% \details{
% }

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{dredge}}
}

\examples{

data(Cement)
Cement$X1 <- cut(Cement$X1, 3)
Cement$X2 <- cut(Cement$X2, 2)

fm1 <- glm(formula = y ~ X + X1 + X2 * X3, data = Cement)
fm2 <- update(fm1, . ~ . - X - X1)
fm3 <- update(fm1, . ~ . - X2 - X3)

# ranked with AICc by default
mod.sel(fm1, fm2, fm3)


# ranked with BIC
mod.sel(fm1, fm2, fm3, rank=AIC, rank.args=alist(k=log(nobs(x))))

}

\keyword{models}
