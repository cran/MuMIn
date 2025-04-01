\name{nested}
\alias{nested}
\encoding{utf-8}

\title{Identify nested models}
\description{
Find models that are \sQuote{nested} within each model in the model 
selection table.
}

\usage{
nested(x, indices = c("none", "numeric", "rownames"), rank = NULL)
}

\arguments{
	\item{x}{a \code{"model.selection"} object (result of \code{dredge} or 
	\code{model.sel}). }
	
	\item{indices}{if omitted or \code{"none"} then the function checks if, for 
		each model, there are any higher ranked models nested within it.
		If \code{"numeric"} or \code{"rownames"}, indices or names of all nested
		models are returned. See \dQuote{Value}.
		}
		
	\item{rank}{the name of the column with the ranking values (defaults to 
	    the one before \dQuote{delta}). Only used if \code{indices} is 
		\code{"none"}. 
		}
}

\value{
	A vector of length equal to the number of models (table rows).
	
	If \code{indices = "none"} (the default), it is a vector of logical 
	values where \emph{i}-th element is \code{TRUE} if any model(s) higher up in
	the table are nested within it (i.e. if simpler models have lower IC pointed
	by \code{rank}).
	
	For \code{indices} other than \code{"none"}, the function returns a list of 
		vectors of numeric indices or names of models nested within each 
		\emph{i}-th model.
}

\details{
	In model comparison, a model is said to be \dQuote{nested} within another 
	model if it contains a subset of parameters of the latter model, but does not
	include other parameters (e.g. model \sQuote{A+B} is nested within 
	\sQuote{A+B+C} but not \sQuote{A+C+D}). 

	This function can be useful in a model selection approach suggested by 
	Richards (2008), in which more complex variants of any model with a lower 
	IC value are excluded from the candidate set. 
	
}

\note{
This function determines nesting based only on fixed model terms, within groups of 
models sharing the same \sQuote{varying} parameters (see \lcode{dredge} and 
example in \lcode{Beetle}).
}


\author{Kamil Barto\enc{ń}{n}}

\references{

Richards, S. A., Whittingham, M. J., Stephens, P. A. 2011 Model selection and 
model averaging in behavioural ecology: the utility of the IT-AIC framework. 
\emph{Behavioral Ecology and Sociobiology} \bold{65}, 77--89.

Richards, S. A. 2008 Dealing with overdispersed count data in applied
ecology. \emph{Journal of Applied Ecology} \bold{45}, 218--227.

}

\seealso{
\lcode{dredge}, \lcode{model.sel}

}

\examples{

fm <- lm(y ~ X1 + X2 + X3 + X4, data = Cement, na.action = na.fail)
ms <- dredge(fm)

# filter out overly complex models according to the 
# "nesting selection rule":
subset(ms, !nested(.)) # dot represents the ms table object

# print model "4" and all models nested within it
nst <- nested(ms, indices = "row")
ms[c("4", nst[["4"]])]

ms$nested <- sapply(nst, paste, collapse = ",")

ms

}

\keyword{models}
