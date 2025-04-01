\name{model.selection.object}
\alias{model.selection.object}
\title{Description of Model Selection Objects}
\description{
  An object of class \code{"model.selection"} holds a table of model 
  coefficients and ranking statistics. It is produced by \lcode{dredge} 
  or \lcode{model.sel}.
}

\value{
  The object is a \code{data.frame} with additional attributes. Each row 
  represents one model. The models are ordered by the information criterion 
  value specified by \code{rank} (lowest on top).

Data frame columns:
  
\item{<model terms>}{For numeric covariates these columns hold coefficent value,
    for factors their presence in the model. If the term is not present in a 
    model, value is \code{NA}. }
\item{<varying arguments>}{Optional. If any arguments differ between the 
    modelling function calls (except for formulas and some other arguments), 
    these will be held in additional columns (of class \code{"factor"}).}
\item{df}{Number of model parameters}
\item{logLik}{Log-likelihood (or quasi-likelihood for \acronym{GEE})}
\item{<rank>}{Information criterion value}
\item{delta}{the \IC difference, i.e. the the relative difference to 
    the best model, 
    \ifelse{latex}{\eqn{\Delta_{IC} = IC_i - IC_{min}}}{\enc{Δ}{Delta}_IC = IC_i - IC_min}, }
\item{weight}{\sQuote{Akaike weights}, i.e. \link[=Weights]{normalized model likelihoods}.}

Attributes:
    \item{model.calls}{A list containing model calls (arranged in 
        the same order as in the table). A model call can be retrieved with
        \code{getCall(*, i)} where \var{i} is a vector of model index or name 
        (if given as character string). }
    \item{global}{The \code{global.model} object }
    \item{global.call}{Call to the \code{global.model} }
    \item{terms}{A character string holding all term names. Attribute 
        \code{"interceptLabel"} gives the name of the intercept term. }
    \item{rank}{The \code{rank} function used }
    \item{beta}{A character string, representing the coefficient standardizing 
        method used. Either \code{"none"}, \code{"sd"} or \code{"partial.sd"}
    }
    \item{coefTables}{List of matrices of class \code{"coefTable"} containing 
        each model's coefficients with std. errors and associated \var{df}s }
    \item{nobs}{Number of observations }
    \item{warnings}{optional (\code{pdredge} only). A list of errors and 
        warnings issued by the modelling function during the fitting, with a model
        number appended to each. }


It is not recommended to directly access the attributes. Instead, use extractor
    functions if possible. These include \code{getCall} for retrieving model
    calls, \code{coefTable} and \code{coef} for coefficients, 
    and \code{nobs}. \code{logLik} extracts list of model log-likelihoods (as 
    \code{"logLik"} objects), and \code{Weights} extracts \sQuote{Akaike 
    weights}.

The object has class \code{c("model.selection", "data.frame")}.    
}

\seealso{
  \lcode{dredge}, \lcode{model.sel}.
}
\keyword{models}