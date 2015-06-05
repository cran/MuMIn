\name{predict_avg}
\alias{predict_avg}
\alias{avg.predictions}
\alias{std_predict}
\alias{std_predict.lm}
\encoding{utf-8}
\title{
 Average predictions (experimental)
}
\description{
 Average predictions (experimental functions under development)
}
\usage{
predict_avg(x, newdata, type = c("response", "invlink", "link", "terms"), 
    se.fit = FALSE, full = TRUE, linkinv = NULL, use.lincomb = FALSE, ...)
	
avg.predictions(yall, w, type = c("atomic", "matrix", "se.fit", "terms"), 
    revised.var = TRUE, full = FALSE, ...) 

%% std_predict(object, ...)
\method{std_predict}{lm}(object, newdata, type, se.fit, ...)

}
\arguments{
  \item{x,object}{model averaging object.}
  \item{yall}{list of predictions.}
  \item{newdata}{new data.}
  \item{type}{type of prediction.}
  \item{se.fit}{return SE? }
  \item{w}{weights}
  \item{revised.var}{revised.var?}
  \item{full}{logical, full or subset prediction }
  \item{linkinv}{inverse link (as function, family, link-glm, or character). 
  \emph{Caution! \code{linkinv = "log"} will transform the prediction on a link scale using
exponential function (via \code{make.link("log")}, whereas \code{linkinv = log} will
be interpreted as a function and log-transform them.} }
  \item{use.lincomb}{logical, whether to predict using linear combination of 
  model-averaged coefficients? Note that this currently does not compute predictions properly
  for models using fitting weights and offset. Component model must have 
  \code{model.matrix} method.
  
  If \code{use.lincomb = FALSE} predictions from component models are averaged.
}
  \item{\dots}{other arguments...}
}
\details{TODO}
\value{TODO}
\author{Kamil Barto\enc{ń}{n}}


\seealso{predict.averaging, predict, predict.glm}
\examples{
# TODO
}
\keyword{manip}

