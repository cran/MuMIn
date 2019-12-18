\name{stdize}
\alias{stdize}
\alias{stdize.default}
\alias{stdize.logical}
\alias{stdize.formula}
\alias{stdize.data.frame}
\alias{stdizeFit}

\encoding{utf-8}

\title{Standardize data}
\description{
\code{stdize} standardizes variables by centring and scaling.

\code{stdizeFit} modifies a model call or existing model to use standardized 
variables.
}

\usage{
\method{stdize}{default}(x, center = TRUE, scale = TRUE, ...)

\method{stdize}{logical}(x, binary = c("center", "scale", "binary", "half", "omit"),
    center = TRUE, scale = FALSE, ...)
## also for two-level factors

\method{stdize}{data.frame}(x, binary = c("center", "scale", "binary", "half", "omit"),
    center = TRUE, scale = TRUE, omit.cols = NULL, source = NULL,
    prefix = TRUE, append = FALSE, ...)

\method{stdize}{formula}(x, data = NULL, response = FALSE,
    binary = c("center", "scale", "binary", "half", "omit"),
    center = TRUE, scale = TRUE, omit.cols = NULL, prefix = TRUE,
    append = FALSE, ...)

stdizeFit(object, data, which = c("formula", "subset", "offset", "weights"),
    evaluate = TRUE, quote = NA)
}

\arguments{
    \item{x}{a numeric or logical vector, factor, numeric matrix, 
         \code{data.frame} or a formula. }
    \item{center, scale}{ either a logical value or a logical or numeric vector
         of length equal to the number of columns of \code{x} (see
         \sQuote{Details}). \code{scale} can be also a function to use for
         scaling. 
         }
    \item{binary}{specifies how binary variables (logical or two-level factors)
        are scaled. Default is to \code{"center"} by subtracting the mean
        assuming levels are equal to 0 and 1; use \code{"scale"} to
        both centre and scale by \acronym{SD}, \code{"binary"} to centre to 0 /
        1, \code{"half"} to centre to -0.5 / 0.5, and \code{"omit"} to leave
        binary variables unmodified.
        This argument has precedence over \code{center} and \code{scale}, unless
        it is set to \code{NA} (in which case binary variables are treated like
        numeric variables).
        }
     \item{source}{a reference \code{data.frame}, being a result of previous  
         \code{stdize}, from which \code{scale} and \code{center} values are 
         taken. Column names are matched. This can be used for scaling new data 
         using statistics of another data.  
          }
    \item{omit.cols}{ column names or numeric indices of columns that should 
         be left unaltered. }
    \item{prefix}{ either a logical value specifying whether the names of
        transformed columns should be prefixed, or a two-element character vector
        giving the prefixes. The prefixes default to \dQuote{z.} for scaled and
        \dQuote{c.} for centred variables. }
        
    \item{append}{ logical, if \code{TRUE}, modified columns are appended to the
        original data frame. }
    \item{response}{logical, stating whether the response should be standardized. By
        default, only variables on the right-hand side of the formula are standardized. }
    \item{data}{ an object coercible to \code{data.frame}, containing the
        variables in \code{formula}. Passed to, and used by 
        \code{\link{model.frame}}.
        
        For \code{stdizeFit}, a \code{stdize}d \code{data.frame} to use.
        }
    \item{\dots}{ for the \code{formula} method, additional arguments passed to
        \code{\link{model.frame}}. For other methods, it is silently ignored. }
    \item{object}{a fitted model object or an expression being a \code{call} to
        the modelling function. }
    \item{which}{ a character string naming arguments which should be modified.
        This should be all arguments which are evaluated in the \code{data}
        environment. Can be also \code{TRUE} to modify the expression as a
        whole. The \code{data} argument is additionally replaced with that
        passed to \code{stdizeFit}.
        }
    \item{evaluate}{ if \code{TRUE}, the modified call is evaluated and the
        fitted model object is returned. }
    \item{quote}{ if \code{TRUE}, avoids evaluating \code{object}. Equivalent to
        \code{stdizeFit(quote(expr), ...)}. Defaults to \code{NA} in which case
        \code{object} being a call to non-primitive function is quoted. 
        }
}

\value{

\code{stdize} returns a vector or object of the same dimensions as \code{x},
where the values are centred and/or scaled. Transformation is carried out
column-wise in \code{data.frame}s and matrices. 

The returned value is compatible with that of \code{\link{scale}} in that the
numeric centring and scalings used are stored in attributes
\code{"scaled:center"} and \code{"scaled:scale"} (these can be \code{NA} if no
centring or scaling has been done).

\code{stdizeFit} returns a modified, unevaluated call where the variable names
are replaced to point the transformed variables, or if \code{evaluate} is
\code{TRUE}, a fitted model object.
}

\details{
\code{stdize} resembles \code{\link{scale}}, but uses special rules
for factors, similarly to \code{\link[arm]{standardize}} in package \pkg{arm}.

\code{stdize} differs from \code{\link[arm]{standardize}} in that it is used on
data rather than on the fitted model object. The scaled data should afterwards
be passed to the modelling function, instead of the original data.

Unlike \code{standardize}, it applies special \sQuote{binary} scaling only to
two-level \code{factor}s and logical variables, rather than to any variable with
two unique values.

Variables of only one unique value are unchanged.

By default, \code{stdize} scales by dividing by standard deviation rather than twice
the \acronym{SD} as \code{standardize} does. Scaling by \acronym{SD} is used 
also on uncentred values, which is different from \code{\link{scale}} where
root-mean-square is used.

If \code{center} or \code{scale} are logical scalars or vectors of length equal
to the number of columns of \code{x}, the centring is done by subtracting the
mean (if \code{center} corresponding to the column is \code{TRUE}), and scaling
is done by dividing the (centred) value by standard deviation (if corresponding
\code{scale} is \code{TRUE}). 
If \code{center} or \code{scale} are numeric vectors with length equal 
to the number of columns of \code{x} (or numeric scalars for vector methods),
then these are used instead. Any \code{NA}s in the numeric vector result in no
centring or scaling on the corresponding column.

Note that \code{scale = 0} is equivalent to no scaling (i.e. \code{scale = 1}).

Binary variables, logical or factors with two levels, are converted to
numeric variables and transformed according to the argument \code{binary},
unless \code{center} or \code{scale} are explicitly given.

}

\author{Kamil Barto\enc{ń}{n}}

\references{
Gelman, A. (2008) Scaling regression inputs by dividing by two standard 
deviations. \emph{Statistics in medicine} 27, 2865-2873.
}


\seealso{
Compare with \code{\link{scale}} and \code{\link[arm]{standardize}} or
\code{\link[arm]{rescale}} (the latter two in package \pkg{arm}).

For typical standardizing, model coefficients transformation may be
easier, see \code{\link{std.coef}}.

\code{\link{apply}} and \code{\link{sweep}} for arbitrary transformations of 
columns in a \code{data.frame}.



}

\examples{
# compare "stdize" and "scale"
nmat <- matrix(runif(15, 0, 10), ncol = 3)

stdize(nmat)
scale(nmat)

rootmeansq <- function(v) {
    v <- v[!is.na(v)]
    sqrt(sum(v^2) / max(1, length(v) - 1L))
}

scale(nmat, center = FALSE)
stdize(nmat, center = FALSE, scale = rootmeansq)

if(require(lme4)) {
# define scale function as twice the SD to reproduce "arm::standardize"
twosd <- function(v) 2 * sd(v, na.rm = TRUE)

# standardize data (scaled variables are prefixed with "z.")
z.CO2 <- stdize(uptake ~ conc + Plant, data = CO2, omit = "Plant", scale = twosd)
summary(z.CO2)


fmz <- stdizeFit(lmer(uptake ~ conc + I(conc^2) + (1 | Plant)), data = z.CO2)
# produces:
# lmer(uptake ~ z.conc + I(z.conc^2) + (1 | Plant), data = z.CO2)


## standardize using scale and center from "z.CO2", keeping the original data:
z.CO2a <- stdize(CO2, source = z.CO2, append = TRUE)
# Here, the "subset" expression uses untransformed variable, so we modify only
# "formula" argument, keeping "subset" as-is. For that reason we needed the
# untransformed variables in "data".
stdizeFit(lmer(uptake ~ conc + I(conc^2) + (1 | Plant),
    subset = conc > 100,
    ), data = z.CO2a, which = "formula", evaluate = FALSE)


# create new data as a sequence along "conc"
newdata <-  data.frame(conc = seq(min(CO2$conc), max(CO2$conc), length = 10))

# scale new data using scale and center of the original scaled data: 
z.newdata <- stdize(newdata, source = z.CO2)

\dontshow{ if(require(graphics)) \{ }
# plot predictions against "conc" on real scale:
plot(newdata$conc, predict(fmz, z.newdata, re.form = NA))
\dontshow{ \} }

# compare with "arm::standardize"
\dontrun{
library(arm)
fms <- standardize(lmer(uptake ~ conc + I(conc^2) + (1 | Plant), data = CO2))
plot(newdata$conc, predict(fms, z.newdata, re.form = NA))
}
}

}

\keyword{manip}
