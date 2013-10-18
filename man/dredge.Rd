\name{dredge}
\alias{dredge}
\alias{print.model.selection}

\encoding{utf-8}
\title{Automated model selection}
\description{
Generate a set of models with combinations (subsets) of the terms in the global
model, with optional rules for model inclusion.
}

\usage{
dredge(global.model, beta = FALSE, evaluate = TRUE, rank = "AICc", 
    fixed = NULL, m.max = NA, m.min = 0, subset, marg.ex = NULL, 
    trace = FALSE, varying, extra, ct.args = NULL, ...)

\method{print}{model.selection}(x, abbrev.names = TRUE, warnings = getOption("warn") != -1L, ...)

}

\arguments{
	\item{global.model}{a fitted \sQuote{global} model object. See
	\sQuote{Details} for a list of supported types. }

	\item{beta}{logical, should standardized coefficients be returned? }

	\item{evaluate}{whether to evaluate and rank the models. If \code{FALSE}, a
	list of model \code{call}s is returned. }

	\item{rank}{optional custom rank function (information criterion) to be used
	instead \code{AICc}, e.g. \code{AIC}, \code{QAIC} or \code{BIC}.
	See \sQuote{Details}. }

	\item{fixed}{optional, either a single sided formula or a character vector
	giving names of terms to be included in all models. See \sQuote{Subsetting}. }

	\item{m.max, m.min}{optionally, the maximum and minimum number of terms in a
	single model (excluding the intercept), \code{m.max} defaults to the number
	of terms in \code{global.model}. See \sQuote{Subsetting}. }

	\item{subset}{ logical expression describing models to keep in the resulting
		set. See \sQuote{Subsetting}. }

	\item{marg.ex}{a character vector specifying names of variables for which
	NOT to check for marginality restrictions when generating model formulas. If
	this argument is set to \code{TRUE}, all combinations of terms are used
	(i.e. no checking). If \code{NA} or missing, the exceptions will be found
	based on the terms of \code{global.model}. See \sQuote{Details}.  }

	\item{trace}{if \code{TRUE}, all calls to the fitting function (i.e. updated
	\code{global.model} calls) are printed before actual fitting takes place. }

	\item{varying}{ optionally, a named list describing the additional arguments
	to vary between the generated models. Item names correspond to the
	arguments, and each item provides a list of choices (i.e. \code{list(arg1 =
	list(choice1, choice2, ...), ...)}). Complex elements in the choice list
	(such as \code{family} objects) should be either named (uniquely) or quoted
	(unevaluated, e.g. using \code{\link{alist}}, see \code{\link{quote}}),
	otherwise it may produce rather unpleasant effects. See example in
	\code{\link{Beetle}}. }

	\item{extra}{ optional additional statistics to include in the result,
	provided as functions, function names or a list of such (best if named
	or quoted). Similarly as in \code{rank} argument, each function must accept
	fitted model object as an argument and return a (value coercible to a)
	numeric vector.
	These can be e.g. additional information criterions or goodness-of-fit
	statistics. The character strings "R^2" and "adjR^2" are treated in a
	special way, and will add a likelihood-ratio based \eqn{R^{^2}}{R^2} and
	modified-\eqn{R^{^2}}{R^2} respectively to the result (this is more
	efficient than using \code{\link{r.squaredLR}} directly).
	}

  \item{x}{
    a \code{model.selection} object, returned by \code{dredge}. }

  \item{abbrev.names}{ should printed variable names be abbreviated?
	(useful with many variables). }

  \item{warnings}{ if \code{TRUE}, errors and warnings issued during the model
	fitting are printed below the table (currently, only with \code{pdredge}).
	To permanently remove the warnings, set the object's attribute "warnings" to
	\code{NULL}. }
	
  \item{ct.args}{ optional list of arguments to be passed to
	\code{\link{coefTable}} (e.g. \code{dispersion} parameter for \code{glm}
		affecting standard errors used in subsequent 
		\code{\link[=model.avg]{model averaging}}).}

  \item{\dots}{ optional arguments for the \code{rank} function. Any can be
    an expression (of mode \code{call}), in which case any \code{x} within it
    will be substituted with a current model. }
}


\details{

Models are fitted one by one through repeated evaluation of modified calls to
the \code{global.model} (in a similar fashion as with \code{update}). This
approach, while robust in that it can be applied to a variety of different model
object types is not very efficient and may be time-intensive.

Note that the number of combinations grows exponentially with number of
predictor variables (\ifelse{latex}{\eqn{2^{N}}}{2^N}, less when interactions
are present, see below). As there can be potentially a large number of models to
evaluate, to avoid memory overflow the fitted model objects are not stored in
the result. To get (a subset of) the models, use \code{\link{get.models}} on the
object returned by \code{dredge}.

For a list of model types that can be used as a \code{global.model} see
\link[=MuMIn-models]{list of supported models}. 
Modelling functions not storing \code{call} in their result should be evaluated 
\emph{via} the wrapper created by \code{\link{updateable}}.

\subsection{Information criterion}{
\code{rank} is found by a call to \code{match.fun} and may be specified as a
function or a symbol (e.g. a back-quoted name) or a character string specifying
a function to be searched for from the environment of the call to \code{dredge}.
The function \code{rank} must accept model object as its first argument and 
always return a scalar. Typical choice for \code{rank} would be "AIC", "BIC", or
"QAIC".
}

\subsection{Interactions}{
\code{dredge} by default respects marginality
constraints, so \dQuote{all possible combinations} do not include models
containing interactions without their respective main effects and all lower order 
terms. This behaviour can be altered by \code{marg.ex} argument, which can be used 
to allow for simple nested designs. For example, with global model of form 
\code{a / (x + z)}, one would use \code{marg.ex = "a"} and \code{fixed = "a"}. 
If \code{global.model} uses such a formula and \code{marg.ex} is missing or 
\code{NA}, it will be adjusted automatically.
}

\subsection{Subsetting}{
There are three ways to constrain the resulting set of models: setting limits to
the number of terms in a model with \code{m.max} and \code{m.min}, binding
term(s) to all models with \code{fixed}, and more complex rules can be applied
using argument \code{subset}. To be included in the selection table, the model
formulation must satisfy all these conditions.

%%Terms in \code{fixed} argument are applied before the combinations are
%%generated, therefore more efficient than \code{subset}.

\code{subset} can take either a form of an \emph{expression} or a \emph{matrix}.
The latter should be a lower triangular matrix with logical values, where
columns and rows correspond to \code{global.model} terms. Value
\code{subset["a", "b"] == FALSE} will exclude any model containing both terms 
\var{a} and \var{b}. Values other than \code{FALSE} (or \code{0}) are taken as
\code{TRUE}.

In the form of \code{expression}, the argument \code{subset} acts in a similar
fashion to that in the function \code{subset} for \code{data.frames}: model
terms can be referred to by name as variables in the expression, with the
difference being that they are always logical (i.e. \code{TRUE} if a term exists
in the model).

There is also \code{.(x)} and \code{.(+x)} notation that indicate respectively
any or all model terms including a \emph{variable} \code{x}. This concerns only
interactions containing a particular main effect \code{x} (e.g. \code{x:z},
\code{v:x:z}), and \emph{not} a variable in a complex expression, such as
\code{x} in \code{log(x)} or \code{I(x^2)}. This is only useful with marginality
 exceptions, see \code{marg.ex} argument.

The expression can contain any of the \code{global.model} terms
(\code{getAllTerms(global.model)} lists them), as well as names of the
\code{varying} argument items. Names of \code{global.model} terms take
precedence when identical to names of \code{varying}, so to avoid ambiguity
\code{varying} variables in \code{subset} expression should be enclosed in
\code{V()} (e.g. \code{subset = V(family) == "Gamma"} assuming that
\code{varying} is something like \code{list(family = c(..., "Gamma"))}).

If element names in \code{varying} are missing, the elements themselves are 
used. Call and symbol elements are represented as character values (via 
\code{deparse}), and everything except numeric, logical, character and 
\code{NULL} values is replaced by item numbers (e.g. \code{varying = 
list(family = list(..., Gamma)} should be referred to as \code{subset = 
V(family) == 2}. This can quickly become confusing, therefore it is recommended 
to use named lists in most cases. \code{demo(dredge.varying)} provides examples.

The \code{subset} expression can also contain variable \code{`*nvar*`} (needs to
be backtick-quoted), which is equal to number of terms in the model (\bold{not}
the number of estimated parameters \var{K}).

To make inclusion of a variable conditional on presence of some other variable,
a function \code{dc} (\dQuote{\bold{d}ependency \bold{c}hain}) can be used in
the \code{subset} expression. \code{dc} takes any number of variables as
arguments, and allows a variable to be included only if all preceding variables
are also present (e.g. \code{subset = dc(a, b, c)} allows for models of form
\code{a}, \code{a+b }and \code{a+b+c} but not \code{b}, \code{c}, \code{b+c} or
\code{a+c}).

\code{subset} expression can have a form of an unevaluated \code{call},
\code{expression} object, or a one sided \code{formula}. See \sQuote{Examples}.
Compound model terms (such as \sQuote{as-is} expressions within \code{I()} or
smooths in \code{gam}) should be treated as non-syntactic names and enclosed in
back-ticks (e.g. \code{subset = `s(x, k = 2)` || `I(log(x))`}, see
\link[base]{Quotes}). Mind the spacing, names must match exactly the term names
in model's formula. To simply keep certain terms in all models, use of argument
\code{fixed} is more efficient. Note that the \code{fixed} formula is 
interpreted in the same manner as model formula and so the terms need not to be 
quoted.

\subsection{\code{subset} expression syntax summary}{
	\tabular{ll}{
	\code{a & b} \tab indicates that model terms \var{a} and \var{b} must be
	present (see \sQuote{\link[=Logic]{Logical Operators}}) \cr 
	\code{V(x)} \tab indicates a \code{varying} variable \var{x} \cr
	\code{.(x)} \tab indicates that at least one term containing variable \var{x} 
		must be present \cr
	\code{.(+x)} \tab indicates that all the terms containing variable \var{x} 
		must be present \cr
	\code{dc(a, b, c,...)} \tab \sQuote{dependency chain}: \var{a} is allowed only
		if \var{b} is present, and \var{b} only if \var{c} is present, etc. \cr
	\code{`*nvar*`} \tab number of variables \cr
	}
}
}

\subsection{Missing values}{
Use of \code{na.action = na.omit} (\R's default) in \code{global.model} should
be avoided, as it results with sub-models fitted to different data sets, if
there are missing values. Warning is given if it is detected.
}

\subsection{Methods}{
There are \code{\link[=subset.model.selection]{subset}} and \code{plot} methods,
the latter produces a graphical representation of model weights and variable
relative importance.
Coefficients can be extracted with \code{coef} or \code{\link{coefTable}}.
}

}


\value{
\code{dredge} returns an object of class \code{model.selection}, being a
\code{data.frame} with models' coefficients (or presence/\code{NA} for factors),
	\emph{df} - number of parameters, log-likelihood, the information criterion
	value, delta-IC and \emph{Akaike weight}. Models are ordered by the value of
	the information criterion specified by \code{rank} (lowest on top).

The attribute \code{"calls"} is a list containing the model calls used (arranged
in the same order as the models). A model call can be retrieved with
\code{getCall(*, i)} where \var{i} is a vector of model index or name (if given
as character string).


Other attributes:
\code{"global"} - the \code{global.model} object,
\code{"rank"} - the \code{rank} function used, \code{"call"} - the matched
call, and
\code{"warnings"} - list of errors and warnings given by the modelling function
	during the fitting, with model number appended to each. 
}

\author{Kamil Barto\enc{ń}{n}}

\note{
Users should keep in mind the hazards that a \dQuote{thoughtless approach}
of evaluating all possible models poses. Although this procedure is in certain
cases useful and justified, it may result in selecting a spurious \dQuote{best}
model, due to the model selection bias.

\emph{\dQuote{Let the computer find out} is a poor strategy and usually reflects
the fact that the researcher did not bother to think clearly about the problem
of interest and its scientific setting} (Burnham and Anderson, 2002).

}

\seealso{

\code{\link{pdredge}} is a parallelized version of this function (uses a
cluster).

\code{\link{get.models}}, \code{\link{model.avg}}. \code{\link{model.sel}} for
manual model selection tables.

Possible alternatives: \code{\link[glmulti]{glmulti}} in package \pkg{glmulti}
and \code{\link[bestglm]{bestglm}} (\pkg{bestglm}).
\code{\link[leaps]{regsubsets}} in package \pkg{leaps} also performs all-subsets
regression.
}


\examples{
# Example from Burnham and Anderson (2002), page 100:
data(Cement)
fm1 <- lm(y ~ ., data = Cement)
dd <- dredge(fm1)
subset(dd, delta < 4)

# Visualize the model selection table:
if(require(graphics))
    plot(dd)


# Model average models with delta AICc < 4
model.avg(dd, subset = delta < 4)

#or as a 95\% confidence set:
model.avg(dd, subset = cumsum(weight) <= .95) # get averaged coefficients

#'Best' model
summary(get.models(dd, 1))[[1]]

\dontrun{
# Examples of using 'subset':
# keep only models containing X3
dredge(fm1, subset = ~ X3) # subset as a formula
dredge(fm1, subset = expression(X3)) # subset as expression object
# the same, but more effective:
dredge(fm1, fixed = "X3")
# exclude models containing both X1 and X2 at the same time
dredge(fm1, subset = !(X1 && X2))
# Fit only models containing either X3 or X4 (but not both);
# include X3 only if X2 is present, and X2 only if X1 is present.
dredge(fm1, subset = dc(X1, X2, X3) && xor(X3, X4))
# the same as above, but without using "dc"
dredge(fm1, subset = (X1 | !X2) && (X2 | !X3) && xor(X3, X4))

# Include only models with up to 2 terms (and intercept)
dredge(fm1, m.max = 2)
}


# Add R^2 and F-statistics, use the 'extra' argument
dredge(fm1, m.max = 1, extra = c("R^2", F = function(x)
    summary(x)$fstatistic[[1]]))

# with summary statistics:
dredge(fm1, m.max = 1, extra = list(
    "R^2", "*" = function(x) {
        s <- summary(x)
        c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
            F = s$fstatistic[[1]])
    })
)


# Add other information criterions (but rank with AICc):
\dontshow{
# there is no BIC in R < 2.13.0, so need to add it:
if(!exists("BIC", mode = "function"))
    BIC <- function(object, ...)
        AIC(object, k = log(length(resid(object))))
}
dredge(fm1, m.max = 1, extra = alist(AIC, BIC, ICOMP, Cp))


}

\keyword{models}
