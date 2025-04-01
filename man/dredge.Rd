\name{dredge}
\alias{dredge}
\alias{dc}
\alias{V}
\alias{print.model.selection}
\encoding{utf-8}

\title{Automated model selection}
\description{
Generate a model selection table of models with combinations (subsets) of fixed 
effect terms in the global model, with optional model inclusion rules.
}

\usage{
dredge(global.model, beta = c("none", "sd", "partial.sd"), evaluate = TRUE,
  rank = "AICc", fixed = NULL, m.lim = NULL, m.min, m.max, subset,
  trace = FALSE, varying, extra, ct.args = NULL, deps = attr(allTerms0, "deps"),
  cluster = NULL,
  ...)

\method{print}{model.selection}(x, abbrev.names = TRUE, warnings = getOption("warn") != -1L, ...)

}

\arguments{
	\item{global.model}{a fitted \sQuote{global} model object. See
	\sQuote{Details} for a list of supported types. }

	\item{beta}{indicates whether and how the coefficients are standardized, and 
    must be one of \code{"none"}, \code{"sd"} or \code{"partial.sd"}. You can 
    specify just the initial letter. \code{"none"} corresponds to unstandardized
    coefficients, \code{"sd"} and \code{"partial.sd"} to coefficients 
    standardized by \acronym{SD} and Partial \acronym{SD}, respectively. For 
    backwards compatibility, logical value is also accepted, \code{TRUE} is 
    equivalent to \code{"sd"} and \code{FALSE} to \code{"none"}. 
    See \lcode{std.coef}.	}

	\item{evaluate}{whether to evaluate and rank the models. If \code{FALSE}, a
	list of unevaluated \code{call}s is returned. }

	\item{rank}{optionally, the rank function returning a sort of an information
    criterion, to be used instead \code{AICc}, e.g. \code{AIC}, \code{QAIC} or
    \code{BIC}.
	See \sQuote{Details}. }

	\item{fixed}{optional, either a single-sided formula or a character vector
	giving names of terms to be included in all models. Not to be confused with
    fixed effects. See \sQuote{Subsetting}. }

	\item{m.lim, m.max, m.min}{optionally, the limits \code{c(lower, upper)}
    for the number of terms in a single model (excluding the intercept). An
    \code{NA} means no limit. See \sQuote{Subsetting}.
    Specifying limits as \code{m.min} and \code{m.max} is allowed for backward
    compatibility. }

	\item{subset}{logical expression or a \code{matrix} describing models to 
        be kept in the resulting set. \code{NULL} or \code{TRUE} 
        disables subsetting. For details, see \sQuote{Subsetting}. }

	\item{trace}{if \code{TRUE} or \code{1}, all calls to the fitting function
	are printed before actual fitting takes place. If \code{trace > 1}, a progress bar
	is displayed. }

	\item{varying}{optionally, a named list describing the additional arguments
	to vary between the generated models. Item names correspond to the
	arguments, and each item provides a list of choices (i.e. \code{list(arg1 =
	list(choice1, choice2, ...), ...)}). Complex elements in the choice list
	(such as \code{family} objects) should be either named (uniquely) or quoted
	(unevaluated, e.g. using \lcode{alist}, see \lcode{quote}),
	otherwise the result may be visually unpleasant. See example in
	\lcode{Beetle}. }

	\item{extra}{optional additional statistics to be included in the result,
	provided as functions, function names or a list of such (preferably named
	or quoted). As with the \code{rank} argument, each function must accept as 
    an argument a fitted model object and return (a value coercible to) a
	numeric vector.
	This could be, for instance, additional information criteria or goodness-of-fit
	statistics. The character strings \code{"R^2"} and \code{"adjR^2"} are
	treated in a special way and add a likelihood-ratio based \Rsq and
	modified-\Rsq to the result, respectively (this is more efficient than using
	\lcode{r.squaredLR} directly). }

  \item{x}{a \code{model.selection} object, returned by \code{dredge}. }

  \item{abbrev.names}{Should term names in the table header be abbreviated when
    printed? This is the default. If full names are required, use \code{print()} 
    explicitly with this argument set to \code{FALSE}.
       }

  \item{warnings}{if \code{TRUE}, errors and warnings issued during the model
	fitting are printed below the table (only with \code{pdredge}).
	To permanently remove the warnings, set the object's attribute
	\code{"warnings"} to \code{NULL}. }

  \item{ct.args}{optional list of arguments to be passed to
	\lcode{coefTable} (e.g. \code{dispersion} parameter for \code{glm}
		affecting standard errors used in subsequent
		\lxcode{model averaging}{=model.avg}).}
        
  \item{deps}{a \dQuote{dependency matrix} as returned by \code{getAllTerms}, 
        attribute \code{"deps"}. Can be used to fine-tune marginality 
        exceptions. }

  \item{cluster}{if a valid \code{"cluster"} object is given, it is used for
       parallel execution. If \code{NULL} or omitted, execution is 
       single-threaded.

       With parallel calculation, an extra argument \code{check} is accepted.
       
       See \lcode{pdredge} for details and examples.
       }

  \item{\dots}{optional arguments for the \code{rank} function. Any can be
    an unevaluated expression, in which case any \code{x} within it will be 
	substituted with the current model. }
}


\details{
Models are fitted through repeated evaluation of the modified call extracted from
the \code{global.model} (in a similar fashion to \code{update}). This
approach, while having the advantage that it can be applied to most model types through the
usual formula interface, can have a considerable computational overhead.

Note that the number of combinations grows exponentially with the number of
predictors (\ifelse{latex}{\eqn{2^{N}}}{\ifelse{html}{\eqn{2^{N}}}{2^N}}, less when
interactions are present, see below).

The fitted model objects are not stored in the result. To get (a subset of)
the models, use \lcode{get.models} on the object returned by \code{dredge}.
Another way to get all the models is to run 
\code{lapply(dredge(..., evaluate = FALSE), eval)}, 
which avoids fitting models twice.

For a list of model types that can be used as a \code{global.model} see
\link[=MuMIn-models]{the list of supported models}. Modelling functions that 
do not store a \code{call} in their result should be run \emph{via} a 
wrapper function created by \lcode{updateable}.

\subsection{Information criterion}{
\code{rank} is found by a call to \code{match.fun} and may be specified as a
function, a symbol, or as a character string specifying a function to be searched 
for from the environment of the call to \code{dredge}. It can be also a 
one-element named list, where the first element is taken as the rank function.
The function \code{rank} must accept a model object as its first argument and
always return a scalar.
}

\subsection{Interactions}{
By default, marginality constraints are respected, so that \dQuote{all possible
combinations} include only those that contain interactions with their
respective main effects and all lower order terms,  unless the
\code{global.model} makes an exception to this principle (e.g. due
to a nested design such as \code{a / b}).

}

\subsection{Subsetting}{

The resulting set of models can be constrained with three methods:
(1) set limits on
the number of terms in a model with \code{m.lim}, (2) bind term(s) to all
models with \code{fixed}, and (3) use \code{subset} for more complex rules.
To be included in the selection table, the formulation of a model must satisfy 
all these conditions.

%%Terms in \code{fixed} argument are applied before the combinations are
%%generated, therefore more efficient than \code{subset}.

\code{subset} can be an \emph{expression} or a \emph{matrix}.
If a matrix, it should be a logical, lower triangular matrix, with 
rows and columns corresponding to \code{global.model} terms. If this matrix has 
\code{dimnames}, they must match the term names (as returned by \code{getAllTerms}).
Unmatched names are silently ignored. Otherwise, if rows or columns are unnamed, 
they are matched positionally to the model terms, and \code{dim(subset)} must 
be equal to the number of terms. 
For example, \code{subset["a", "b"] == FALSE} excludes models 
with  both \var{a} and \var{b} terms; and if unnamed, \code{subset}, 
\code{subset[2, 3] == FALSE} will prevent the second and third terms of the 
global model from being both in the same model. \cr

\code{demo(dredge.subset)} has examples of using the \code{subset} matrix in 
conjunction with correlation matrices to exclude models containing collinear 
predictors.

In the form of an \code{expression}, the argument \code{subset} acts similarly
to that of \code{subset()} for \code{data.frame}s. Model terms can be referred 
to by name as variables in the expression, except that they are interpreted as 
logical values indicating the presence of a term in the model.

The expression can contain any of the \code{global.model} term names, as well as
 names of the \code{varying} list items. \code{global.model} term names take
precedence when identical to names of \code{varying}, so to avoid ambiguity
\code{varying} variables in \code{subset} expression should be enclosed in
\code{V()} (e.g. \code{V(family) == "Gamma"}) assuming that
\code{varying} is something like \code{list(family =} \code{c("Gamma", ...))}).

If elements of \code{varying} are unnamed, they are coerced into names. Calls 
and symbols are represented as character values (via "deparse"), and everything 
except numeric, logical, character and NULL values is represented by element 
numbers (e.g. \code{subset = V(family) == 2} points to \code{Gamma} family in 
\code{varying =}\code{list(family =}\code{list(gaussian, Gamma)}). This can 
easily become obscure, so using named lists in \code{varying} is recommended.
Examples can be found in \code{demo(dredge.varying)}.

Term names appearing in \code{fixed} and \code{subset} must be given exactly 
as they are returned by \code{getAllTerms(global.model)}, which may differ
from the original term names (e.g. the interaction term components are ordered
alphabetically).

The \code{with(x)} and \code{with(+x)} notation indicates, respectively, any and
all interactions including the main effect term \code{x}. This is only effective
with marginality exceptions. The extended form \code{with(x, order)} allows to
specify the order of interaction of terms of which \code{x} is a part. For
instance, \code{with(b, 2:3)} selects models with at least one second- or
third-order interaction of variable \code{b}. The second (positional)
argument is coerced to an integer vector. The \dQuote{dot} notation \code{.(x)} is
an alias for \code{with}.

The special variable \ifelse{latex}{\bq{*nvar*}}{\code{`*nvar*`}} 
(backtick-quoted), in the \code{subset} expression is equal to the number of 
terms in the model (\bold{not} the number of parameters).

To include a model term conditionally on the presence of another term, 
use \code{dc} (\dQuote{\bold{d}ependency \bold{c}hain}) in
the \code{subset} expression. \code{dc} takes any number of term names as
arguments, and allows a term to be included only if all preceding ones
are also present (e.g. \code{subset = dc(a, b, c)} allows for models \code{a},
\code{a+b} and \code{a+b+c} but not \code{b}, \code{c}, \code{b+c} or
\code{a+c}).

\code{subset} expression can have a form of an unevaluated \code{call},
\code{expression} object, or a one-sided \code{formula}. See \sQuote{Examples}.

Compound model terms (such as interactions, \sQuote{as-is} expressions within
\code{I()} or smooths in \code{gam}) should be enclosed within curly brackets
 (e.g. \code{{s(x,k=2)}}), or \link[=Quotes]{backticks} (like non-syntactic
 names, e.g. \ifelse{latex}{
 \bq{s(x, k = 2)}
 }{
 \code{`s(x, k = 2)`}
 }), except when they are arguments to \code{with} or \code{dc}.
 Backtick-quoted names must match exactly (including whitespace) the term names
 as returned by \code{getAllTerms}.

\subsection{\code{subset} expression syntax summary}{

	\describe{
	\item{\code{a & b}}{ indicates that model terms \var{a} and \var{b} must be
        present (see \link[=Logic]{Logical Operators}) }
    \item{\code{{log(x,2)}} or \bq{log(x, 2)}}{ represent a complex
		model term \code{log(x, 2)}}
	\item{\code{V(x)}}{ represents a \code{varying} item \var{x} }
	\item{\code{with(x)}}{ indicates that at least one term containing the main 
	    effect term \var{x} must be present }
	\item{\code{with(+x)}}{ indicates that all the terms containing the main 
	    effect term \var{x} must be present }
	\item{\code{with(x, n:m)}}{ indicates that at least one term containing an 
	    \var{n}-th to \var{m}-th order interaction term of \var{x} must be present }
	\item{\code{dc(a, b, c,...)}}{ \sQuote{dependency chain}: \var{b} is allowed only
		if \var{a} is present, and \var{c} only if both \var{a} and \var{b} are
		present, etc. }
	\item{\code{`*nvar*`}}{ the number of terms in the model. }
	}
}

To simply keep certain terms in all models, it is much more efficient to use 
the \code{fixed} argument. The \code{fixed} formula is interpreted in the same
 manner as model formula, so the terms must not be quoted.
}

\subsection{Missing values}{
Use of \code{na.action = "na.omit"} (\R's default) or \code{"na.exclude"}  in
\code{global.model} must be avoided, as it results with sub-models fitted to
different data sets if there are missing values. An error is thrown if it is
detected.

It is a common mistake to give \code{na.action} as an argument in the call
  to \code{dredge} (typically resulting in an error from the \code{rank}
  function to which the argument is passed through \sQuote{\dots}), while the 
  correct way
  is either to pass \code{na.action} in the call to the global model or to set
  it as a \link[=options]{global option}.
}

\subsection{Intercept}{ 
If present in the \code{global.model}, the intercept will be included in all 
sub-models.
}

\subsection{Methods}{
There are \lxcode{subset}{=subset.model.selection} and
\lxcode{plot}{=plot.model.selection} methods, the latter creates a
graphical representation of model weights and per-model term sum of weights.
Coefficients can be extracted with \code{coef} or \lcode{coefTable}.

}
}


\value{
An object of class \code{c("model.selection", "data.frame")}, being a
\code{data.frame}, where each row represents one model.
See \lcode{model.selection.object} for its structure.
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

\lcode{get.models}, \lcode{model.avg}. \lcode{model.sel} for
manual model selection tables.

Possible alternatives: \code{glmulti} in package \pkg{glmulti}
and \code{bestglm} (\pkg{bestglm}).
%% \lxcode{regsubsets}{leaps} in package \pkg{leaps} also performs all-subsets
\code{regsubsets} in package \pkg{leaps} also performs all-subsets
regression.

Variable selection through regularization provided by various packages,
    e.g. \pkg{glmnet}, \pkg{lars} or \pkg{glmmLasso}.
}


\examples{
# Example from Burnham and Anderson (2002), page 100:

#  prevent fitting sub-models to different datasets
\dontshow{oop <- }
options(na.action = "na.fail")

fm1 <- lm(y ~ ., data = Cement)
dd <- dredge(fm1)
subset(dd, delta < 4)

# Visualize the model selection table:
\dontshow{ if(require(graphics)) \{ }
par(mar = c(3,5,6,4))
plot(dd, labAsExpr = TRUE)
\dontshow{ \} }

# Model average models with delta AICc < 4
model.avg(dd, subset = delta < 4)

#or as a 95\% confidence set:
model.avg(dd, subset = cumsum(weight) <= .95) # get averaged coefficients

#'Best' model
summary(get.models(dd, 1)[[1]])

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
# the same as above, without "dc"
dredge(fm1, subset = (X1 | !X2) && (X2 | !X3) && xor(X3, X4))

# Include only models with up to 2 terms (and intercept)
dredge(fm1, m.lim = c(0, 2))
}

# Add R^2 and F-statistics, use the 'extra' argument
dredge(fm1, m.lim = c(NA, 1), extra = c("R^2", F = function(x)
    summary(x)$fstatistic[[1]]))

# with summary statistics:
dredge(fm1, m.lim = c(NA, 1), extra = list(
    "R^2", "*" = function(x) {
        s <- summary(x)
        c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
            F = s$fstatistic[[1]])
    })
)

# Add other information criteria (but rank with AICc):
dredge(fm1, m.lim = c(NA, 1), extra = alist(AIC, BIC, ICOMP, Cp))
\dontshow{options(oop)}
}

\keyword{models}
