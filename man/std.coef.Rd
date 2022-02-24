\name{std.coef}
\alias{std.coef}
\alias{beta.weights}
\alias{partial.sd}

\encoding{utf-8}
\newcommand{\myequation}{\ifelse{latex}{\eqn{#1}{#3}}{\ifelse{html}{\enc{\eqn{#2}}{#3}}{#3}}}

\title{Standardized model coefficients}
\description{
Standardize model coefficients by Standard Deviation or Partial Standard Deviation.
}

\usage{

std.coef(x, partial.sd, ...)

partial.sd(x)

# Deprecated:
beta.weights(model)
}
\arguments{
  \item{x, model}{a fitted model object. }
  \item{partial.sd}{logical, if set to \code{TRUE}, model coefficients are
	multiplied by partial \acronym{SD}, otherwise they are multiplied by the
	ratio of the standard deviations of the independent variable and dependent
	variable. }
 \item{\dots}{additional arguments passed to \code{\link{coefTable}}, e.g.
    \code{dispersion}. }
}

\details{
	Standardizing model coefficients has the same effect as centring and
	scaling the input variables. \dQuote{Classical} standardized coefficients
	are calculated as 
    \myequation{\beta^{*}_i = \beta_i\frac{s_{X_{i}}}{s_{y}}
        }{\betaᵢ* = \betaᵢ (sₓᵢ / Sᵧ)
        }{Bi* = Bi * (SXi / Sy)}
    , where
	\myequation{\beta}{\beta}{B} is the unstandardized coefficient, 
    \myequation{s_{X_{i}}}{sₓᵢ}{SXi}    
    is the standard deviation of associated dependent variable 
    \myequation{X_i}{Xᵢ}{Xi}  and
	\myequation{s_{y}}{Sᵧ}{Sy}
    is \acronym{SD} of the response variable.
    

	If variables are intercorrelated, the standard deviation of
    \myequation{X_i}{Xᵢ}{Xi}
	used in computing the standardized coefficients 
    \myequation{\beta_i^{*}}{\betaᵢ*}{Bi*} should be
	replaced by the partial standard deviation of 
    \myequation{X_i}{Xᵢ}{Xi} which is adjusted for
	the multiple correlation of 
    \myequation{X_i}{Xᵢ}{Xi} with the other \eqn{X} variables
	included in the regression equation. The partial standard deviation is
	calculated as 
    \myequation{s_{X_{i}}^{*}=s_{X_{i}} \var{VIF}(X_i)^{-0.5} (\frac{n-1}{n-p} )^{0.5}
        }{s*ₓᵢ = sₓᵢ √(VIFₓᵢ⁻¹) √((n-1)/(n-p))
        }{sXi* = SXi * sqrt(1/VIF(Xi)) * sqrt((n-1)/(n-p))}, 
    where \var{VIF} is the variance inflation factor,
	\var{n} is the number of observations and \var{p}, the number of predictors in
	the model. The coefficient is then transformed as 
    \myequation{\beta^{*}_i = \beta_i s_{X_{i}}^{*}
        }{\betaᵢ* = \betaᵢ s*ₓᵢ
        }{Bi* = Bi * SXi*}.

}

\value{
    A matrix with at least two columns for the standardized coefficient estimate
    and its standard error. Optionally, the third column holds degrees of freedom
    associated with the coefficients.
}


\author{Kamil Barto\enc{ń}{n}. Variance inflation factors calculation is based
on function \code{vif} from package \pkg{car} written by Henric Nilsson and John
Fox.
}

\references{

Cade, B.S. (2015) Model averaging and muddled multimodel inferences. 
    \emph{Ecology} 96, 2370-2382.

Afifi A., May S., Clark V.A. (2011) \emph{Practical Multivariate Analysis},
Fifth Edition. CRC Press.

Bring, J. (1994). How to standardize regression coefficients. \emph{The American
Statistician} 48, 209-213.
}

\seealso{
\code{partial.sd} can be used with \code{\link{stdize}}.

\code{\link{coef}} or \code{\link{coeffs}} and \code{\link{coefTable}} for
unstandardized coefficients.
}
\examples{

# Fit model to original data:
fm  <- lm(y ~ x1 + x2 + x3 + x4, data = GPA)

# Partial SD for the default formula:   y ~ x1 + x2 + x3 + x4
psd <- partial.sd(lm(data = GPA))[-1] # remove first element for intercept

# Standardize data:
zGPA <- stdize(GPA, scale = c(NA, psd), center = TRUE)
# Note: first element of 'scale' is set to NA to ignore the first column 'y'

# Coefficients of a model fitted to standardized data:
zapsmall(coefTable(stdizeFit(fm, newdata = zGPA)))
# Standardized coefficients of a model fitted to original data:
zapsmall(std.coef(fm, partial = TRUE))


# Standardizing nonlinear models:
fam <- Gamma("inverse")
fmg <- glm(log(y) ~ x1 + x2 + x3 + x4, data = GPA, family = fam)

psdg <- partial.sd(fmg)
zGPA <- stdize(GPA, scale = c(NA, psdg[-1]), center = FALSE)
fmgz <- glm(log(y) ~ z.x1 + z.x2 + z.x3 + z.x4, zGPA, family = fam)

# Coefficients using standardized data:
coef(fmgz) # (intercept is unchanged because the variables haven't been
           #  centred)
# Standardized coefficients:
coef(fmg) * psdg

}


\keyword{models}
