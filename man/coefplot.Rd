\encoding{utf-8}
\name{coefplot}
\alias{coefplot}
\alias{plot.averaging}
\title{Plot model coefficients}
\usage{
coefplot(
  x,
  lci,
  uci,
  labels = NULL,
  width = 0.15,
  shift = 0,
  horizontal = TRUE,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  xlim = NULL,
  ylim = NULL,
  labAsExpr = TRUE,
  mar.adj = TRUE,
  lab.line = 0.5,
  lty = par("lty"),
  lwd = par("lwd"),
  pch = 21,
  col = par("col"),
  bg = par("bg"),
  dotcex = par("cex"),
  dotcol = col,
  staplelty = lty,
  staplelwd = lwd,
  staplecol = col,
  zerolty = "dotted",
  zerolwd = lwd,
  zerocol = "gray",
  las = 2,
  ann = TRUE,
  axes = TRUE,
  add = FALSE,
  type = "p",
  ...
)

\method{plot}{averaging}(
  x,
  full = TRUE,
  level = 0.95,
  intercept = TRUE,
  parm = NULL,
  labels = NULL,
  width = 0.1,
  shift = max(0.2, width * 2.1 + 0.05),
  horizontal = TRUE,
  xlim = NULL,
  ylim = NULL,
  main = "Model-averaged coefficients",
  xlab = NULL,
  ylab = NULL,
  add = FALSE,
  ...
)
}
\arguments{
\item{x}{either a (possibly named) vector of coefficients (for \code{coefplot}),
or an \code{\link[=model.avg]{"averaging"}} object.}

\item{lci, uci}{vectors of lower and upper confidence intervals. Alternatively
a two-column matrix with columns containing confidence intervals, in
which case \code{uci} is ignored.}

\item{labels}{optional vector of coefficient names. By default, names of \code{x}
are used for labels.}

\item{width}{width of the staples (= end of whisker).}

\item{shift}{the amount of perpendicular shift for the dots and whiskers.
Useful when adding to an existing plot.}

\item{horizontal}{logical indicating if the plots should be horizontal;
defaults to \code{TRUE}.}

\item{main}{an overall title for the plot: see \code{title}.}

\item{xlab, ylab}{x- and y-axis annotation. Can be suppressed by \code{ann=FALSE}.}

\item{xlim, ylim}{optional, the \emph{x} and \emph{y} limits of the plot.}

\item{labAsExpr}{logical indicating whether the coefficient names should
be transformed to expressions to create prettier labels (see
\code{\link{plotmath}})}

\item{mar.adj}{logical indicating whether the (left or lower) margin should be
expanded to fit the labels}

\item{lab.line}{margin line for the labels}

\item{lty, lwd, pch, col, bg}{default line type, line width, point character,
foreground colour for all elements, and background colour for open symbols.}

\item{dotcex, dotcol}{dots point size expansion and color.}

\item{staplelty, staplelwd, staplecol}{staple line type, width, and colour.}

\item{zerolty, zerolwd, zerocol}{zero-line type, line width, colour.
Setting \code{zerolty}    to \code{NA} suppresses the line.}

\item{las}{the style of labels for coefficient names. See \code{\link{par}}.}

\item{ann}{\code{logical} indicating if axes should be annotated (by \code{xlab} and \code{ylab}).}

\item{axes}{a logical value indicating whether both axes should be drawn on
the plot.}

\item{add}{logical, if true \emph{add} to current plot.}

\item{type}{if \code{"n"}, the plot region is left empty, any other value causes
the plot being drawn.}

\item{...}{additional arguments passed to \code{coefplot} or more
\link[=par]{graphical parameters}.}

\item{full}{a logical value specifying whether the \dQuote{full}
model-averaged coefficients are plotted. If \code{FALSE}, the
\dQuote{subset}-averaged coefficients are plotted, and both types if
\code{NA}. See \code{\link{model.avg}}.}

\item{level}{the confidence level required.}

\item{intercept}{logical indicating if intercept should be included in the
plot}

\item{parm}{a specification of which parameters are to be plotted, either a
vector of numbers or a vector of names. If missing, all parameters are
considered.}
}
\value{
An invisible \code{matrix} containing coordinates of points and whiskers, or,
a two-element list of such, one for each coefficient type in
\code{plot.averaging} when \code{full} is \code{NA}.
}
\description{
Produce dot-and-whisker plot of the model(-averaged) coefficients,
with confidence intervals
}
\details{
Plot model(-averaged) coefficients with confidence intervals.
}
\examples{
fm <- glm(Prop ~ dose + I(dose^2) + log(dose) + I(log(dose)^2),
  data = Beetle, family = binomial, na.action = na.fail)
ma <- model.avg(dredge(fm))

# default coefficient plot:
plot(ma, full = NA, intercept = FALSE)

# Add colours per coefficient type
# Replicate each colour n(=number of coefficients) times
clr <- c("black", "red2")
i <- rep(1:2, each = length(coef(ma)) - 1)
plot(ma, full = NA, intercept = FALSE,
   pch = 22, dotcex = 1.5,
   col = clr[i], bg = clr[i],
   lwd = 6, lend = 1, width = 0, horizontal = 0)

# Use `type = "n"` and `add` argument to e.g. add grid beneath the figure
plot(ma, full = NA, intercept = FALSE,
   width = 0, horizontal = FALSE, zerolty = NA, type = "n")
grid()
plot(ma, full = NA, intercept = FALSE,
   pch = 22, dotcex = 1.5,
   col = clr[i], bg = clr[i],
   lwd = 6, lend = 1, width = 0, horizontal = FALSE, add = TRUE)

}
\author{
Kamil Barto\enc{ń}{n}
}
\keyword{hplot}
