\newcommand{\href}{#2}
\newcommand{\bq}{\verb{`}\code{#1}\verb{`}}

\newcommand{\mydequation}{\ifelse{latex}{\deqn{#1}{#3}}{\ifelse{html}{\enc{\deqn{#1}}{#3}}{#3}}}
\newcommand{\myequation}{\ifelse{latex}{\eqn{#1}{#3}}{\ifelse{html}{\enc{\eqn{#1}}{#3}}{#3}}}

%% \newcommand{\CRANpkg}{\href{https://cran.r-project.org/package=#1}{\pkg{#1}}} 
\newcommand{\unCRANpkg}{\code{#1}}

\newcommand{\lcode}{\code{\link{#1}}}

\newcommand{\lxcode}{\code{\link[#2]{#1}}}
\newcommand{\unlxcode}{\code{#1}\ifelse{FALSE}{#2}{}}



\newcommand{\logLik}{\myequation{\log\mathcal{L}(#1)}{XXXX}{logL(#1)}}
\newcommand{\Rsq}{\myequation{R^{2}}{XXXX}{R^2}}
\newcommand{\Rsq}{\myequation{R^{2}}{XXXX}{R^2}}
\newcommand{\Rsqx}{\myequation{R_{#1}^{2}}{XXXX}{R_#1^2}}


\newcommand{\AICc}{\acronym{AIC\eqn{_{c}}{c}}}
\newcommand{\QAICc}{\acronym{QAIC\eqn{_{c}}{c}}}
\newcommand{\QICu}{\acronym{QIC\eqn{_{u}}{u}}}
\newcommand{\AIC}{\acronym{AIC}}
\newcommand{\QIC}{\acronym{QIC}}
\newcommand{\QAIC}{\acronym{QAIC}}
\newcommand{\BIC}{\acronym{BIC}}
\newcommand{\IC}{\acronym{IC}}

\newcommand{\DeltaIC}{\myequation{\Delta_{#1}}{}{Delta_#1}}

