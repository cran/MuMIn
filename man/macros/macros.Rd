\newcommand{\href}{#2}
\newcommand{\bq}{\verb{`}\code{#1}\verb{`}}

\newcommand{\mydequation}{\ifelse{latex}{\deqn{#1}{#3}}{\ifelse{html}{\enc{\deqn{#1}}{#3}}{#3}}}
\newcommand{\myequation}{\ifelse{latex}{\eqn{#1}{#3}}{\ifelse{html}{\enc{\eqn{#1}}{#3}}{#3}}}

%% \newcommand{\CRANpkg}{\href{https://cran.r-project.org/package=#1}{\pkg{#1}}} 
\newcommand{\fckCRANpkg}{\bold{#1}}
\newcommand{\fckpkg}{\bold{#1}}
\newcommand{\lxcode}{\code{\link[#2]{#1}}}
\newcommand{\lcode}{\code{\link{#1}}}

\newcommand{\logLik}{\myequation{\log\mathcal{L}(#1)}{XXXX}{logL(#1)}}
\newcommand{\Rsq}{\myequation{R^{2}}{XXXX}{R^2}}
\newcommand{\Rsq}{\myequation{R^{2}}{XXXX}{R^2}}
%% \newcommand{\Rsq}{\ifelse{latex}{\eqn{R^{2}}{R^2}}{\ifelse{html}{\enc{R²}{R^2}}{R^2}}}
\newcommand{\Rsqx}{\myequation{R_{#1}^{2}}{XXXX}{R_#1^2}}
