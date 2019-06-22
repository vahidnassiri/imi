\name{imi.make.plots}
\alias{imi.make.plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
imi.make.plots}
\description{
This function take the simulation output available via data of type list with the names: sim.ss, sim.sm, sim.sl, sim.ms, sim.mm, sim.ml, sim.ls, sim.lm, and sim.ll as the input and produce the convergence plots for estimated parameters as well as the p-values of the t-tests. These are similar to the Figure 1 and Figure 2 in the paper 'Iterative multiple imputation:	a framework to determine the number of imputed datasets'
}
\usage{
imi.make.plots(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list which should be selected from the available datasets in the package: sim.ss, sim.sm, sim.sl, sim.ms, sim.mm, sim.ml, sim.ls, sim.lm, and sim.ll
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
sim.ss, sim.sm, sim.sl, sim.ms, sim.mm, sim.ml, sim.ls, sim.lm, and sim.ll }
\examples{
imi.make.plots(sim.ss)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
