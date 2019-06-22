\name{combine.mi}
\alias{combine.mi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{combine.mi
%%  ~~function to do ... ~~
}
\description{
This function combines vectors of estimated parameters and theircovariance matrices from different imputed sets of data using Rubin's rule (Rubin, 1978).
}
\usage{
combine.mi(mi.param.est,mi.param.cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mi.param.est}{
A matrix of size number of imputed datasets by number of parameters in the model. Each row of this matrix should present the estimated parameters from each imputed dataset.
}
  \item{mi.param.cov}{
A list of length number of imputed datasets. Each element of this list is a covariance matrix of the size of the number of parameters in the model.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{param.est }{The combined vector of estimated parameters.}
 \item{param.cov }{The combined covariance matrix of the estimated parameters.}
  \item{between.cov }{The between-imputation covariance matrix.}
 \item{within.cov}{The within-imputation covariance matrix.}

%% ...
}
\references{
Rubin, Donald B. Multiple imputation for nonresponse in surveys. John Wiley & Sons, 1978.}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (x)
%{
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
