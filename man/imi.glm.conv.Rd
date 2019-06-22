\name{imi.glm.conv}
\alias{imi.glm.conv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
imi.glm.conv}
\description{
This function uses iterative multiple imputation idea to compute the distance between estimated parameters and their precisions from imputed datasets when fitting a generalized linear model to see whether the already generated sets of imputed data are sufficient for such analysis or not.
}
\usage{
imi.glm.conv(data.imp,family=binomial(link='logit'),epsilon,resp,regressors,
                         conv.plot=TRUE, dis.method='mahalanobis',
                         mah.scale='combined',successive.valid=3,max.iter.glm=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data.imp}{A list imputed sets of data as its components.}
  \item{family}{Indicating the error distribution and the link function (see documentation of function glm in R stats package.)}
   \item{epsilon}{The threshold for difference between two iterations.}
   \item{resp}{A string value input with the name of the response variable. Note that this should match the name of one of the columns in the data.miss.}
  \item{regressors}{A vector of string values with the names of the predictors. Note that they should match the names of the variables in data.miss.}
  \item{conv.plot}{A logitical value, if 'TRUE' then a convergence plot will be generated, if 'FALSE' no plot will be provided.}
  \item{dis.method}{A string takes its values among 'euclidean', 'inf.norm', and 'mahalanobis' which specifies the distance measure between two iterations. Note that our suggestion is to use 'mahalanobis', other options are provided for research purposes.}
   \item{mah.scale}{A string takes its values among 'within', 'between', and 'combined' which specifies the scale matrix in Mahalanobis distance. Note that our suggestion is to use 'combined', other options are provided for research purposes.}
  \item{successive.valid}{An integer with minimum 1 which specifies the number of successive steps the stopping rule should be validated so the procedure could terminate.}
  \item{max.iter.glm}{The maximum number of iterations for the glm algorithm.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{dis.steps }{A vector with computed distance between iterations.}
 \item{sufficient.M }{An integer indicating the minumum number of sufficent imputed datasets, in case of insufficiency it will take the value 'Not sufficient!'}

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
imi.glm
}
\examples{
# In the setting of the example for imi.glm, we have:
> imi.glm.conv(out.glm$data.imp,family=binomial(link='logit'),epsilon=0.05,resp,regressors,
+              conv.plot=TRUE,
+              dis.method='mahalanobis',mah.scale='within',successive.valid=3,max.iter.glm=1000)
$dis.steps
 [1] 0.46355012 0.20284904 0.14844815 0.02788952 0.09077971 0.12970693 0.04556116 0.02601013 0.05598301
[10] 0.06236743 0.04832095 0.02365007 0.01836882

$sufficient.M
[1] 12

> imi.glm.conv(out.glm$data.imp,family=binomial(link='logit'),epsilon=0.01,resp,regressors,
+              conv.plot=TRUE,
+              dis.method='mahalanobis',mah.scale='within',successive.valid=3,max.iter.glm=1000)
$dis.steps
 [1] 0.46355012 0.20284904 0.14844815 0.02788952 0.09077971 0.12970693 0.04556116 0.02601013 0.05598301
[10] 0.06236743 0.04832095 0.02365007 0.01836882

$sufficient.M
[1] "Not sufficient!"

>
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
