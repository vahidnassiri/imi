\name{imi.t.test.conv}
\alias{imi.t.test.conv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
imi.t.test.conv}
\description{
This function uses iterative multiple imputation idea to compute the distance between p-values of t-tests from imputed datasets to see whether the already generated sets of imputed data are sufficient for such analysis or not.}
\usage{
imi.glm.conv(data.imp,epsilon,x, y = NULL, alternative='two.sided',
              mu,paired = FALSE, var.equal = FALSE, conf.level = 0.95,
              conv.plot=TRUE,successive.valid=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data.imp}{A list imputed sets of data as its components.}
 \item{epsilon}{The threshold for difference between two iterations.}
  \item{x}{Name of the data column which the t-test should be performed on it.}
  \item{y}{In case of two sample t-test this input variable specifiess the name of the second column for paired-test, or an indicator variable showing different populations for indepdendent two sample t-test.}
  \item{alternative}{A character string specifying the alternative hypothesis, it takes the values "two.sided" (default), "greater" or "less".}
 \item{mu}{A real number specifying the test value.}
  \item{paired}{A logical value specifying whether a paired t-test should be performed or not.}
  \item{var.equal}{A logical value specifying whether variance equality should be assumed or not.}
    \item{conf.level}{Confidence level of the test.}
  \item{conv.plot}{A logitical value, if TRUE then a convergence plot will be generated, if FALSE no plot will be provided.}
  \item{successive.valid}{An integer with minimum 1 which specifies the number of successive steps the stopping rule should be validated so the procedure could terminate.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{dis.steps }{A vector with computed distance between iterations.}
 \item{sufficient.M }{An integer indicating the minumum number of sufficent imputed datasets, in case of insufficiency it will take the value 'Not sufficient!'}

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
imi.t.test}
\examples{
# In th setting of the example of imi.t.test
imi.t.test.conv (chol.ttest$data.imp,epsilon=0.05/10,x=names(cholesterol)[3], 
+     y = NULL, alternative='two.sided',mu,paired = FALSE, var.equal = FALSE,
+                        conf.level = 0.95,conv.plot=TRUE,successive.valid=3)
$dis.steps
 [1] 0.3117110189 0.1578016324 0.0295781152 0.0487635084 0.0288667328 0.0318036141 0.0173372375 
 [8] 0.0090263689 0.0575828782
[10] 0.0007220647 0.0145497285 0.0417895357 0.0044211378 0.0124385165 0.0203312525 
[16] 0.0050927153 0.0058802488 0.0117815975
[19] 0.0022407713 0.0026592932 0.0007328110

$sufficient.M
[1] 22

>
> imi.t.test.conv (chol.ttest$data.imp,epsilon=0.05/100,x=names(cholesterol)[3], 
+        y = NULL, alternative='two.sided',mu,paired = FALSE, var.equal = FALSE,
+                  conf.level = 0.95,conv.plot=TRUE,successive.valid=3)
$dis.steps
 [1] 0.3117110189 0.1578016324 0.0295781152 0.0487635084 0.0288667328 0.0318036141 0.0173372375
 [8] 0.0090263689 0.0575828782
[10] 0.0007220647 0.0145497285 0.0417895357 0.0044211378 0.0124385165 0.0203312525 0.0050927153 
[17] 0.0058802488 0.0117815975
[19] 0.0022407713 0.0026592932 0.0007328110

$sufficient.M
[1] "Not sufficient!"

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
