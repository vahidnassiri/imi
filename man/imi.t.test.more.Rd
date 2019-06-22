\name{imi.t.test.more}
\alias{imi.t.test.more}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
imi.t.test.more}
\description{
This function uses the function mi.t.test in R package MKmisc to perform one or two sample Student's t-test on an incomplete dataset using multiple imputation and determines the sufficient number of imputations using iterative multiple imputation (imi) procedure.
}
\usage{
imi.t.test(data.miss,data.imp0,max.M=500,epsilon,method='pmm',
                       x, y = NULL, alternative='two.sided',mu,paired = FALSE, var.equal = FALSE,
                       conf.level = 0.95,conv.plot=TRUE,successive.valid=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data.miss}{A data frame with the variable in the model as its columns. Note that the missing values should be indicated by NA.}
    \item{data.imp0}{A list with already imputed sets of data as its components.}
   \item{max.M}{The maximum number of iterations which the algorithm should terminate afterwards in case of non-covergence.}
  \item{epsilon}{The threshold for difference between two iterations.}
  \item{method}{Specifying string value 'mvn' would impute the data using a multivariate normal predictive model in R package amelia2, any other specification will impute the data using fully conditional specification approach in R package mice. One can see the method in documentation of function mice in R package mice. Specifying 'auto' will selected the predictive model based on the measurement level of each variable.}
  \item{x}{Name of the data column which the t-test should be performed on it.}
  \item{y}{In case of two sample t-test this input variable specifiess the name of the second column for paired-test, or an indicator variable showing different populations for indepdendent two sample t-test.}
  \item{alternative}{A character string specifying the alternative hypothesis, it takes the values "two.sided" (default), "greater" or "less".}
 \item{mu}{A real number specifying the test value.}
  \item{paired}{A logical value specifying whether a paired t-test should be performed or not.}
  \item{var.equal}{A logical value specifying whether variance equality should be assumed or not.}
    \item{conf.level}{Confidence level of the test.}
  \item{conv.plot}{A logitical value, if TRUE then a convergence plot will be generated, if FALSE no plot will be provided.}
  \item{successive.valid}{An integer with minimum 1 which specifies the number of successive steps the stopping rule should be validated so the procedure could terminate.}
   \item{print.progress}{A logical variable, if TRUE it prints the progress of imputation.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{test.result }{A list with the final MI-based outcome of the t-test.}
  \item{data.imp }{A list with imputed datasets as its components.}
  \item{dis.steps }{A vector with computed distance between iterations.}
  \item{conv.status }{If 1 then convergence is achieved, if 0 with max.M iterations, still the convergence could not be achieved.}
  \item{M }{The selected number of imputed datasets.}
}
\references{
https://stat.ethz.ch/R-manual/R-devel/library/stats/html/t.test.html}
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
\examples{
# To illustrate the use of this function we use the cholestrol dataset in R package norm2.
> library(norm2)
> data(cholesterol)
>
> chol.ttest=imi.t.test (cholesterol,M0='manual',max.M=500,epsilon=0.05/10,method='mvn',
+                           x=names(cholesterol)[3], y = NULL, alternative='two.sided',mu=220,
+                           paired = FALSE, var.equal = FALSE,
+                           conf.level = 0.95,conv.plot=TRUE,successive.valid='manual')
-- Imputation 1 --

  1  2  3  4  5  6

[1] "The time it takes (in seconds) to imput the data once and fit the model to it is: 0.01"
What is your choice of initial number of imputations?2
What is your choice for successive steps validation?3
-- Imputation 1 --

  1  2  3  4  5  6  7  8

-- Imputation 1 --

  1  2  3  4  5  6  7

# it continues till the convergence happens.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
