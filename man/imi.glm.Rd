\name{imi.glm}
\alias{imi.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
imi.glm
}
\description{
This function fits a generalized linear model for a given set of predictors and a response variable to an incomplete dataset using multiple imputation and determines the sufficient number of imputed datasets using iterative multiple imputation (imi) procedure.
}
\usage{
imi.glm(data.miss,family=binomial(link='logit'),M0=5,max.M=500,epsilon,method='pmm',
          resp,regressors,conv.plot=TRUE,dis.method='mahalanobis',
          mah.scale='combined',successive.valid=3,max.iter.glm=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.miss}{A data frame with the variable in the model as its columns. Note that the missing values should be indicated by NA.}
  \item{family}{Indicating the error distribution and the link function (see documentation of function glm in R stats package).}
  \item{M0}{The initial number of imputations, it should be an integer with minimum 2. It also can take the string value 'auto' which lets the user to decide about the initial number of imputation after the software reports the time it takes to impute the data and fit the model to it twice.}
  \item{max.M}{The maximum number of iterations which the algorithm should terminate afterwards in case of non-covergence.}
  \item{epsilon}{The threshold for difference between two iterations.}
  \item{method}{Specifying string value 'mvn' would impute the data using a multivariate normal predictive model in R package amelia2, any other specification will impute the data using fully conditional specification approach in R package mice. One can see the method in documentation of function mice in R package mice. Specifying 'auto' will selected the predictive model based on the measurement level of each variable.}
  \item{resp}{A string value input with the name of the response variable. Note that this should match the name of one of the columns in the data.miss.}
  \item{regressors}{A vector of string values with the names of the predictors. Note that they should match the names of the variables in data.miss.}
  \item{conv.plot}{A logitical value, if TRUE then a convergence plot will be generated, if FALSE no plot will be provided.}
  \item{dis.method}{A string takes its values among 'euclidean', 'inf.norm', and 'mahalanobis' which specifies the distance measure between two iterations. Note that our suggestion is to use 'mahalanobis', other options are provided for research purposes.}
   \item{mah.scale}{A string takes its values among 'within', 'between', and 'combined' which specifies the scale matrix in Mahalanobis distance. Note that our suggestion is to use 'combined', other options are provided for research purposes.}
  \item{successive.valid}{An integer with minimum 1 which specifies the number of successive steps the stopping rule should be validated so the procedure could terminate.}
  \item{max.iter.glm}{The maximum number of iterations for the glm algorithm.}
   \item{print.progress}{A logical variable, if TRUE it prints the progress of imputation.}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{mi.param }{A list with the final MI-based estimated model parameters, their covariance matrix, as well as the within and between imputation covariance matrices.}
  \item{data.imp }{A list with imputed datasets as its components.}
  \item{dis.steps }{A vector with computed distance between iterations.}
  \item{conv.status }{If 1 then convergence is achieved, if 0 with max.M iterations, still the convergence could not be achieved.}
  \item{M }{The selected number of imputed datasets.}

%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
https://www.rdocumentation.org/packages/mice/versions/2.30/topics/mice

https://cran.r-project.org/web/packages/Amelia/
}
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
> # specifying sample size and number of predictors
> sample.size=100
> num.var.log=2
> # creating a correlated set of predictors
> x.orig=matrix(rnorm(num.var.log*sample.size),sample.size,num.var.log)
> x=cbind(scale(x.orig,
+               center=TRUE,scale=FALSE))
> # creating the compound-symmetry structured covariance matrix
> sigma2=4
> tau=1
> cov.mat=diag(sigma2,num.var.log)+tau
> # making the data correlated
> chol.cov=chol(cov.mat)
> for (i in 1:sample.size){
+   x[i,]=t(chol.cov)%*%x[i,]
+ }
> x=t(t(x)+apply(x.orig,2,mean))
> # specifying model parameters
> beta=c( 0.2, -1.0,  0.5)
> z=beta[1]+(x%*%beta[2:3])
> # computing inverse logit transformation
> pr = 1/(1+exp(-z))
> # generating response variable
> y = rbinom(sample.size,1,pr)
> # creating complete the data
> data = data.frame(y,x)
> # specifying the regressors and predictors
> resp='y'
> regressors=c('X1','X2')
> # creating missing values in the dataset
> require('mice')
> x.miss=ampute(x,prop=0.1,mech='MAR')$amp
> data.miss=data.frame(y,x.miss)
> # Determining number of imputations, impute the incomplete data and fit the model to it
> out.glm=imi.glm (data.miss,family=binomial(link='logit'),M0='manual',max.M=500,
+          epsilon=0.05,method='mvn',resp,regressors, conv.plot=TRUE,
+          dis.method='mahalanobis',mah.scale='within',
+          successive.valid='manual',max.iter.glm=1000)
[1] "The time it takes (in seconds) to imput the data two times and fit the model to them is: 0.02"
What is your choice of initial number of imputations?2
What is your choice for successive steps validation?3
We are working on the initial imputation.
We are working on imputed dataset number  3
We are working on imputed dataset number  4
We are working on imputed dataset number  5
We are working on imputed dataset number  6
We are working on imputed dataset number  7
We are working on imputed dataset number  8
We are working on imputed dataset number  9
We are working on imputed dataset number  10
We are working on imputed dataset number  11
We are working on imputed dataset number  12
We are working on imputed dataset number  13
We are working on imputed dataset number  14
We are done! The convergence is acheived and the sufficient number of imputations is  14
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
