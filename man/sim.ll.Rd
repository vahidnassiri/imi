\name{sim.ll}
\alias{sim.ll}
\docType{data}
\title{
sim.ll}
\description{
This dataset is a list contains the simulation results from the paper 'Iterative multiple imputation: a framework to determine the number of imputed datasets' for the scenario where \eqn{\sigma^2=64} and \eqn{\rho=0.8}. It contains the averaged distance between two steps (over 100 replications) for \eqn{M=500}. Also the selected final M for different \eqn{\epsilon}'s and \eqn{k_0}'s.}
\usage{data("sim.ll")}
\format{
 The format is: List of 3 lists: the first list contains the averaged simulation results over 100 replications for parameter estimates, the second list contains the same but for the p-values of one sample and paired t-test, and the third list contains the final selected M for eahc scenario.

}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(sim.ll)
}
\keyword{datasets}
