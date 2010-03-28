\name{ddf}
\alias{ddf}
\title{Linear Programming for the Directional Distance Function}
\description{Solve the Additive Model under the VRS assumption}
\usage{
ddf(base = NULL, frontier = NULL, noutput = 1, direction = NULL)
}
\arguments{
  \item{base}{A data set for DMUs to be evaluated. A data frame with
    J1*(M+N) dimention, where J1 is the number of DMUs, M for the number
    of inputs, and N for the number of outputs.}
  \item{frontier}{A data set for DMUs to be used in constructing a
    production possibility set (PPS). A data frame with J2*(M+N)
    dimention, where J2 is the number of DMUs, M for the number of
    inputs, and N for the number of outputs.}
  \item{noutput}{The number of outputs (N).}
  \item{direction}{A directional vector for inputs and outputs.}
}
\details{
  The DDF under the VRS assumption is calculated. For model
  specification, take a look at Cooper et al. (2007).
}

\value{
A data frame with J1*(J1+M+N), of which has efficiency scores, optimal
weightes and optimal slacks. Take a look at the example below.
}

\references{
  Cooper, W., Seiford, L. and Tone, K. (2007). \emph{Data envelopment
  analysis: a comprehensive text with models, applications, references
  and DEA-solver software} (2nd ed.). Springer Verlag, New York.

  Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
    Envelopment Analysis}. Press (in Korean). 
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{direc.dea}}}

\examples{
## Simple Example of one input and one output.
my.dat <- data.frame(y = c(1, 2, 4, 6, 7, 9, 9),
                     x = c(3, 2, 6, 4, 8, 8, 10))
(re <- ddf(my.dat, noutput = 1, direction = c(1, 1)))
}

\keyword{Data Envelopment Analysis} 