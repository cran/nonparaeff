\name{cost.dea}
\alias{cost.dea}
\title{Linear Programming for Cost Minimization}
\description{Solve the Cost Minimization Probem with Given Input Prices}
\usage{
cost.dea(base = NULL, frontier = NULL, noutput = 1, input.price = NULL)
}
\arguments{
  \item{base}{A data set for DMUs to be evaluated. A data frame with
    J1*(M+N) dimention, where J1 is the number of DMUs, M for the number
    of inputs, and N for the number of outputs.}
  \item{frontier}{A data set for DMUs to be used in constructing a
    production possibility set (PPS). A data frame with J2*(M+N)
    dimention, where J2 is the number of DMUs, M for the number of
    inputs, and N for the number of outputs.}
  \item{noutput}{The number of outputs (M).}
  \item{input.price}{A vector for market prices of input factors.}
}
\details{
  The cost minimization problem under the CRS assumption is
  calculated. For model specification, take a look at Cooper et al. (2007).
  }
\value{
A data frame with J1*(M+6), which has optimal M input factors,
minimized cost when overally efficient, minimized cost when
technically-efficient, revealed cost, overall efficiency, allocative
efficiency, and technical efficiency. 
}

\references{
  Cooper, W., Seiford, L. and Tone, K. (2007). \emph{Data envelopment
  analysis: a comprehensive text with models, applications, references
  and DEA-solver software} (2nd ed.). Springer Verlag, New York.

Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
  Envelopment Analysis}. Press (in Korean).
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{revenue.dea}}}

\examples{

dat.io <- data.frame(y = c(1, 1, 1, 1, 1, 1, 1),
                          x1 = c(2, 3, 5, 9, 6, 3, 8),
                          x2 = c(8, 6, 3, 2, 7, 9, 4))
dat.wm<- c(w1 = 1, w2 = 2)      ## market prices of input factors
(re <- cost.dea(base = dat.io, noutput = 1, input.price = dat.wm))

}

\keyword{Data Envelopment Analysis}