\name{dual.dea}
\alias{dual.dea}
\title{Linear Programming for the Dual Data Envelopment Analysis}
\description{Solve the Dual DEA}
\usage{
dual.dea(base = NULL, frontier = NULL, noutput = 1, orientation=1, rts = 1)
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
    \item{orientation}{Orientation of measurement. 1 for the input-oriented
    measure, and 2 for the output-oriented measure.}
  \item{rts}{Returns to scale. 1 for the CRS assumption, and 2 for the
    VRS assumption.}
}
\details{
  The input-oriented dual DEA under the CRS assumption is
  calculated. For model specification, take a look at Cooper et
  al. (2007).
}
\value{
A data frame with J1*(1+M+N) dimension, of which has efficiency scores, optimal
virtual prices for inputs and outputs.
}

\references{
  Cooper, W., Seiford, L. and Tone, K. (2007). \emph{Data envelopment
  analysis: a comprehensive text with models, applications, references
  and DEA-solver software} (2nd ed.). Springer Verlag, New York.

  Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
  Envelopment Analysis}. Press (in Korean).
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{dea}}, \code{\link{ar.dual.dea}}}

\examples{
## An output-oriented primal problem with 1 input and 2 outputs
tab5.1.dat <- data.frame(y1 = c(4, 8, 8, 4, 3, 1),
                               y2 = c(9, 6, 4, 3, 5, 6),
                               x = c(1, 1, 1, 1, 1, 1))
(re <- dea(tab5.1.dat, noutput = 2, orientation = 2, rts = 1))

## An output-oriented dual problem with 1 input and 2 outputs
re <- dual.dea(tab5.1.dat, noutput = 2, orientation = 2, rts = 1)
}

\keyword{Data Envelopment Analysis}