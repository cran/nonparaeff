\name{ar.dual.dea}
\alias{ar.dual.dea}
\title{Assurance Region Data Envelopment Aanlysis (AR-DEA)}
\description{Solve the AR-DEA}
\usage{
ar.dual.dea(base = NULL, frontier = NULL,
           noutput = 1, orientation=1, rts = 1, ar.l = NULL,
           ar.r = NULL, ar.dir = NULL, dual = FALSE)
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
  \item{ar.l}{A data frame for the assurance region of which is the
    left-hand.}
  \item{ar.r}{A vector for the assurance region of which is the
    right-hand.}
  \item{ar.dir}{A vector for the assurance region of which is the
    direction.}
  \item{dual}{Logical.}
}
\details{
  The AR model under the CRS assumption is calculated. For model
  specification, take a look at Cooper et al. (2007).
}

\value{
A data frame with J1*(M+N), which has efficiency scores, optimal
virtual prices. Take a look at the example below.
}

\references{
  Cooper, W., Seiford, L. and Tone, K. (2007). \emph{Data envelopment
  analysis: a comprehensive text with models, applications, references
  and DEA-solver software} (2nd ed.). Springer Verlag, New York.

  Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
    Envelopment Analysis}. Press (in Korean).
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{dea}}, \code{\link{dual.dea}}}

\examples{

## AR constraint of 0.25 <= v2/v1 <= 1.
library(Hmisc)
library(lpSolve)
ar.dat <- data.frame(y = c(1, 1, 1, 1, 1, 1),
                     x1 = c(2, 3, 6, 3, 6, 6),
                     x2 = c(5, 3, 1, 8, 4, 2))
(re <-
ar.dual.dea(ar.dat, noutput = 1, orientation = 1, rts = 1, ar.l =
matrix(c(0, 0, 0.25, -1, -1, 1), nrow = 2, ncol = 3), ar.r = c(0, 0),
ar.dir = c("<=", "<=")))

}

\keyword{Data Envelopment Analysis}