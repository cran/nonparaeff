\name{sbm.tone}
\alias{sbm.tone}
\title{Linear Programming for the Slacks-based Model under the CRS}
\description{Solve Slacks-based Model under the CRS (Tone, 2001)}
\usage{
sbm.tone(base= NULL, frontier = NULL, noutput = 1)
}
\arguments{
  \item{base}{A data set for DMUs to be evaluated. A data frame with
    J1*(M+N) dimention, where J1 is the number of DMUs, M for the number
    of inputs, and N for the number of outputs.}
  \item{frontier}{A data set for DMUs to construct a production
    possibility set (PPS). A data frame with
    J2*(M+N) dimention, where J2 is the number of DMUs, M for the number
    of inputs, and N for the number of outputs.}
  \item{noutput}{The number of outputs (N).}
}
\details{
  The SBM under the CRS assumption is calculated. See Tone (2001).
}

\value{
  A data frame with (1+J1+M+N), which shows efficiency scores, optimal
  weights, and optiaml input and output slacks.
}

\references{
  Cooper, W., Seiford, L. and Tone, K. (2007). \emph{Data envelopment
  analysis: a comprehensive text with models, applications, references
  and DEA-solver software} (2nd ed.). Springer Verlag, New York.

Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
    Envelopment Analysis}. Press (in Korean).

  Tone, K. (2001). A slacks-based measure of efficiency in data
  envelopment analysis. \emph{European Journal of Operational Research},
  130(3):498-509.
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{sbm.vrs}}}

\examples{
tab7.6.dat <- data.frame(y = c(1, 1, 1, 1, 1, 1),
                              x1 = c(1, 3, 6, 2, 5, 9),
                              x2 = c(4, 1, 1, 8, 5, 2))
(re <- sbm.tone(tab7.6.dat, noutput = 1))
}

\keyword{Data Envelopment Analysis}