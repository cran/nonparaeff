\name{revenue.dea}
\alias{revenue.dea}
\title{Linear Programming for Revenue Maximization}
\description{Solve the Revenue Maximization Probem with Given Output Prices}
\usage{
revenue.dea(base = NULL, frontier = NULL, noutput = 1, output.price = NULL)
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
  \item{output.price}{A vector for market prices of outputs.}
}
\details{
  The revenue maximization problem under the CRS assumption is
  calculated. See Cooper et al. (2007).
  }
\value{
A data frame with J1*(N+6), which has optimal N output factors,
maximized revenue when overally efficient, maximized revenue when
technically-efficient, revealed revenue, overall efficiency, allocative
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
\seealso{\code{\link{cost.dea}}}

\examples{
tab8.3 <- data.frame(y1 = c(1, 3, 6, 6, 3, 9),
                          y2 = c(6, 6, 3, 5, 4, 1),
                          x = c(1, 1, 1, 1, 1, 1))
tab8.3.ps.f <- c(p1 = 2, p2 = 2)
(ex8.3 <- revenue.dea(base = tab8.3,
                    noutput = 2, output.price = tab8.3.ps.f))
}

\keyword{Data Envelopment Analysis}