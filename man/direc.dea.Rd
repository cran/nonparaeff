\name{direc.dea}
\alias{direc.dea}
\title{Linear Programming for the Directional Distance Function with
  Undesirable Outputs}
\description{Solve the DDF with undesirable outputs. The directional
  vecor is (y's, b's).}
\usage{
  direc.dea(base = NULL, frontier = NULL, ngood = 1, nbad = 1)
}
\arguments{
  \item{base}{A data set for DMUs to be evaluated. A data frame with
    J1*(M+P+Q) dimention, where J1 is the number of DMUs, M for the number
    of inputs, P for the number of good outputs, and Q for the
    undesirable outputs.}
  \item{frontier}{A data set for DMUs to be used in constructing a
    production possibility set (PPS). A data frame with J2*(M+P+Q)
    dimention, where J2 is the number of DMUs, M for the number of
    inputs, P for the number of good outputs, and Q for the
    undesirable outputs}
  \item{ngood}{The number of good outputs (P).}
  \item{nbad}{The number of bad outputs (Q).}
}
\details{
  The DDF with undesirable outputs under the CRS assumption is
  calculated. For model specification, take a look at Chung et al. (1997).
}

\value{
A J1 vector of which is inefficiency score.
}

\references{
  Chung, Y. Fare, R. and Grosskopf, S. (1997). Productivity and
  undesirable outputs: A directional distance function
  approach. \emph{Journal of Environmental Management} 51(3):229-240.

  Cooper, W., Seiford, L. and Tone, K. (2007). \emph{Data envelopment
  analysis: a comprehensive text with models, applications, references
  and DEA-solver software} (2nd ed.). Springer Verlag, New York.

  Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
  Envelopment Analysis}. Press (in Korean).
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{ddf}}}

\examples{
## Simple Example of one input, one good output, and one bad output.
my.dat <- data.frame(yg = c(2, 5, 7, 8, 3, 4, 6),
                     yb = c(1, 2, 4, 7, 4, 5, 6),
                     x = c(1, 1, 1, 1, 1, 1, 1))
direc.dea(my.dat, ngood = 1, nbad = 1)
}

\keyword{Data Envelopment Analysis} 