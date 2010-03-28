\name{orderm}
\alias{orderm}
\title{Efficiency Measures with the order-m Method.}
\description{Calculate order-m efficiency scores}
\usage{
orderm(base = NULL, frontier = NULL, noutput = 1, orientation=1, M = 25, B = 500)
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
  \item{M}{The number of elements in each of the bootstrapped samples.}
  \item{B}{The number of bootstap replicates}
}
\details{
  See Simar (2003).
}

\value{
  A data frame with J1*1 dimention, which shows efficiency scores. 
}

\references{
  Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
    Envelopment Analysis}. Press (in Korean).

  Simar, L. (2003). Detecting outliers in frontier models: A simple
  approach. \emph{Journal of Productivity Analysis}, 20(3):391-424.
  }

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{fdh}}}

\examples{
x <- abs(runif(200, min = 0.1, max = 4)) 
y <- 3*x*abs(rnorm(200))
dat.orderm <- data.frame(y = y, x = x)  
dat.orderm.out <- rbind(dat.orderm, c(4, 0.1)) ## add one outlier.
(eff <- orderm(dat.orderm.out, noutput = 1, M = 25, B = 20))
}

\keyword{Data Envelopment Analysis}