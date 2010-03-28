\name{int.dea}
\alias{int.dea}
\title{Linear Programming for the Data Envelopment Analysis with
  Integer-valued Inputs.}
\description{Solve input-oriented DEA under the CRS}
\usage{
int.dea(base = NULL, frontier = NULL, noutput = 1, intinput = 1,
orientation=1, epsilon = 1e-06)
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
  \item{intinput}{The number of integer inputs.}
  \item{orientation}{Orientation of measurement. 1 for the input-oriented
    measure, and 2 for the output-oriented measure. This argument is ignored.}
  \item{epsilon}{Non-Armechidean number. Use the default value.}
}
\details{
  The input-oriented IDEA under the CRS assumption is calcualted. See
  Kuosmanen and Matin (2009). 
}

\value{
  A data frame of J1*(1+J1+N+M+Q+Q), which shows efficiency scores,
  optimal weightes, optimal slacks for outputs and inputs, optiaml
  slacks for integer-valued inputs, and optimal integer inputs.
}

\references{
  Kuomanen, T. and Matin, R. (2009). Theory of integer-valued data
  envelopment analysis. \emph{European Journal of Operational Research}
  192(2):658-667

  Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
    Envelopment Analysis}. Press (in Korean).
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{dea}}}

\examples{
int.dat <- data.frame(y = c(1, 1, 1, 1, 1),
                      x1 = c(2, 7, 3, 7, 9),
                      x2 = c(4, 1, 4, 2, 4))
int.dea(int.dat, noutput = 1, intinput = 1)
}

\keyword{Data Envelopment Analysis}