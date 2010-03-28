\name{fdh}
\alias{fdh}
\title{Linear Programming for the Free Disposable Hull}
\description{Solve input(output)-oriented FDH}
\usage{
fdh(base = NULL, frontier = NULL, noutput = 1, orientation=1)
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
}
\details{
  The input (output) -oriented FDH is calculated. 
}

\value{
  A data frame of J1*1 dimention which shows efficiency scores.
}

\references{
    Cooper, W., Seiford, L. and Tone, K. (2007). \emph{Data envelopment
  analysis: a comprehensive text with models, applications, references
  and DEA-solver software} (2nd ed.). Springer Verlag, New York.

  Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
    Envelopment Analysis}. Press (in Korean).
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{dea}}, \code{\link{orderm}}}

\examples{
## input-oriented FDH with 1 input and 1 output.
tab7.1.dat <- data.frame(y = c(1, 2, 4, 6, 7, 9, 9),
                              x = c(3, 2, 6, 4, 8, 8, 10))
(re <- fdh(tab7.1.dat, noutput = 1, orientation = 1))

## input-oriented FDH with 2 input and 1 output.
tab7.10.dat <- data.frame(y = c(1, 1, 1, 1, 1, 1),
                               x1 = c(2, 3, 6, 3, 6, 6),
                               x2 = c(5, 3, 1, 8, 4, 2))
(re <- fdh(tab7.10.dat, noutput = 1, orientation = 1))
}

\keyword{Data Envelopment Analysis}