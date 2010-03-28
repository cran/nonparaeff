\name{effdea.b.f}
\alias{effdea.b.f}
\title{Linear Programming for the Data Envelopment Analysis}
\description{Solve input(output)-oriented DEA under the CRS (VRS) with
  convexhull. Do not use when the total number of inputs and outputs are
  greater than eight. If used, it may take more than hundreds day to get
  results.}
\usage{
effdea.b.f(base = NULL, frontier = NULL, noutput = 1,
                       orientation=1, rts = 1, convhull = TRUE)
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
\item{convhull}{Logical. If this is TRUE, very efficient calculation of
efficiency score is used. However, when the total number of inputs and
outputs is larger than eight, it is very slow for this option. In cases
when the total number of inputs and outputs is larger than eight, use
FALSE for this argument.}
}
\details{
  This function uses the convhull function in geometry package. After
  finding convex hull of \emph{frontier} by using the convhull
  function. points on the convex hull are used in constructing
  the second production possibility set (PPS). Then efficiency scores in
  \emph{base} are calculated based on the second PPS.
}

\value{A data frame with J1*1 dimension, which shows efficiency scores. 
}

\references{
  Cooper, W., Seiford, L. and Tone, K. (2007). \emph{Data envelopment
    analysis: a comprehensive text with models, applications, references
    and DEA-solver software} (2nd ed.). Springer Verlag, New York.

  Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
    Envelopment Analysis}. Press (in Korean)
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{dual.dea}}}

\examples{
## input-oriented DEA under the CRS assumption (1 input and 1 output)
tab3.1.dat <- data.frame(y = c(1, 2, 4, 6, 7, 9, 9), 
                         x = c(3, 2, 6, 4, 8, 8, 10))
(re <- effdea.b.f(base = tab3.1.dat, noutput = 1, orientation = 1, rts =
1, convhull = TRUE))
}

\keyword{Data Envelopment Analysis}