\name{dea}
\alias{dea}
\title{Linear Programming for the Data Envelopment Analysis}
\description{Solve input(output)-oriented DEA under the CRS (VRS)}
\usage{
dea(base = NULL, frontier = NULL, noutput = 1, orientation=1, rts = 1, onlytheta = FALSE)
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
  \item{onlytheta}{Logical. If onlytheta is TRUE, then only efficiency
    scores are obtained. If it is FALSE, then optimal lambda's and
    slacks are also obtained.}
}
\details{
  The input (output) -oriented DEA under the CRS (VRS) assumption are
  calcuated. For model specification, take a look at Cooper et al. (2007).
}

\value{
  If onlytheta is TRUE, then a (J1*1) data.frame is obtained. If
  onlytheta if FALSE, then a data frame with a J1*(J1+M+N) dimension is
  obtained, in which optimal weights, input slacks and output slacks are
  presented.
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
(re <- dea(base = tab3.1.dat, noutput = 1, orientation = 1, rts = 1,
onlytheta = FALSE))

## input-oriented DEA under the CRS assumption (2 inputs and 1 output)
tab3.3.dat <- data.frame(y = c(1, 1, 1, 1, 1, 1),
                              x1 = c(1, 3, 6, 2, 5, 9),
                              x2 = c(4, 1, 1, 8, 5, 2))
re <- dea(base=tab3.3.dat, noutput = 1, orientation = 1, rts = 1)
## finding references points
(ref <- data.frame(y = c(tab3.3.dat$y + re$slack.y1),
x1 = c(tab3.3.dat$x1 * re$eff - re$slack.x1),
x2 = c(tab3.3.dat$x2 * re$eff - re$slack.x2)))


## output-oriented DEA under the CRS assumption (1 input and 2 outputs)
tab5.1.dat <- data.frame(y1 = c(4, 8, 8, 4, 3, 1),
                         y2 = c(9, 6, 4, 3, 5, 6),
                         x = c(1, 1, 1, 1, 1, 1))
(re <- dea(tab5.1.dat, noutput = 2, orientation = 2, rts = 1))

## input-oriented DEA under the VRS assumption (1 input and 1 output)
tab6.1.dat <- data.frame(y = c(1, 2, 4, 6, 7, 9, 9),
                              x = c(3, 2, 6, 4, 8, 8, 10))
(re <- dea(tab6.1.dat, noutput = 1, orientation = 1, rts = 2))

## output-oriented DEA under the VRS assumtion (1 input and 1 output)
(re <- dea(tab6.1.dat, noutput = 1, orientation = 2, rts = 2))

## scale efficiency
re.crs <-
    dea(tab6.1.dat, noutput = 1, orientation = 1, rts = 1,onlytheta = TRUE)
re.vrs<-
    dea(tab6.1.dat, noutput = 1, orientation = 1, rts = 2,
         onlytheta = TRUE)
scale.eff <- re.crs/re.vrs

## finding DRS, IRS, CRS
dat6.1 <- data.frame(y = c(1, 2, 4, 6, 7, 9, 9),
     x = c(3, 2, 6, 4, 8, 8, 10))
re <- dea(dat6.1, noutput = 1, rts = 1)
lambdas <- re[, 2:8]
apply(lambdas, 1, sum) 
}

\keyword{Data Envelopment Analysis}