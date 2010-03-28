\name{additive}
\alias{additive}
\title{Linear Programming for the Additive Model}
\description{Solve the Additive Model under the VRS assumption}
\usage{
additive(base = NULL, frontier = NULL, noutput = 1)
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
}
\details{
  The additive model under the VRS assumption is as follows:

  \eqn{
    \theta^{k*}_{ADD} = \max_{\lambda, s^+, s^-} \left( \sum_{m=1}^M
    s^-_m + \sum_{n=1}^N s^+_n \right)}{}
  
  \eqn{s.t.~~ x^k_m = \sum_{j=1}^J x_m^j \lambda^j + s^-_n ~ (m = 1, 2,
    \cdots, M);}{}
  
  \eqn{~~~~~~ y^k_n = \sum_{j=1}^J y_n^j \lambda^j + s^-_n ~ (n = 1, 2,
    \cdots, N);}{}
  
  \eqn{~~~~~~\sum_{j=1}^J \lambda^j = 1;}{}
  \eqn{~~~~~~\lambda^j \geq 0~ (j = 1, 2, \cdots, J);~ s^-_m \geq 0~ (m = 1, 2, \cdots, M);~ s^+_n \geq 0~ (n = 1, 2, \cdots, N).
  }{
    theta^k* = max sum(s^-_m) + sum(s^+_n)

    s.t.   x^k_m = sum(x_m^j) l^j + s^-_m  (m = 1, 2, ..., M);
    
    y^k_n = sum(y_n^j) l^j + s^+_n  (n = 1, 2, ..., N);
    
    sum(l^j) = 1;
    
    l^j >= 0, s^-_m >= 0, s^+_n >= 0.
  }
}

\value{
A data frame with J1*(J1+M+N), which has efficiency scores, optimal
weightes and optimal slacks. Take a look at the example below.
}

\references{
  Cooper, W., Seiford, L. and Tone, K. (2007). \emph{Data envelopment
  analysis: a comprehensive text with models, applications, references
  and DEA-solver software} (2nd ed.). Springer Verlag, New York.

Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
  Envelopment Analysis}. Press (in Korean).
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{sbm.tone}}, \code{\link{sbm.vrs}}}

\examples{

## Simple Example
my.dat <- data.frame(y = c(1, 2, 4, 6, 7, 9, 9),
                     x = c(3, 2, 6, 4, 8, 8, 10))
(re <- additive(my.dat, noutput = 1))

## Property of the Additive Model
dat1 <- data.frame(y = c(1, 1, 1, 1, 1, 1),
                        x1 = c(2, 3, 6, 3, 6, 6),
                        x2 = c(5, 3, 1, 8, 4, 2))
dat2 <- dat1 
dat2$x1 <- dat2$x1 * 10 
dat3 <- dat1 
dat3$x1 <- dat3$x1 + 10 
(re1 <- additive(dat1, noutput = 1))
(re2 <- additive(dat2, noutput = 1))
(re3 <- additive(dat3, noutput = 1))

}

\keyword{Data Envelopment Analysis}