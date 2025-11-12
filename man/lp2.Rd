\name{lp2}
\alias{lp2}
\title{Linear Programming with Free Variables}
\description{Solve LP with free variables}
\usage{
lp2(direction = "min", objective.in, const.mat, const.dir, 
    const.rhs, free.var = NULL)
}
\arguments{
  \item{direction}{Character string giving direction of optimization: "min"
    (default) or "max."}
  \item{objective.in}{Numeric vector of coefficients of objective
    function}
  \item{const.mat}{Matrix of numeric constraint coefficients, one row per
    constraint, one column per variable (unless
    transpose.constraints = FALSE; see below).}
  \item{const.dir}{Vector of character strings giving the direction of the
    constraint: each value should be one of "<," "<=," "=," "==,"
    ">," or ">=". (In each pair the two values are identical.)}
  \item{const.rhs}{Vector of numeric values for the right-hand sides of the
    constraints.}
  \item{free.var}{Vector of numeric values for indicating free
    variables. If this argument is NULL, no free variables is included.}
}
\details{
  lp2 extends lpSolve::lp() to incorporate free variables easily.
}

\value{
  An lp object. See 'lp.object' for details.
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}

\seealso{\code{\link[lpSolve:lp]{lp}}}

\examples{
     # Set up problem: maximize
     #   x1 + 9 x2 +   x3 subject to
     #   x1 + 2 x2 + 3 x3  <= 9
     # 3 x1 + 2 x2 + 2 x3 <= 15
     #
     f.obj <- c(1, 9, 3)
     f.con <- matrix (c(1, 2, 3, 3, 2, 2), nrow=2, byrow=TRUE)
     f.dir <- c("<=", "<=")
     f.rhs <- c(9, 15)
     #
     # Now run.
     #
     lp2("max", f.obj, f.con, f.dir, f.rhs)
     lp2("max", f.obj, f.con, f.dir, f.rhs, free.var = c(0, 1, 0))
}
\keyword{Linear Programming}
