\name{faremalm2}
\alias{faremalm2}
\title{Linear Programming for the Malmquist Productivity Growth Index}
\description{Calculate productivity growth index under the DEA framework.}
\usage{
faremalm2(dat = NULL, noutput = 1, id = "id", year = "year")
}
\arguments{
  \item{dat}{A data frame to be evaluated. The format of this data frame
  is data.frame(id, year, outputs, inputs). This data frame should have a
  balanced panel data form.}
  \item{noutput}{The number of outputs.}
  \item{id}{A column name for the producer index.}
  \item{year}{A column name for the time index.}
}
\details{
  The Malmquist productivity growth index is calculated. For model
  specification, take a look at Fare et al. (1994).
}

\value{A data frame with (
  id: the id index of the original data.
  time: the time index of the original data.
  y's: original outputs
  x's: original inputs
  Dt2t2: \eqn{D^{t+1} (x^{t+1}, y^{t+1})}{D^{t+1} (x^{t+1}, y^{t+1})}
  Dtt2: \eqn{D^{t} (x^{t+1}, y^{t+1})}{D^{t} (x^{t+1}, y^{t+1})}
  Dt2t: \eqn{D^{t+1} (x^t, y^t)}{D^{t+1} (x^t, y^t)}
  ec: efficiency change
  tc: technical change
  pc: productivity change
}

\references{
  Cooper, W., Seiford, L. and Tone, K. (2007). \emph{Data envelopment
  analysis: a comprehensive text with models, applications, references
  and DEA-solver software} (2nd ed.). Springer Verlag, New York.

Fare, R., Grosskopf, S., Norris, M. and Zhang, Z. (1994). Productivity
growth, technical progress and efficiency change in industrialized
countries. \emph{American Economic Review}, 84(1):66-83.

Lee, J. and Oh, D. (forthcoming). \emph{Efficiency Analysis: Data
  Envelopment Analysis}. Press (in Korean).
}

\author{Dong-hyun Oh, \email{oh.donghyun77@gmail.com}}
\seealso{\code{\link{dea}}}

\examples{
malm.dat <- data.frame(id = rep(LETTERS[1:3], 3),
       time = rep(1:3, each = 3),
                       y = c(1, 2, 2, 3, 2, 8, 3, 2, 5),
                       x = c(2, 3, 7, 3, 5, 6, 8, 9, 6))

malm.re1 <- faremalm2(malm.dat, noutput = 1, id = "id", year = "time")


## Malmquist productivity growth index of OECD countries
library(pwt)       ## Use Penn World Table
my.dat <- pwt5.6    
head(my.dat)        
my.oecd.ctry <- c("AUS", "AUT", "BEL", "CAN", "CHE", "DNK", "ESP",
                    "FIN", "FRA", "GBR", "GER", "GRC", "IRL", "ISL",
                    "ITA", "JPN", "KOR", "LUX", "MEX", "NLD", "NOR",
                    "NZL", "PRT", "SWE", "TUR", "USA", "DEU")
my.dat <- my.dat[my.dat$wbcode \%in\% my.oecd.ctry,]
my.dat <- my.dat[my.dat$year \%in\% 1980:1990,]
my.dat$rgdpl <- as.numeric(my.dat$rgdpl) ## GDP per capita
my.dat$pop <- as.numeric(my.dat$pop) ## total population (1000)
my.dat$rgdpwok <- as.numeric(my.dat$rgdpwok) ## GDP per labor
my.dat$kapw <- as.numeric(my.dat$kapw)  ## Capital stock per labor
my.dat$gdp <- my.dat$rgdpl * my.dat$pop ## Total GDP of a country
my.dat$labor <- with(my.dat, gdp/rgdpwok) ## Total labor force
my.dat$capital <- with(my.dat, kapw * labor) ## Toal capital stock
oecd <- my.dat[, c("wbcode", "year", "gdp", "labor", "capital")] 
re.oecd <- faremalm2(dat = oecd, noutput = 1, id = "wbcode", year =
"year")
## productivity growth for each country
pc.c <- tapply(re.oecd$pc, re.oecd$wbcode, geometric.mean)
## a trend of productivity growth of OECD countries
pc.y <- tapply(re.oecd$pc, re.oecd$year, geometric.mean)
## efficiency change for each country
ec.c <- tapply(re.oecd$ec, re.oecd$wbcode, geometric.mean)
## a trend of efficiency change of OECD countries
ec.y <- tapply(re.oecd$ec, re.oecd$year, geometric.mean)
}

\keyword{Data Envelopment Analysis}