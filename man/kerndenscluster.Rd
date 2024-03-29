\name{kerndenscluster}
\alias{kerndenscluster}

\title{Aggregated distance to elliptical unimodal density over clusters}

\description{
  This calls \code{\link{kerndensp}} for computing and aggregating
  density- and
  principal components-based distances between
  multivariate data and a unimodal
  elliptical distribution about the data mean for all clusters in a
  mixture-based clustering as generated by \code{\link{otrimle}} or
  \code{\link{rimle}}. For use in \code{\link{otrimleg}}.
} 

\usage{
  kerndenscluster(x,fit,maxq=qnorm(0.9995),kernn=100)
}


\arguments{
  \item{x}{something that can be coerced into a matrix. Dataset.}
  \item{fit}{output object of \code{\link{otrimle}} or \code{\link{rimle}}.}
  \item{maxq}{positive numeric. One-dimensional densities are evaluated
    between \code{mean(x)-maxq} and \code{mean(x)+maxq}.}
  \item{kernn}{integer. Number of points at which the one-dimensional
    density is evaluated, input parameter \code{n} of
    \code{\link{density}}. This should be even.}
}

\details{
  See Hennig and Coretto (2021), Sec. 4.2. \code{kerndenscluster} calls
  \code{\link{kerndensp}} for all clusters and aggregates the resulting
  measures as root sum of squares. 
}
  
\value{
  A list with components \code{ddpi, ddpm, measure}.

  \item{ddpi}{list of outputs of \code{\link{kerndensp}} for all
    clusters.}
  \item{ddpm}{vector of \code{measure}-components of \code{ddpi}.}
  \item{measure}{Final aggregation result.}
}

\references{
  Hennig, C. and P.Coretto (2021). An adequacy approach for deciding the
  number of clusters for OTRIMLE robust Gaussian mixture based
  clustering. To appear in \emph{Australian and New Zealand Journal of
  Statistics}, \url{https://arxiv.org/abs/2009.00921}.
  }

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}

   
\seealso{
  \code{\link{kerndensp}}, \code{\link{kerndensmeasure}},
  \code{\link{otrimle}},  \code{\link{rimle}}
}

\examples{
   data(banknote)
   selectdata <- c(1:30,101:110,117:136,160:161)
   set.seed(555566)
   x <- banknote[selectdata,5:7]
   ox <- otrimle(x, G=2, ncores=1)
   kerndenscluster(x,ox)$measure
}

\keyword{htest}
\keyword{cluster}


