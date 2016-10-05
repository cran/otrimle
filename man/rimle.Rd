
\name{rimle}

\alias{rimle}

\alias{print.rimle}

\title{Robust Improper Maximum Likelihood Clustering}

\description{
  \code{rimle} searches for \code{G} approximately Gaussian-shaped
  clusters with/without noise/outliers. The method's tuning controlling
  the noise
  level is fixed and is to be provided by the user or will be guessed by
  the function in
  a rather quick and dirty way (\code{\link{otrimle}} performs a more
  sophisticated data-driven choice).
}


\usage{
rimle(data, G, initial=NULL, logicd=NULL,
      npr.max=0.5, erc=50, det.min=.Machine$double.eps, cmstep=TRUE,
      em.iter.max=500, em.tol=1e-6, monitor=1)
}


\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Rows correspond
    to observations and columns correspond to variables. Categorical
    variables and \code{NA} values are not allowed.
  }
  \item{G}{
    An integer specifying the number of clusters.
  }
  \item{initial}{An integer vector specifying the initial cluster
    assignment  with \code{0=noise/outliers}. If \code{NULL} (default),
    initialization is performed using \code{\link{InitClust}}.
  }
  \item{logicd}{
    A number \code{log(icd)}, where \code{0 <= icd < Inf}
    is the value of the improper constant density (\code{icd}). This is
    the RIMLE's tuning for controlling the size of the noise.
    If \code{logicd=NULL} (default), an \code{icd} value is guessed
    based on the data. A pure Gaussian Mixture Model fit is obtained
    with \code{log(icd)=-Inf}.
  }
  \item{npr.max}{
    A number in \code{(0,1)} specifying the maximum proportion of
    noise/outliers. This defines the \emph{noise proportion constraint}.
  }
  \item{erc}{
    A number \code{>=1} specifying the maximum allowed ratio between
    within-cluster covariance matrix eigenvalues. This defines the
    \emph{eigenratio constraint}. \code{erc=1} enforces spherical clusters
    with equal covariance matrices. A large \code{erc} allows for large
    between-cluster covariance discrepancies. In order to facilitate the
    setting of \code{erc}, it is suggested to scale the columns of
    \code{data} (see \code{\link[base]{scale}}) whenever measurement units
    of the different variables are grossly incompatible.
  }
  \item{det.min}{
    Lower bound for the minimum determinant of covariance matrices. This
    is only active if \code{cmstep=FALSE} (see \emph{Details}).
  }
  \item{cmstep}{
    A logical value. When set to \code{TRUE} the \emph{eigenratio constraint} is
    enforced at each M-step of the underlying EM-algorithm (see
    \emph{Details}).
  }
  \item{em.iter.max}{
    An integer value specifying the maximum number of iterations allowed
    in the underlying EM-algorithm.
  }
  \item{em.tol}{
    Stopping criterion for the the underlying EM-algorithm. An EM iteration
    stops if two successive improper log-likelihood values are within
    \code{em.tol}.
  }
  \item{monitor}{Set the verbosity level of tracing messages.  Possible values
    are \code{monitor=0} (no messages), \code{monitor=1} and
    \code{monitor=2} for increased verbosity.
  }
}


\details{
  The \code{rimle} function allows to approximate the RIMLE solution
  with two  different versions of the underlying EM-type algorithm.

  \describe{
    \item{\emph{ECM-algorithm}: \code{cmstep=TRUE}}{
      The RIMLE solution is obtained based on the \emph{ECM-algorithm}
      proposed in Coretto and Hennig (2016). In this case both the
      \emph{eigenratio constraint} and the \emph{noise proportion
      constraint} are enforced in each conditional M-step of the
      algorithm.
    }

    \item{\emph{Approximate EM-algorithm}: \code{cmstep=FALSE}}{
      This corresponds to the algorithm proposed in Coretto and Hennig
      (2015).  In this case covariance matrices are regularized in each
      step based on \code{det.min}, and the \emph{eigenratio constraint}
      is applied only at the end of the EM iteration.}
  }

  The \emph{ECM-algorithm} is the default choice. The \emph{Approximate
  EM-algorithm} is often slower than the \emph{ECM-algorithm} by a
  factor of two. Furthermore the \emph{Approximate EM-algorithm} is
  more prone to lead to problems indicated by \code{code=0} (see
  \emph{Value}-section below) because of numerical degeneracies
  connected to a low value of \code{min.det}. \cr

  There may be datasets for which the function does not provide a
  solution based on default arguments. This corresponds to
  \code{code=0} and \code{flag=1} or 2 or 3 in the output (see
  \emph{Value}-section below).  This usually happens when some (or all) of the
  following circumstances occur: (i) \code{log(icd)} is too
  large; (ii) \code{erc} is too large; (iii) npr.max is too large;
  (iv) choice of the initial partition. In these cases it is suggested
  to find a suitable interval of \code{icd} values by using the
  \code{\link{otrimle}} function. The \emph{Details} section of
  \code{\link{otrimle}} suggests several actions to take
  whenever a \code{code=0} non-solution occurs.
}



\value{
  An S3 object of class \code{'rimle'}. Output components are as follows:

  \item{code}{
    An integer indicator for the convergence.
    \code{code=0} if no solution is found (see \emph{Details});
    \code{code=1} if the  EM-algorithm did not converge within \code{em.iter.max};
    \code{code=2} convergence is fully achieved.
  }
  \item{flag}{
    A character string containing one or more flags related to
    the EM iteration at the optimal icd.
    \code{flag=1}  if it was not possible to prevent the numerical
    degeneracy of improper posterior probabilities (\code{tau} value
    below).
    \code{flag=2} if enforcement of the \emph{noise proportion constraint}
    failed for numerical reasons.
    \code{flag=3} if enforcement of the \emph{eigenratio constraint} failed for
    numerical reasons.
    \code{flag=4} if the  \emph{noise proportion constraint} has been
    successfully applied at least once.
    \code{flag=5} if the  \emph{eigenratio constraint} has been
    successfully applied at least once.
  }
  \item{iter}{
    Number of iterations performed in the underlying  EM-algorithm.
  }
  \item{logicd}{
    Value of the \code{log(icd)}.
  }
  \item{iloglik}{
    Value of the improper likelihood.
  }
  \item{criterion}{
    Value of the OTRIMLE criterion.
  }
  \item{npr}{
    Estimated expected noise proportion.
  }
  \item{cpr}{
    Vector of estimated expected cluster proportions (notice that \code{sum(cpr)=1-npr}).
  }
  \item{mean}{
    A matrix of dimension \code{ncol(data) x G} containing the
    mean parameters  of each cluster (column-wise).
  }
  \item{cov}{
    An array of size \code{ncol(data) x ncol(data) x G}
    containing the covariance matrices of each cluster.
  }
  \item{tau}{A matrix of dimension \code{nrow(data) x {1+G}} where
    \code{tau[i, 1+j]} is the estimated (improper) posterior probability that
    the \emph{i}th observation belongs to the \emph{j}th cluster.
    \code{tau[i,1]} is the estimated (improper) posterior probability that
    \emph{i}th observation belongs to the noise component.
  }
  \item{smd}{
    A matrix of dimension  \code{nrow(data) x G } where \code{smd[i,j]}
    is the squared Mahalanobis distance of \code{data[i,]} from
    \code{mean[,j]} according to \code{cov[,,j]}.
  }
  \item{cluster}{
    A vector of integers  denoting cluster assignments for each
    observation. It's \code{0} for observations assigned to  noise/outliers.
  }
  \item{size}{
    A vector of integers with sizes (counts) of each cluster.
  }
}







\section{References}{
  Coretto, P.  and C. Hennig (2015).
  Robust improper maximum likelihood: tuning, computation, and a comparison with other methods for robust Gaussian clustering.
  To appear on the \emph{Journal of the American Statistical Association}.
  arXiv preprint at \href{http://arxiv.org/abs/1406.0808}{arXiv:1406.0808}
  with (\href{http://arxiv.org/src/1406.0808v4/anc/supplement.pdf}{supplement}).


  Coretto, P.  and C. Hennig (2016).
  Consistency, breakdown robustness, and algorithms for robust improper maximum likelihood clustering.
  arXiv preprint at \href{http://arxiv.org/abs/1309.6895}{arXiv:1309.6895}.
}







\seealso{
  \code{\link{plot.rimle}},
  \code{\link{InitClust}},
  \code{\link{otrimle}},
}







\examples{
## Load  Swiss banknotes data
data(banknote)
x <- banknote[,-1]

## -----------------------------------------------------------------------------
## EXAMPLE 1:
## Perform RIMLE with default inputs
## -----------------------------------------------------------------------------
set.seed(1)
a <- rimle(data=x, G=2)
print(a)

## Plot clustering
plot(a, data=x, what="clustering")

## P-P plot of the clusterwise empirical weighted squared Mahalanobis
## distances against the target distribution pchisq(, df=ncol(data))
plot(a, what="fit")
plot(a, what="fit", cluster=1)



## -----------------------------------------------------------------------------
## EXAMPLE 2:
## Compare solutions for different choices of logicd
## -----------------------------------------------------------------------------
set.seed(1)

## Case 1: noiseless solution, that is fit a pure Gaussian Mixture Model
b1 <- rimle(data=x, G=2, logicd=-Inf)
plot(b1, data=x, what="clustering")
plot(b1, what="fit")

## Case 2: low noise level
b2 <- rimle(data=x, G=2, logicd=-100)
plot(b2, data=x, what="clustering")
plot(b2, what="fit")

## Case 3: medium noise level
b3 <- rimle(data=x, G=2, logicd=-10)
plot(b3, data=x, what="clustering")
plot(b3, what="fit")

## Case 3: large noise level
b3 <- rimle(data=x, G=2, logicd=5)
plot(b3, data=x, what="clustering")
plot(b3, what="fit")

}

