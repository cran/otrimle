\name{otrimle}

\alias{otrimle}

\alias{print.otrimle}

\title{Optimally Tuned Robust Improper Maximum Likelihood Clustering}

\description{
  \code{otrimle} searches for \code{G} approximately Gaussian-shaped
  clusters with/without noise/outliers. The method's tuning controlling  the noise
  level is adaptively chosen based on the data.
}


\usage{
otrimle(data, G, initial=NULL, logicd=NULL,
        npr.max=0.5, erc=50, det.min=.Machine$double.eps, beta=0,
        opt.selector=FALSE, cmstep=TRUE,
        iter.max=100, tol=1e-4, em.iter.max=500, em.tol=1e-6, monitor=1)
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
    assignment  with \code{0} denoting noise/outliers. If \code{NULL} (default)
    initialization is performed using \code{\link{InitClust}}.
  }
  \item{logicd}{
    A vector \code{c(log(icd.min), log(icd.max))}, where \code{icd.min}
    and \code{icd.max}  are the lower and the upper bound of the improper
    constant density (icd) with \code{.Machine$double.xmin <= icd.min <
      icd.max < Inf}.
    If \code{logicd=NULL} (default) the bounds are
    computed based on the data. A pure Gaussian Mixture Model fit
    (obtained when \code{log(icd)=-Inf}) is always included in the search path.
  }
  \item{npr.max}{
    A number in \code{(0,1)} specifying the maximum proportion of
    noise/outliers. This defines the \emph{noise proportion
      constraint}.
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
\item{beta}{
  A number \code{>=0} that specifies the penalty parameter for the noise
  level. When \code{beta=0} (default) the penalty term is inactive.
}
\item{opt.selector}{
  A logical value. When set to \code{TRUE} solutions on the
  border of the parameter space are only explored if an interior
  solution is not available (see \emph{Details}).
}
\item{cmstep}{
  A logical value. When set to \code{TRUE} the \emph{eigenratio constraint} is
  enforced at each M-step of the underlying EM-algorithm (see
  \emph{Details}).
}
\item{iter.max}{
  An integer value specifying the maximum number of iterations for the
  OTRIMLE criterion optimization.
}
\item{tol}{
  Stopping criterion for the OTRIMLE criterion optimization. The
  optimization stops when the difference between two
  successive values of \code{log(icd)} is smaller than \code{tol}.
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
  The \code{otrimle} function allows to
  approximate the OTRIMLE solution with two different versions of the underlying EM-type
  algorithm.

  \describe{
    \item{\emph{ECM-algorithm}: \code{opt.selector=FALSE,
	cmstep=TRUE}}{
      The OTRIMLE search computes the RIMLE (see
      \code{\link{rimle}}) based on the \emph{ECM-algorithm} proposed in Coretto  and
      Hennig (2016). In this case both the \emph{eigenratio constraint}, and the
      \emph{noise proportion constraint} are enforced in each conditional
      M-step of the algorithm.
    }

    \item{\emph{Approximate EM-algorithm}: \code{opt.selector=TRUE,
	cmstep=FALSE}}{
      This corresponds to the algorithm
      proposed in Coretto and Hennig (2015). In this case
      covariance matrices are regularized in each step based on
      \code{det.min}, and the \emph{eigenratio constraint} is applied
	only at the end of
      the EM iteration.
    }
  }

  The \emph{ECM-algorithm} is the default choice. The \emph{Approximate
    EM-algorithm} is often slower than the \emph{ECM-algorithm} by a factor of two. Furthermore the \emph{Approximate
    EM-algorithm}  is
  more prone to lead to problems indicated by \code{code=0} (see
  \emph{Value}-section
  below) because of  numerical degeneracies connected to a low value of
  \code{min.det}. \cr

  There may be datasets for which the function does not provide a
  solution based on default arguments. This corresponds to
  \code{code=0} and \code{flag=1} or 2 or 3 in the output (see
  \emph{Value}-section below).  This usually happens when some (or all) of the
  following circumstances occur: (i) \code{erc} is too
  large; (ii) npr.max is too large; (iii) choice of the initial
  partition. Regarding (i) and (ii) it is not possible to give numeric
  references because whether these numbers are too large/small
  strongly depends on the sample size and the dimensionality of the
  data. References given below explain the relationship between
  these quantities. \cr

  It is suggested to try the following whenever a \code{code=0}
  non-solution occurs. Set \code{logicd} wide enough
  (e.g. \code{logicd=c(-500,-5)}), choose \code{erc=1}, and a low choice
  of \code{npr.max} (e.g. \code{npr.max=0.02}).  Monitor the solution
  with the criterion profiling plot
  (\code{\link{plot.otrimle}}). According to the criterion profiling
  plot change \code{logicd}, and increase \code{erc} and \code{npr.max}
  up to the point when a "clear" minimum in the criterion profiling plot
  is obtained. If this strategy does not work it is suggested to
  experiment with a different initial partitions (see \code{initial} above).
}



\value{
  An S3 object of class \code{'otrimle'} providing the optimal (according to
  the OTRIMLE criterion) clustering. Output components are as follows:

  \item{code}{
    An integer indicator for the convergence.
    \code{code=0} if no solution is found (see \emph{Details});
    \code{code=1} if at the optimal icd value the corresponding EM-algorithm did not
    converge within \code{em.iter.max};
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
    Number of iterations performed in the underlying  EM-algorithm at the
    optimal \code{icd}.
  }
  \item{logicd}{
    Resulting value of the optimal \code{log(icd)}.
  }
  \item{iloglik}{
    Resulting value of the improper likelihood.
  }
  \item{criterion}{
    Resulting value of the OTRIMLE criterion.
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
  \item{optimization}{
    A data.frame with  the OTRIMLE optimization profiling. For each
    value of \code{log(icd)} explored by the algorithm the data.frame
    stores \code{criterion, iloglik, code} and  \code{flag}.
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
  \code{\link{plot.otrimle}},
  \code{\link{InitClust}},
  \code{\link{rimle}},
}







\examples{
## Load  Swiss banknotes data
data(banknote)
x <- banknote[,-1]

## Perform otrimle clustering with default arguments
set.seed(1)
a <- otrimle(data=x, G=2)
print(a)

## Plot clustering
plot(a, data=x, what="clustering")

## Plot OTRIMLE criterion profiling
plot(a, what="criterion")

## Plot Improper log-likelihood profiling
plot(a, what="iloglik")

## P\code{-}P plot of the clusterwise empirical weighted squared Mahalanobis
## distances against the target distribution pchisq(, df=ncol(data))
plot(a, what="fit")
plot(a, what="fit", cluster=1)

}

