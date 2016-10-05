\name{InitClust}

\alias{InitClust}

\title{Robust Initilization for Model-based Clustering Methods}

\description{
  Computes the initial cluster assignment based on a combination of
  nearest neighbor based clutter/noise detection,
  and agglomerative hierarchical clustering based on maximum likelihood
  criteria for Gaussian mixture models.
}

\usage{
 InitClust(data, G, cpr.min={ncol(data)+1}/nrow(data),
           K=5, nstart.km=50, modelName="VVV", monitor=FALSE)
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
   \item{cpr.min}{
     The minimum cluster proportion allowed in the initial clustering.
   }
   \item{K}{
     An integer specifying the number of considered nearest neighbors per point
     used for the denoising step (see \emph{Details}).
   }
   \item{nstart.km}{
     An integer specifying the number of  random starts  for the k-means step.
   }
   \item{modelName}{
     A character string indicating the covariance model to be used. Possible models are: \cr
     \code{"E"}:   equal variance  (one-dimensional) \cr
     \code{"V"} :  spherical, variable variance (one-dimensional) \cr
     \code{"EII"}: spherical, equal volume \cr
     \code{"VII"}: spherical, unequal volume \cr
     \code{"EEE"}: ellipsoidal, equal volume, shape, and orientation \cr
     \code{"VVV"}: ellipsoidal, varying volume, shape, and orientation (default).\cr
     See \emph{Details}.
   }
   \item{monitor}{
     A logical value; \code{TRUE} means that tracing messages
     will be produced.
   }

 }


 \details{
   The initialization is described in the supplementary material of
   Coretto and Hennig (2015). Noise/outliers are removed based on nearest neighbor based clutter/noise
   detection  (NNC) of Byers and Raftery (1998). This step is performed
   with  \code{\link[prabclus]{NNclean}}. The input argument \code{K}
   is passed as \code{k} to
   \code{\link[prabclus]{NNclean}}. Based on
   this step a denoised version of \code{data} is obtained. The initial
   clustering is then obtained based on the following steps. Note
   that these steps are reported in the \code{code} element of the output
   list (see \emph{Value}).\cr

   Clustering steps:\cr

   \emph{Step 1}: perform the model-based hierarchical clustering (MBHC)
   proposed in Fraley (1998). This step is performed using
   \code{\link[mclust]{hc}}. The input argument \code{modelName} is passed
   to \code{\link[mclust]{hc}}. See \emph{Details} of
   \code{\link[mclust]{hc}} for more details.

   \emph{Step 2}: if too small clusters (cluster proportions
   \code{<cpr.min}) are found in the previous step,  assign small clusters
   to noise and perform MBHC again on the denoised data.

   \emph{Step 3}: if too small clusters are found in the previous step,
   assign small clusters to noise and perform k-means on the denoised data.

   \emph{Step 4}: if too small clusters are found in the previous step, then a
   completely random partition that satisfies  \code{cpr.min} is returned.
 }


 \value{
   A \code{list} with the following components:

   \item{code}{
     An integer indicating the step at which the initial clustering has been
     found (see \emph{Details}).
   }
   \item{cluster}{
     A vector of integers  denoting cluster assignments for each
     observation. \code{cluster=0} for observations assigned to  noise/outliers.
   }
 }







 \section{References}{
   Fraley, C.  (1998).
   Algorithms for model-based Gaussian hierarchical clustering.
   \emph{SIAM Journal on Scientific Computing} 20:270-281.

   Byers, S. and A. E. Raftery  (1998).
   Nearest-Neighbor Clutter Removal for Estimating Features in Spatial
   Point Processes,
   \emph{Journal of the American Statistical Association}, 93, 577-584.

   Coretto, P.  and C. Hennig (2015).
   Robust improper maximum likelihood: tuning, computation, and a comparison with other methods for robust Gaussian clustering.
   To appear on the \emph{Journal of the American Statistical Association}.
   arXiv preprint at \href{http://arxiv.org/abs/1406.0808}{arXiv:1406.0808}
   with (\href{http://arxiv.org/src/1406.0808v4/anc/supplement.pdf}{supplement}).
 }


 \seealso{
   \link[prabclus]{NNclean},
   \link[mclust]{hc}
 }


 \examples{
 ## Load  Swiss banknotes data
 data(banknote)
 x <- banknote[,-1]

 ## Initial clusters with default arguments
 init1 <- InitClust(data=x, G=2)
 print(init1)

 ## Perform otrimle
 a <- otrimle(data=x, initial=init1$cluster)
 plot(a, what="clustering", data=x)
 }

