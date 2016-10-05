InitClust <- function(data, G, cpr.min = {
  ncol(data) + 1
}/nrow(data), K = 5, nstart.km = 50, modelName = "VVV", monitor = FALSE) {
  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  }
  if (!is.matrix(data)) {
    data <- data.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("data must be numeric.")
  }
  if (cpr.min < 2/nrow(data)) {
    stop("cpr.min > 2/nrow(data) is required.")
  }
  ndistinct <- .CountUniqueRows(data)
  if (ndistinct < {
    2 * G
  }) {
    stop("data contains ", ndistinct, " distinct points, the maximum value for G is ", 
      floor(ndistinct/2))
  }
  N <- nrow(data)
  P <- ncol(data)
  min.nc <- max(2, ceiling(cpr.min * N))
  cl <- rep(0, N)
  nj <- mj <- dj <- rep(0, G)
  do.random <- FALSE
  do.next <- FALSE
  code <- 1
  if (monitor) {
    message("\n<> Perform nearest neighbor based noise detection")
  }
  nncl <- NNclean(data, k = K)$z
  idx.noise <- which(nncl == 0)
  idx.reg <- which(nncl != 0)
  ndistinct <- .CountUniqueRows(data[idx.reg, ])
  if (ndistinct < {
    2 * G
  }) {
    do.random <- TRUE
  }
  if (!do.random) {
    if (monitor) {
      message("<> Perform model-based hierarchical clustering [MBHC]")
    }
    Hclust <- hc(modelName = modelName, data = data[idx.reg, ], minclus = G)
    cl[idx.reg] <- as.vector(hclass(Hclust, G))
    for (j in 1:G) {
      nj[j] <- sum(cl == j)
      mj[j] <- .CountUniqueRows(data[cl == j, ])
    }
    if (any(mj < 2) | any(nj < min.nc)) {
      do.next <- TRUE
    }
  }
  if (do.next) {
    if (monitor) {
      message("<> Small clusters found in the previous step. Assign small clusters to noise and perform MBHC on denoised data")
    }
    code <- 2
    badcl <- which({
      nj < min.nc
    } | {
      mj < 2
    } | {
      dj <= 0
    })
    badx <- which(cl == badcl[1])
    if (length(badcl) > 1) {
      for (i in 2:length(badcl)) {
        badx <- c(badx, which(cl == badcl[i]))
      }
    }
    idx.noise <- c(badx, idx.noise)
    cl[idx.noise] <- 0
    idx.reg <- idx.reg[-badx]
    ndistinct <- .CountUniqueRows(data[idx.reg, ])
    if (ndistinct < {
      2 * G
    }) {
      do.random <- TRUE
    }
    else {
      Hclust <- hc(modelName = modelName, data = data[idx.reg, ], minclus = G)
      cl[idx.reg] <- as.vector(hclass(Hclust, G))
      for (j in 1:G) {
        nj[j] <- sum(cl == j)
        mj[j] <- .CountUniqueRows(data[cl == j, ])
      }
    }
    if (any(mj < 2) | any(nj < min.nc)) {
      do.next <- TRUE
    }
  }
  if (do.next) {
    if (monitor) {
      message("<> Small clusters found in the previous step. Assign small cluster to noise and perform k-means on the denoised data")
    }
    code <- 3
    badcl <- which({
      nj < min.nc
    } | {
      mj < 2
    } | {
      dj <= 0
    })
    badx <- which(cl == badcl[1])
    if (length(badcl) > 1) {
      for (i in 2:length(badcl)) {
        badx <- c(badx, which(cl == badcl[i]))
      }
    }
    idx.noise <- c(badx, idx.noise)
    cl[idx.noise] <- 0
    idx.reg <- idx.reg[-badx]
    ndistinct <- .CountUniqueRows(data[idx.reg, ])
    if (ndistinct < {
      2 * G
    }) {
      do.random <- TRUE
    }
    else {
      cl[idx.reg] <- as.vector(kmeans(x = data[idx.reg, ], centers = G, nstart = nstart.km)$cluster)
      for (j in 1:G) {
        nj[j] <- sum(cl == j)
        mj[j] <- .CountUniqueRows(data[cl == j, ])
      }
    }
    if (any(mj < 2) | any(nj < min.nc)) {
      do.random <- TRUE
    }
  }
  if (do.random) {
    if (monitor) {
      message("<> Small clusters found in previous steps. Perform random initialization")
    }
    code <- 4
    StopWhile <- FALSE
    iter <- 1L
    while (!StopWhile) {
      if (iter > 1000) {
        stop("impossible to define G valid initial clusters. Try to reduce G and/or cpr.min.")
      }
      cl <- sample(1:G, size = N, replace = TRUE)
      for (j in 1:G) {
        nj[j] <- sum(cl == j)
        mj[j] <- .CountUniqueRows(data[cl == j, ])
      }
      if (all(nj >= min.nc) & all(mj >= 2)) {
        StopWhile <- TRUE
      }
      iter <- 1L + iter
    }
  }
  ans <- list()
  ans$code <- code
  ans$cluster <- cl
  return(ans)
}
.Cluster2Assign <- function(cluster) {
  cl <- sort(unique(cluster[cluster != 0]))
  G <- length(cl)
  N <- length(cluster)
  A <- matrix(0, nrow = N, ncol = 1 + G)
  if (sum(cluster == 0) > 0) {
    A[cluster == 0, 1] <- 1
  }
  for (j in 1:G) {
    A[cluster == cl[j], j + 1] <- 1
  }
  return(A)
}
