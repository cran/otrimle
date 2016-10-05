rimle <- function(data, G, initial = NULL, logicd = NULL, npr.max = 0.5, erc = 50, 
  det.min = .Machine$double.eps, cmstep = TRUE, em.iter.max = 500, em.tol = 1e-06, 
  monitor = 1) {
  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  }
  if (!is.matrix(data)) {
    data <- data.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("\"data\" must be numeric.")
  }
  if (!is.null(initial) & {
    !is.numeric(initial) | any(initial < 0) | {
      length(initial) != nrow(data)
    }
  }) {
    stop("\"initial\" must be an integer vector of length equal to nrow(data).")
  }
  if (!is.null(logicd) & {
    !is.numeric(logicd) | length(logicd) != 1
  }) {
    stop("\"logicd\" must be numeric, with length(logicd)=1.")
  }
  if (npr.max <= 0 | npr.max >= 1) {
    stop("\"npr.max\" must belog to the interval (0,1). Tipically npr.max<=0.5.")
  }
  if (erc < 1) {
    stop("\"ecr\" must be larger or equal to 1.")
  }
  if (det.min <= 0) {
    stop("\"det.min\" must be positive.")
  }
  if (em.iter.max < 1) {
    stop("\"em.iter.max\" must be larger or equal to 1.")
  }
  if (em.tol < 0) {
    stop("\"em.tol\" must be positive.")
  }
  if (monitor > 0) {
    message("\nRIMLE:")
  }
  if (cmstep) {
    doEM <- .EMb
  }
  else {
    doEM <- .EMa
  }
  if (is.null(initial)) {
    if (monitor >= 1) {
      message(".......: computing initial clusters")
    }
    initial <- .Cluster2Assign(InitClust(data = data, G = G)$cluster)
  }
  else {
    initial <- .Cluster2Assign(initial)
    G <- ncol(initial) - 1
    nj <- rep(0, times = G)
    for (j in 1:G) {
      nj[j] <- .CountUniqueRows(data[initial[, 1 + j] == 1, ])
    }
    if (any(nj < 2)) {
      stop("some of the initial clusters contain less then 2 distinct points. Initial parameters cannot be computed.")
    }
  }
  if (is.null(logicd)) {
    N <- nrow(data)
    P <- ncol(data)
    G <- ncol(initial) - 1
    cH <- -{
      P * log(2 * pi) + qchisq(0.5, df = P)
    }/2
    cL <- max(-{
      P * log(2 * pi) + qchisq({
        1 - 0.1
      }^{
        1/N
      }, df = P)
    }/2, log(.Machine$double.xmin))
    L <- matrix(0, nrow = P, ncol = G)
    for (j in 1:G) {
      idx <- initial[, 1 + j] == 1
      L[, j] <- eigen(cov(as.matrix(data[idx, ])), symmetric = TRUE, only.values = TRUE)$values
    }
    if (any(L <= 0) & erc == Inf) {
      L <- .GssERC(values = L, erc = 1e+12, Lmin = min(L), Lmax = max(L), sumtau = nj, 
        P = P, G = G)
    }
    else if (any(L <= 0) & erc < Inf) {
      L <- .GssERC(values = L, erc = erc, Lmin = min(L), Lmax = max(L), sumtau = nj, 
        P = P, G = G)
    }
    logdetj <- apply(log(L), 2, sum)
    a <- min(c(cL - logdetj/2, cH - logdetj/2))
    b <- max(c(cL - logdetj/2, cH - logdetj/2))
    logicd <- {
      a + b
    }/2
  }
  if (monitor >= 1) {
    message(".......: computing solution")
  }
  ans <- doEM(data = data, initial = initial, logicd = logicd, npr.max = npr.max, 
    erc = erc, det.min = det.min, beta = 0, em.iter.max = em.iter.max, em.tol = em.tol, 
    opt.selector = FALSE, gss.operation = FALSE)
  class(ans) <- "rimle"
  return(ans)
}
print.rimle <- function(x, ...) {
  if (x$code > 0) {
    G <- ncol(x$mean)
    if (x$size[1] > 0) {
      cat("\nRIMLE with fixed logicd=", round(x$logicd, 6), " discovered G=", 
        G, " clusters plus noise of size:", "\n\n", sep = "")
      print(x$size)
      cat("\nAvailable components:\n")
      print(names(x))
    }
    else {
      cat("\nRIMLE with fixed logicd=", round(x$logicd, 6), " discovered G=", 
        G, " clusters (without noise) of size:", "\n\n", sep = "")
      print(x$size[-1])
      cat("\nAvailable components:\n")
      print(names(x))
    }
  }
  else {
    cat("\nRIMLE with fixed logicd=", round(x$logicd, 6), " could not find a solution based on the current settings. Please check flags below, and adjust settings as suggested in \"Details\" of help(\"rimle\").", 
      sep = "")
    cat(paste("\n\nFlags: ", x$flag, sep = ""))
  }
}
