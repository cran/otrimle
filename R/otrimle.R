otrimle <- function(data, G, initial = NULL, logicd = NULL, npr.max = 0.5, erc = 50, 
  det.min = .Machine$double.eps, beta = 0, opt.selector = FALSE, cmstep = TRUE, 
  iter.max = 100, tol = 1e-04, em.iter.max = 500, em.tol = 1e-06, monitor = 1) {
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
    !is.numeric(logicd) | length(logicd) != 2
  }) {
    stop("\"logicd\" must be a vector of length=2 with logicd[1]<logicd[2].")
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
  if (iter.max < 2) {
    stop("\"iter.max\" must be larger or equal to 2.")
  }
  if (em.iter.max < 1) {
    stop("\"em.iter.max\" must be larger or equal to 1.")
  }
  if (tol < 0) {
    stop("\"tol\" must be positive.")
  }
  if (em.tol < 0) {
    stop("\"em.tol\" must be positive.")
  }
  if (monitor > 0) {
    message("\nOTRIMLE:")
  }
  if (cmstep) {
    doEM <- .EMb
  }
  else {
    doEM <- .EMa
  }
  if (erc == Inf) {
    opt.selector <- FALSE
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
  UB1 <- 1.01 + beta * npr.max
  grc <- {
    sqrt(5) - 1
  }/2
  niter <- iter.max - 2
  optimize.again <- FALSE
  if (is.null(logicd)) {
    N <- nrow(data)
    P <- ncol(data)
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
    x1 <- a + {
      1 - grc
    } * {
      b - a
    }
    x2 <- a + grc * {
      b - a
    }
    logicd.init <- c(a, x1, x2, b)
  }
  else {
    a <- logicd[1]
    b <- logicd[2]
    x1 <- a + {
      1 - grc
    } * {
      b - a
    }
    x2 <- a + grc * {
      b - a
    }
    logicd.init <- c(a, x1, x2, b)
  }
  if (opt.selector) {
    if (monitor >= 1) {
      message(".......: searching an interior solution")
    }
    a <- logicd.init[1]
    x1 <- logicd.init[2]
    x2 <- logicd.init[3]
    b <- logicd.init[4]
    iter <- 0L
    STOP <- FALSE
    A <- data.frame(logicd = rep(0, times = iter.max), criterion = rep(0, times = iter.max), 
      iloglik = rep(0, times = iter.max), code = rep(0L, times = iter.max), 
      flag = rep("2", times = iter.max), row.names = NULL, stringsAsFactors = FALSE)
    if (monitor == 1) {
      message(".......: searching for an interior solution")
    }
    if (monitor >= 2) {
      message(".......: evaluating criterion at log(icd) = ", -Inf)
    }
    iter <- iter + 1L
    x <- -Inf
    A[iter, ] <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
      erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
      em.tol = em.tol, opt.selector = TRUE, gss.operation = TRUE)
    if (A[iter, 2] >= UB1) {
      A[iter, c(2, 3)] <- NA
    }
    iter <- iter + 1L
    x <- x1
    if (monitor >= 2) {
      message(".......: evaluating criterion at log(icd) = ", format(x, scientific = FALSE, 
        digits = 10))
    }
    A[iter, ] <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
      erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
      em.tol = em.tol, opt.selector = TRUE, gss.operation = TRUE)
    f1 <- A[iter, 2]
    if (A[iter, 2] >= UB1) {
      A[iter, c(2, 3)] <- NA
    }
    iter <- iter + 1L
    x <- x2
    if (monitor >= 2) {
      message(".......: evaluating criterion at log(icd) = ", format(x, scientific = FALSE, 
        digits = 10))
    }
    A[iter, ] <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
      erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
      em.tol = em.tol, opt.selector = TRUE, gss.operation = TRUE)
    f2 <- A[iter, 2]
    if (A[iter, 2] >= UB1) {
      A[iter, c(2, 3)] <- NA
    }
    while ({
      abs(b - a) > tol
    } & {
      iter < niter
    }) {
      if (f1 > f2) {
        a <- x1
        x1 <- x2
        x2 <- a + grc * {
          b - a
        }
        f1 <- f2
        iter <- iter + 1L
        x <- x2
        if (monitor >= 2) {
          message(".......: evaluating criterion at log(icd) = ", format(x, 
          scientific = FALSE, digits = 10))
        }
        A[iter, ] <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
          erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
          em.tol = em.tol, opt.selector = TRUE, gss.operation = TRUE)
        f2 <- A[iter, 2]
        if (A[iter, 2] >= UB1) {
          A[iter, c(2, 3)] <- NA
        }
      }
      else {
        b <- x2
        x2 <- x1
        x1 <- a + {
          1 - grc
        } * {
          b - a
        }
        f2 <- f1
        iter <- iter + 1L
        x <- x1
        if (monitor >= 2) {
          message(".......: evaluating criterion at log(icd) = ", format(x, 
          scientific = FALSE, digits = 10))
        }
        A[iter, ] <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
          erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
          em.tol = em.tol, opt.selector = TRUE, gss.operation = TRUE)
        f1 <- A[iter, 2]
        if (A[iter, 2] >= UB1) {
          A[iter, c(2, 3)] <- NA
        }
      }
    }
    A <- A[1:iter, ]
    A <- A[order(A[, 1]), ]
    rownames(A) <- 1:iter
    if (all(is.na(A[1:iter, 2]))) {
      optimize.again <- TRUE
      if (monitor >= 1) {
        message("\n.......: failed to find an interior solution\n")
        message(".......: searching for a border solution")
      }
    }
    else {
      x <- A[which.min(A[1:iter, 2]), 1]
      if (monitor >= 1) {
        message(".......: computing the optimal interior solution")
      }
      ans <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
        erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
        em.tol = em.tol, opt.selector = TRUE, gss.operation = FALSE)
      ans$optimization <- A
    }
  }
  if (!opt.selector | optimize.again) {
    a <- logicd.init[1]
    x1 <- logicd.init[2]
    x2 <- logicd.init[3]
    b <- logicd.init[4]
    iter <- 0L
    STOP <- FALSE
    A <- data.frame(logicd = rep(0, times = iter.max), criterion = rep(0, times = iter.max), 
      iloglik = rep(0, times = iter.max), code = rep(0L, times = iter.max), 
      flag = rep("2", times = iter.max), row.names = NULL, stringsAsFactors = FALSE)
    if (monitor >= 2) {
      message(".......: evaluating criterion at log(icd) = ", -Inf)
    }
    iter <- iter + 1L
    x <- -Inf
    A[iter, ] <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
      erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
      em.tol = em.tol, opt.selector = FALSE, gss.operation = TRUE)
    if (A[iter, 2] >= UB1) {
      A[iter, c(2, 3)] <- NA
    }
    iter <- iter + 1L
    x <- x1
    if (monitor >= 2) {
      message(".......: evaluating criterion at log(icd) = ", format(x, scientific = FALSE, 
        digits = 10))
    }
    A[iter, ] <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
      erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
      em.tol = em.tol, opt.selector = FALSE, gss.operation = TRUE)
    f1 <- A[iter, 2]
    if (A[iter, 2] >= UB1) {
      A[iter, c(2, 3)] <- NA
    }
    iter <- iter + 1L
    x <- x2
    if (monitor >= 2) {
      message(".......: evaluating criterion at log(icd) = ", format(x, scientific = FALSE, 
        digits = 10))
    }
    A[iter, ] <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
      erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
      em.tol = em.tol, opt.selector = FALSE, gss.operation = TRUE)
    f2 <- A[iter, 2]
    if (A[iter, 2] >= UB1) {
      A[iter, c(2, 3)] <- NA
    }
    while ({
      abs(b - a) > tol
    } & {
      iter < niter
    }) {
      if (f1 > f2) {
        a <- x1
        x1 <- x2
        x2 <- a + grc * {
          b - a
        }
        f1 <- f2
        iter <- iter + 1L
        x <- x2
        if (monitor >= 2) {
          message(".......: evaluating criterion at log(icd) = ", format(x, 
          scientific = FALSE, digits = 10))
        }
        A[iter, ] <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
          erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
          em.tol = em.tol, opt.selector = FALSE, gss.operation = TRUE)
        f2 <- A[iter, 2]
        if (A[iter, 2] >= UB1) {
          A[iter, c(2, 3)] <- NA
        }
      }
      else {
        b <- x2
        x2 <- x1
        x1 <- a + {
          1 - grc
        } * {
          b - a
        }
        f2 <- f1
        iter <- iter + 1L
        x <- x1
        if (monitor >= 2) {
          message(".......: evaluating criterion at log(icd) = ", format(x, 
          scientific = FALSE, digits = 10))
        }
        A[iter, ] <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
          erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
          em.tol = em.tol, opt.selector = FALSE, gss.operation = TRUE)
        f1 <- A[iter, 2]
        if (A[iter, 2] >= UB1) {
          A[iter, c(2, 3)] <- NA
        }
      }
    }
    A <- A[1:iter, ]
    A <- A[order(A[, 1]), ]
    rownames(A) <- 1:iter
    if (all(is.na(A[1:iter, 2]))) {
      ans <- list()
      ans$code <- 0
      ans$optimization <- A
    }
    else {
      x <- A[which.min(A[1:iter, 2]), 1]
      if (monitor >= 1) {
        message(".......: computing the optimal solution")
      }
      ans <- doEM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
        erc = erc, det.min = det.min, beta = beta, em.iter.max = em.iter.max, 
        em.tol = em.tol, opt.selector = TRUE, gss.operation = FALSE)
      ans$optimization <- A
    }
  }
  class(ans) <- "otrimle"
  return(ans)
}
print.otrimle <- function(x, ...) {
  if (x$code > 0) {
    G <- ncol(x$mean)
    if (x$size[1] > 0) {
      cat("\nOTRIMLE with optimal logicd=", round(x$logicd, 6), " discovered G=", 
        G, " clusters plus noise of size:", "\n\n", sep = "")
      print(x$size)
      cat("\nAvailable components:\n", sep = "")
      print(names(x))
    }
    else {
      cat("\nOTRIMLE with optimal logicd=", round(x$logicd, 6), " discovered G=", 
        G, " clusters (without noise) \nof size:", "\n\n", sep = "")
      print(x$size[-1])
      cat("\nAvailable components:\n")
      print(names(x))
    }
  }
  else {
    cat("\nOTRIMLE could not find a solution based on the current settings. Please check flags below and adjust settings as suggested in \"Details\" of help(\"otrimle\").", 
      sep = "")
    cat("\n\nOptimization profiling:\n ", sep = "\n")
    print(x$optimization)
  }
}
