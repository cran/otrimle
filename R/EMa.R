.EMa <- function(data, initial, logicd, npr.max, erc, det.min, beta, em.iter.max, 
  em.tol, opt.selector = FALSE, gss.operation = TRUE) {
  N <- dim(data)[1L]
  P <- dim(data)[2L]
  G <- dim(initial)[2L] - 1
  icd <- exp(logicd)
  lmt <- log(det.min)
  tau.old <- initial
  sumtau.old <- .colSums(tau.old, m = N, n = {
    G + 1
  }, na.rm = TRUE)
  gausscost <- {
    2 * pi
  }^{
    -0.5 * P
  }
  UB1 <- 1.01 + beta * npr.max
  UB2 <- 2.02 + beta * npr.max
  STOP <- em.failed <- FALSE
  criterion <- 0
  flag <- rep(FALSE, times = 5)
  iter <- 0L
  iloglik.old <- iloglik.new <- -.Machine$double.xmax
  kd <- rep(0, G)
  pr <- sumtau.new <- rep(0, times = {
    1 + G
  })
  psi <- tau.new <- wts <- rep(0, times = N)
  smd <- matrix(0, nrow = N, ncol = G, byrow = TRUE, dimnames = NULL)
  phi <- matrix(0, nrow = N, ncol = {
    G + 1
  }, byrow = TRUE, dimnames = NULL)
  phi[, 1] <- icd
  Xpp <- matrix(0, nrow = P, ncol = P, byrow = TRUE, dimnames = NULL)
  Delta <- data
  L <- M <- matrix(0, nrow = P, ncol = G, byrow = TRUE, dimnames = NULL)
  V <- array(0, dim = c(P, P, G), dimnames = NULL)
  E <- list(values = L[, 1], vectors = V[, , 1])
  if (gss.operation) {
    ans <- list(logicd = logicd, criterion = UB2, iloglik = iloglik.old, code = 3, 
      flag = NA)
  }
  else {
    ans <- list(code = 3, flag = NA)
  }
  while (!STOP) {
    pr <- sumtau.old/N
    for (j in 1:G) {
      wts <- tau.old[, 1 + j]/sumtau.old[1 + j]
      M[, j] <- .colSums(wts * data, m = N, n = P, na.rm = TRUE)
      Delta <- {
        data - matrix(M[, j], nrow = N, ncol = P, byrow = TRUE, dimnames = NULL)
      }
      Xpp <- crossprod(sqrt(wts) * Delta, y = NULL)
      E <- eigen(x = Xpp, symmetric = TRUE)
      L[, j] <- E$values
      V[, , j] <- E$vectors
      if ({
        prod(L[, j]) < det.min
      } | any(L[, j] <= 0) | any(!is.finite(L[, j]))) {
        if (length(unique(L[, j])) == 1) {
          L[, j] <- rep(det.min^(1/P), P)
        }
        else {
          L[, j][L[, j] < .Machine$double.xmin] <- .Machine$double.xmin
          k <- 1
          while ((k <= (P - 1)) & (prod(L[, j]) < det.min)) {
          tmp <- lmt - sum(log(L[, j][1:(P - k)])) - k * log(L[, j][P - 
            k + 1])
          tmp <- exp(tmp/k)
          if (tmp > 0 & tmp <= 1) {
            L[, j][(P - k + 1):P] <- tmp * L[, j][P - k + 1]
          }
          k <- k + 1
          }
          if (prod(L[, j]) < det.min) {
          L[, j] <- rep(det.min^(1/P), P)
          }
        }
      }
      Xpp <- V[, , j] %*% {
        {
          1/L[, j]
        } * t(V[, , j])
      }
      smd[, j] <- .rowSums({
        Delta %*% Xpp
      } * Delta, m = N, n = P, na.rm = TRUE)
      phi[, 1 + j] <- pr[1 + j] * gausscost * {
        {
          prod(L[, j])
        }^{
          -0.5
        }
      } * {
        exp(-0.5 * smd[, j])
      }
    }
    phi[, 1] <- pr[1] * icd
    psi <- .rowSums(phi, m = N, n = {
      G + 1
    }, na.rm = TRUE)
    tau.new <- phi * matrix(1/psi, nrow = N, ncol = {
      1 + G
    }, byrow = FALSE, dimnames = NULL)
    sumtau.new <- .colSums(tau.new, m = N, n = {
      G + 1
    }, na.rm = TRUE)
    if (any(!is.finite(tau.new)) | any(sumtau.new[-1] == 0)) {
      em.failed <- TRUE
      flag[1] <- TRUE
    }
    if (!em.failed) {
      tmp <- (.rowSums(tau.new, m = N, n = {
        G + 1
      }, na.rm = TRUE) == 0)
      if (any(tmp)) {
        tau.new[tmp, ] <- 0
        if (icd > 0) {
          tau.new[tmp, 1] <- 1
        }
        if (icd == 0) {
          i_tmp <- which(tmp)
          n_tmp <- length(i_tmp)
          ED <- matrix(0, nrow = n_tmp, ncol = G, byrow = TRUE, dimnames = NULL)
          for (j in 1:G) {
          tmpDelta <- data[i_tmp, ] - matrix(M[, j], byrow = TRUE, nrow = n_tmp, 
            ncol = P, dimnames = NULL)
          ED[, j] <- .rowSums(tmpDelta * tmpDelta, m = n_tmp, n = P, na.rm = TRUE)
          }
          cl_tmp <- 1 + apply(-ED, 1, .rwhich.max)
          for (k in 1:n_tmp) {
          tau.new[i_tmp[k], cl_tmp[k]] <- 1
          }
        }
        sumtau.new <- .colSums(tau.new, m = N, n = {
          G + 1
        }, na.rm = TRUE)
      }
    }
    if (!em.failed & icd > 0 & {
      sumtau.new[1]/N > npr.max
    }) {
      z <- NULL
      try(z <- uniroot(f = .EqNPC, icd = icd, sumtau.old = sumtau.old, pr = pr, 
        psi = psi, N = N, npr.max = npr.max, lower = 0, upper = npr.max, 
        tol = .Machine$double.eps, maxiter = 500, trace = 0), silent = TRUE)
      if (is.null(z)) {
        flag[2] <- TRUE
        em.failed <- TRUE
      }
      else {
        flag[4] <- TRUE
        for (j in 1:G) {
          phi[, 1 + j] <- phi[, 1 + j]/pr[1 + j]
        }
        pr <- c(z$root, {
          {
          1 - z$root
          }/{
          N - sumtau.old[1]
          }
        } * sumtau.old[-1])
        for (j in 1:G) {
          phi[, 1 + j] <- pr[1 + j] * phi[, 1 + j]
        }
        phi[, 1] <- pr[1] * icd
        psi <- .rowSums(phi, m = N, n = {
          G + 1
        }, na.rm = TRUE)
        tau.new <- phi * matrix(1/psi, nrow = N, ncol = {
          1 + G
        }, byrow = FALSE, dimnames = NULL)
        sumtau.new <- .colSums(tau.new, m = N, n = {
          G + 1
        }, na.rm = TRUE)
        anpr <- mean({
          apply(tau.new, 1, .rwhich.max) - 1
        } == 0)
        if (anpr > {
          npr.max + 0.01
        }) {
          flag[2] <- TRUE
          em.failed <- TRUE
        }
        if (any(!is.finite(tau.new)) | any(sumtau.new[-1] == 0)) {
          em.failed <- TRUE
          flag[1] <- TRUE
        }
      }
    }
    if (em.failed) {
      STOP <- TRUE
    }
    else {
      iter <- iter + 1L
      iloglik.new <- sum(log(psi))
      delta.iloglik <- abs(iloglik.new - iloglik.old)
      if (is.finite(delta.iloglik) & {
        delta.iloglik < em.tol
      }) {
        STOP <- TRUE
      }
      if (iter >= em.iter.max) {
        STOP <- TRUE
      }
      iloglik.old <- iloglik.new
      tau.old <- tau.new
      sumtau.old <- sumtau.new
    }
  }
  if ({
    pr[1] * icd > 0
  } & {
    erc < Inf
  } & !em.failed) {
    Lmin <- min(L)
    Lmax <- max(L)
    EigenRatio <- Lmax/Lmin
    if ({
      EigenRatio > erc
    } | !is.finite(EigenRatio)) {
      EqLow <- .EqERC(Lmin, data = data, icd = icd, M = M, L = L, V = V, pr = pr, 
        erc = erc, npr.max = npr.max, N = N, G = G, P = P)
      if (!is.finite(EqLow)) {
        em.failed <- TRUE
        flag[3] <- TRUE
      }
      if (!em.failed & EqLow > 0) {
        qEq <- quantile(L, probs = seq(0.1, 1, by = 0.1), names = FALSE)
        k <- 1
        EqUp <- EqLow
        Lupper <- NULL
        while ({
          k <= 10
        } & is.null(Lupper)) {
          EqUp <- .EqERC(x = qEq[k], data = data, icd = icd, M = M, L = L, 
          V = V, pr = pr, erc = erc, npr.max = npr.max, N = N, G = G, P = P)
          if (is.finite(EqUp) & {
          EqUp <= 0
          }) {
          Lupper <- qEq[k]
          }
          else {
          k <- k + 1L
          }
        }
        if (is.null(Lupper)) {
          em.failed <- TRUE
          flag[3] <- TRUE
        }
        else {
          z <- NULL
          try(z <- uniroot(f = .EqERC, lower = Lmin, upper = Lupper, tol = .Machine$double.eps, 
          maxiter = 250, trace = 0, data = data, icd = icd, L = L, V = V, 
          pr = pr, erc = erc, npr.max = npr.max, N = N, G = G, P = P), 
          silent = TRUE)
          if (is.null(z)) {
          em.failed <- TRUE
          flag[3] <- TRUE
          }
          else {
          Lmin <- z$root
          }
        }
      }
      if (!em.failed) {
        L[L <= Lmin] <- Lmin
        L[L >= Lmin * erc] <- Lmin * erc
        for (j in 1:G) {
          Xpp <- V[, , j] %*% {
          {
            1/L[, j]
          } * t(V[, , j])
          }
          Delta <- {
          data - matrix(M[, j], nrow = N, ncol = P, byrow = TRUE, dimnames = NULL)
          }
          smd[, j] <- .rowSums({
          Delta %*% Xpp
          } * Delta, m = N, n = P, na.rm = TRUE)
          phi[, 1 + j] <- pr[1 + j] * gausscost * {
          {
            prod(L[, j])
          }^{
            -0.5
          }
          } * {
          exp(-0.5 * smd[, j])
          }
        }
        psi <- .rowSums(phi, m = N, n = {
          G + 1
        }, na.rm = TRUE)
        iloglik.old <- sum(log(psi))
        tau.old <- phi * matrix(1/psi, nrow = N, ncol = {
          1 + G
        }, byrow = FALSE, dimnames = NULL)
        sumtau.old <- .colSums(tau.old, m = N, n = {
          G + 1
        }, na.rm = TRUE)
        if (any(!is.finite(tau.old)) | any(sumtau.old[-1] == 0) | !is.finite(iloglik.old)) {
          em.failed <- TRUE
          flag[3] <- TRUE
        }
        else {
          flag[5] <- TRUE
        }
      }
    }
  }
  if ({
    pr[1] * icd == 0
  } & {
    erc < Inf
  } & !em.failed) {
    Lmin <- min(L)
    Lmax <- max(L)
    EigenRatio <- Lmax/Lmin
    if ({
      EigenRatio > erc
    } | !is.finite(EigenRatio)) {
      L[L <= Lmin] <- Lmin
      L[L >= Lmin * erc] <- Lmin * erc
      for (j in 1:G) {
        Xpp <- V[, , j] %*% {
          {
          1/L[, j]
          } * t(V[, , j])
        }
        Delta <- {
          data - matrix(M[, j], nrow = N, ncol = P, byrow = TRUE, dimnames = NULL)
        }
        smd[, j] <- .rowSums({
          Delta %*% Xpp
        } * Delta, m = N, n = P, na.rm = TRUE)
        phi[, 1 + j] <- pr[1 + j] * gausscost * {
          {
          prod(L[, j])
          }^{
          -0.5
          }
        } * {
          exp(-0.5 * smd[, j])
        }
      }
      psi <- .rowSums(phi, m = N, n = {
        G + 1
      }, na.rm = TRUE)
      iloglik.old <- sum(log(psi))
      tau.old <- phi * matrix(1/psi, nrow = N, ncol = {
        1 + G
      }, byrow = FALSE, dimnames = NULL)
      sumtau.old <- .colSums(tau.old, m = N, n = {
        G + 1
      }, na.rm = TRUE)
      if (any(!is.finite(tau.old)) | any(sumtau.old[-1] == 0) | !is.finite(iloglik.old)) {
        em.failed <- TRUE
        flag[3] <- TRUE
      }
      else {
        flag[5] <- TRUE
      }
    }
  }
  if (!any(flag[c(1, 2, 3)])) {
    for (j in 1:G) {
      kd[j] <- max(abs(.wecdf(x = smd[, j], weights = tau.old[, 1 + j]) - pchisq(smd[, 
        j], df = P)))
    }
    criterion <- (sum(kd * pr[-1])/sum(pr[-1])) + {
      beta * pr[1]
    }
  }
  if (any(flag[c(1, 2, 3)])) {
    ans$code <- 0
  }
  else if (iter == em.iter.max) {
    ans$code <- 1
  }
  else {
    ans$code <- 2
  }
  if (any(flag)) {
    ans$flag <- paste(as.character(which(flag)), collapse = ",")
  }
  else {
    ans$flag <- "None"
  }
  if (gss.operation) {
    if (ans$code == 0) {
      ans$criterion <- UB2
      ans$iloglik <- NA
    }
    else if (opt.selector & {
      flag[4] | flag[5]
    }) {
      ans$criterion <- criterion + UB1
      ans$iloglik <- iloglik.old
    }
    else if (opt.selector & !flag[4] & !flag[5]) {
      ans$criterion <- criterion
      ans$iloglik <- iloglik.old
    }
    else if (!opt.selector) {
      ans$criterion <- criterion
      ans$iloglik <- iloglik.old
    }
  }
  else if (ans$code > 0) {
    ans$iter <- iter
    ans$logicd <- logicd
    ans$iloglik <- iloglik.old
    ans$criterion <- criterion
    ans$npr <- pr[1]
    ans$cpr <- pr[-1]
    ans$mean <- M
    ans$cov <- V
    ans$tau <- tau.old
    ans$smd <- smd
    ans$cluster <- apply(tau.old, 1, .rwhich.max) - 1
    ans$size <- rep(0, 1 + G)
    ans$size[1] <- sum(ans$cluster == 0)
    for (j in 1:G) {
      ans$cov[, , j] <- V[, , j] %*% {
        L[, j] * t(V[, , j])
      }
      ans$size[1 + j] <- sum(ans$cluster == j)
    }
    VarNames <- colnames(data)
    ClustName <- paste("Cluster.", 1:G, sep = "")
    names(ans$cpr) <- colnames(ans$smd) <- colnames(ans$mean) <- ClustName
    dimnames(ans$cov)[[3]] <- as.list(ClustName)
    names(ans$size) <- colnames(ans$tau) <- c("Noise", ClustName)
    if (!is.null(VarNames)) {
      rownames(ans$mean) <- dimnames(ans$cov)[[1]] <- dimnames(ans$cov)[[2]] <- VarNames
    }
  }
  else if (ans$code == 0) {
    ans$iter <- iter
  }
  return(ans)
}
