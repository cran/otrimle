otrimle <- function(data, G, initial = NULL, logicd = NULL, npr.max = 0.5, erc = 50, 
  iter.max = 500, tol = 1e-06, ncores = NULL, monitor = TRUE) {
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
  G <- as.integer(G)
  if (G < 1) {
    stop("G must be an integer greater or equal to one.")
  }
  if (!is.null(initial)) {
    if (G != max(initial)) {
      stop("G is not equal to max(initial).")
    }
  }
  if (!is.null(logicd) & {
    !is.numeric(logicd) | length(logicd) <= 1
  }) {
    stop("\"logicd\" must be a vector with length(logicd)>=2.")
  }
  if (npr.max <= 0 | npr.max >= 1) {
    stop("\"npr.max\" must belog to the interval (0,1). Tipically npr.max<=0.5.")
  }
  if (erc < 1) {
    stop("\"ecr\" must be larger or equal to 1.")
  }
  if (iter.max < 2) {
    stop("\"iter.max\" must be an integer larger or equal to 2.")
  }
  if (tol < 0) {
    stop("\"tol\" must be positive.")
  }
  if (is.null(ncores)) {
    DetectedCores <- detectCores()
    if (DetectedCores <= 2) {
      ncores <- 1
    }
    else {
      ncores <- {
        DetectedCores - 1
      }
    }
  }
  else {
    ncores <- as.integer(ncores)
    if (ncores <= 0) {
      stop("\"ncores\" must be an integer larger or equal to 1.")
    }
  }
  if (!is.null(logicd)) {
    logicd <- sort(logicd)
  }
  else {
    logicd <- c(-Inf, seq(-700, -100, by = 50), seq(-95, -55, by = 5), seq(-50, 
      -10, by = 2.5), seq(-9, 0, by = 1))
  }
  nGrid <- length(logicd)
  if (monitor) {
    message("\nOTRIMLE:")
  }
  if (is.null(initial)) {
    if (monitor) {
      message("........computing initial clusters")
    }
    initial <- .Cluster2Assign(InitClust(data = data, G = G, knnd.trim = npr.max))
  }
  else {
    initial <- .Cluster2Assign(initial)
  }
  G <- ncol(initial) - 1
  nj <- rep(0, times = G)
  for (j in 1:G) {
    nj[j] <- .CountUniqueRows(data[initial[, 1 + j] == 1, ])
  }
  if (any(nj < 2)) {
    Suggest <- paste("initial = sample(0:", G, ",  size=", nrow(data), ", replace=TRUE)", 
      sep = "")
    StopMsg <- paste("\n\nSome of the initial clusters contain less then 2 distinct points. Initial parameters cannot be computed. Try a random initial partition setting: \n\n", 
      Suggest, sep = "")
    stop(StopMsg)
  }
  A <- data.frame(logicd = logicd, criterion = NA, iloglik = NA, npr = NA, code = NA, 
    flag = NA, row.names = NULL, stringsAsFactors = FALSE)
  rownames(A) <- 1:nGrid
  registerDoParallel(ncores)
  if (monitor) {
    message("........searching for the optimal log(icd)")
  }
  B <- foreach(k = 1:nGrid) %dopar% {
    .ECM(data = data, initial = initial, logicd = A[k, 1], npr.max = npr.max, 
      erc = erc, iter.max = iter.max, tol = tol, grid.operation = TRUE)
  }
  for (k in 1:nGrid) {
    A[k, ] <- B[[k]]
  }
  if (all(is.na(A[, 2]))) {
    ans <- list()
    ans$code <- 0
    ans$optimization <- A
  }
  else {
    x <- A[which.min(A[, 2]), 1]
    if (monitor) {
      message("........computing the optimal solution\n")
    }
    ans <- .ECM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
      erc = erc, iter.max = iter.max, tol = tol, grid.operation = FALSE)
    ans$optimization <- A
  }
  class(ans) <- "otrimle"
  if (monitor) {
    cat(print(ans))
  }
  return(ans)
}
