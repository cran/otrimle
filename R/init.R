InitClust <- function(data, G, k = 3, knnd.trim = 0.5, modelName = "VVV") {
  n <- nrow(data)
  distances <- dist(data)
  kthNND <- rep(0, n)
  Labels <- 1:(n - 1)
  kthNND[1] <- sort(distances[Labels])[k]
  Labels[(2):(n - 1)] <- Labels[(2):(n - 1)] + (n - 1 - 1)
  for (i in 2:n) {
    kthNND[i] <- sort(distances[Labels])[k]
    Labels[1:(i - 1)] <- Labels[1:(i - 1)] + 1
    Labels[(i + 1):(n - 1)] <- Labels[(i + 1):(n - 1)] + (n - i - 1)
  }
  kthNND <- kthNND[1:n]
  idx.reg <- {
    kthNND <= sort(kthNND)[n - floor(n * knnd.trim)]
  }
  H <- hc(data = data[idx.reg, ], minclus = G, modelName = modelName)
  Hcluster <- as.vector(hclass(H, G))
  cluster <- rep(0, n)
  cluster[idx.reg] <- Hcluster
  return(cluster)
}
