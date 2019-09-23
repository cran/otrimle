print.otrimle <- function(x, ...) {
  if (x$code > 0) {
    G <- ncol(x$mean)
    if (x$size[1] > 0) {
      cat("\nOTRIMLE with optimal logicd=", round(x$logicd, 6), " discovered G=", 
        G, " clusters plus noise of size:", "\n\n", sep = "")
      print(x$size)
      cat("\nAvailable components:\n")
      cat(paste(as.character(names(x)), collapse = ", "))
    }
    else {
      cat("\nOTRIMLE with optimal logicd=", round(x$logicd, 6), " discovered G=", 
        G, " clusters (without noise) \nof size:", "\n\n", sep = "")
      print(x$size[-1])
      cat("\nAvailable components:\n")
      cat(paste(as.character(names(x)), collapse = ", "))
    }
  }
  else {
    cat("\nOTRIMLE could not find a solution based on the current settings. Please check flags below and adjust settings as suggested in \"Details\" of help(\"otrimle\").", 
      sep = "")
    cat("\n\nOptimization profiling:\n ", sep = "\n")
    print(x$optimization)
  }
}
print.rimle <- function(x, ...) {
  if (x$code > 0) {
    G <- ncol(x$mean)
    if (x$size[1] > 0) {
      cat("\nRIMLE with fixed logicd=", round(x$logicd, 6), " discovered G=", 
        G, " clusters plus noise of size:", "\n\n", sep = "")
      print(x$size)
      cat("\nAvailable components:\n")
      cat(paste(as.character(names(x)), collapse = ", "))
    }
    else {
      cat("\nRIMLE with fixed logicd=", round(x$logicd, 6), " discovered G=", 
        G, " clusters (without noise) of size:", "\n\n", sep = "")
      print(x$size[-1])
      cat("\nAvailable components:\n")
      cat(paste(as.character(names(x)), collapse = ", "))
    }
  }
  else {
    cat("\nRIMLE with fixed logicd=", round(x$logicd, 6), " could not find a solution based on the current settings. Please check flags below, and adjust settings as suggested in \"Details\" of help(\"rimle\").", 
      sep = "")
    cat(paste("\n\nFlags: ", x$flag, sep = ""))
  }
}
