summary.wAF<-function(x, ...){
  cat("Method:\n")
  cat(paste(x$method, "\n"))
  cat("\n")
  cat("P-value:\n")
  print(x$pv)
  cat("\n")
  cat("Weights:\n")
  if (x$weight_method == "flat") {
    cat("Flat Weights")
  } else {
    cat(paste(x$weight_method, "weights: \n"))
    print(x$weight_vector)
  }
}

print.wAF<-function(x, ...){
  cat("Method:\n")
  cat(paste(x$method, "\n"))
  cat("\n")
  cat("P-value:\n")
  print(x$pv)
  cat("\n")
  cat("Weights:\n")
  if (x$weight_method == "flat") {
    cat("Flat Weights")
  } else {
    cat(paste(x$weight_method, "weights: \n"))
    print(x$weight_vector)
  }
}