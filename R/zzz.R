Rcpp::loadModule(module = "ds", what = TRUE)

.onUnload <- function(libpath) {
  library.dynam.unload("DirichletStudy", libpath)
}