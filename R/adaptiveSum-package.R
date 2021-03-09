#' @useDynLib adaptiveSum, .registration = TRUE


.onUnload <- function (libpath) {
  library.dynam.unload("adaptiveSum", libpath)
}
