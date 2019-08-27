"bread.ppml" <- function(object, ...) {
  class(object) <- c("glm", "lm")
  NextMethod()
}
