"ppml" <- function(formula = formula, data = data, ...){

  ## use glm function to estimate gravity model
  model <- glm(formula, family = quasipoisson(link = "log"), data = data, ...)

  class(model) <- c("ppml", "glm")
  return(model)
}

