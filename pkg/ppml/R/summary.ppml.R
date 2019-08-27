"summary.ppml" <- function(object, ...){
  robust <- coeftest(object, sandwich(object))
  summary <- summary.glm(object)
  summary$coefficients <- robust 
  return(summary)
}
