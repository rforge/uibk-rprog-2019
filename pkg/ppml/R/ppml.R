"ppml" <- function(formula = formula, data = data, robust = TRUE, ...){
  
  ## use glm function to estimate gravity model
  model <- glm(formula, family = quasipoisson(link = "log"), data = data)
  ## extract results of the glm object with summary
  results <- summary(model)
  
  ## calculate robust standard errors (HC1)
  if(robust == TRUE){
    
    vcovR <- vcovHC(model, type = "HC1")
    robust <- coeftest(model, vcovR)
    
    ## overwrite non robust standard errors
    results$coefficients <- robust 
    
  }
  
  ## return the estiamted results
  return(results)
  
}
