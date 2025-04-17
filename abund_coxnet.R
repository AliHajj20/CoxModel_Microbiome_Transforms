abund_coxnet <- function (x, time, status, covar = NULL, lambda = "lambda.1se", 
          nvar = NULL, alpha = 0.9, nfolds = 10, showPlots = TRUE, 
          coef_threshold = 0) 
{
  x <- impute_zeros(x)
  y = Surv(time, status)
  
  if (is.null(covar)) {
    cvfit <- glmnet::cv.glmnet(x, y, family = "cox", type.measure = "C", 
                               alpha = alpha, nfolds = nfolds, keep = TRUE)
  }
  else {
    df0 <- data.frame(as.matrix(y), covar)
    model0 <- coxph(Surv(time, status) ~ ., data = df0)
    x0 <- predict(model0)
    cvfit <- glmnet::cv.glmnet(x, y, family = "cox", type.measure = "C", 
                               nfolds = nfolds, alpha = alpha, keep = TRUE, offset = x0)
  }
  
  if (showPlots == TRUE) {
    plot(cvfit)
  }
  
  if (!is.null(nvar)) {
    rowlasso <- max(which(cvfit$glmnet.fit$df <= nvar))
    lambda <- cvfit$glmnet.fit$lambda[rowlasso]
  }
  
  lambdavalue <- lambda
  if (is.character(lambda)) {
    if (lambda == "lambda.1se") 
      lambdavalue <- cvfit$lambda.1se
    if (lambda == "lambda.min") 
      lambdavalue <- cvfit$lambda.min
  }
  
  idrow <- max(which(cvfit$glmnet.fit$lambda >= lambdavalue))
  coefs <- as.vector(coef(cvfit, s = lambda))
  selected_vars <- which(coefs != 0)
  coef_values <- coefs[selected_vars]
  
  names.select <- colnames(x)[selected_vars]
  sign <- ifelse(coef_values > 0, 1, 0)
  sign <- factor(sign, levels = c(0, 1), labels = c("negative", "positive"))
  
  if (is.null(covar)) {
    predictions <- as.numeric(predict(cvfit, x, s = lambdavalue))
  }
  else {
    predictions <- as.numeric(predict(cvfit, x, s = lambdavalue, newoffset = x0))
  }
  
  # Normalize coefficients (optional)
  coef_values <- 2 * coef_values/sum(abs(coef_values))
  
  if (length(selected_vars) == 0) {
    Cindex_signature <- 0.5
  }
  else {
    Cindex_signature <- glmnet::Cindex(pred = predictions, y)
  }
  
  mcvCindex <- cvfit$cvm[idrow]
  sdcvCindex <- cvfit$cvsd[idrow]
  plot1 <- NULL
  plot2 <- NULL
  
  if (length(selected_vars > 0)) {
    plot1 <- plot_riskscore(predictions, x, time, status, showPlots = showPlots)
    plot2 <- plot_signature(names.select, coef_values, showPlots = showPlots)
  }
  else {
    print("No variables are selected. The risk score plot and the signature plot are not displayed.")
  }
  
  results <- list(taxa.num = selected_vars, 
                  taxa.name = names.select, 
                  coefficients = coef_values, 
                  risk.score = predictions, 
                  `apparent Cindex` = Cindex_signature, 
                  `mean cv-Cindex` = mcvCindex, 
                  `sd cv-Cindex` = sdcvCindex, 
                  `risk score plot` = plot1, 
                  `signature plot` = plot2)
  
  return(results)
}
