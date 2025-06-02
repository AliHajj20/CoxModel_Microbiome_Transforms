```{r}
#risk score
# Ensure coef_vec is numeric and named
coef_vec <- as.numeric(coef(coxmodel_logratio))
names(coef_vec) <- names(coef(coxmodel_logratio))

# Make sure x_test is a numeric matrix
x_test <- as.matrix(x_test)

# Align x_test columns to match the coefficient names
x_test <- x_test[, names(coef_vec), drop = FALSE]

# Compute the risk scores
risk_score_test <- as.vector(x_test %*% coef_vec)
#Harrell's C-index
library(survival)
concordance <- survConcordance(Surv(time_test, status_test) ~ risk_score_test)
print(concordance$concordance)





#code to calculate convert the assay to all pairwise log ratio
logratios <- function(x, type = "all") {
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- ncol(x)
  combs <- combn(n, 2)
  lr_mat <- apply(combs, 2, function(idx) {
    log(x[, idx[1]] / x[, idx[2]])
  })
  colnames(lr_mat) <- apply(combs, 2, function(idx) {
    paste0("log(", colnames(x)[idx[1]], "/", colnames(x)[idx[2]], ")")
  })
  return(as.data.frame(lr_mat))
}
