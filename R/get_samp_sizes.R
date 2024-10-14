#' Calculate sample sizes by subgroup
#'
#' @param Y a vector of responses
#' @param X matrix of covariates
#' @param grouping_col name of a column in X that defines subgroups
#' @param preds vector of predictions
#' @param effect_size desired minimal detectable difference
#' @param alpha desired Type II error rate
#' @param beta desired Type II error rate
#'
#' @return A data frame with 7 columns: `def`, `num`, `MSE`, `variance`,`resid_var`,
#' `samp_size`, and `samp_size_without`.
#' `def` is the definition of the subgroup, which may be an interval or a single
#' factor level, depending on the values in `grouping_col`. `num` is the number of
#' observations in that subgroup. `MSE` is the mean squared error of the predictions
#' for that subgroup. `variance` is the variance Y of the subgroup. `resid_var` is
#' the variance of the residuals calculated from `preds` and `Y`. `resid_var` is assumed
#' to be \eqn{\sigma^2} and is passed through the `samp_size` function. `samp_size`
#' is the resulting estimated sample size needed if the RCT population resembles this
#' subgroup. `samp_size_without` is the estimated sample size needed if the RCT
#' population resembles this subgroup and auxiliary data is not used to improve precision.
#' @export
#'
get_samp_sizes <- function(Y, X, grouping_col, preds, effect_size, alpha, beta) {

  #effect size needs to be a vector passed through based on if it's standard deviation based or not
  #but if it's not, it should just be that value repeated??
  #also why are we always getting 785 for without???
  #or, pass through a method because otherwise we don't have access to variance
  #effect_size_units = "sd" or effect_size_units = "raw" and then create a vector from this and pass it through?
  #first figure out 785 issue
  #785 issue is because the sigma squared estimate and the effect size are cancelling out??

  if (!grouping_col %in% colnames(X)) {
    stop("grouping_col must be the name of a column in X")
  }


  df <- cbind(Y, X, preds)

  groups <- unique(X[[grouping_col]])

  results <- data.frame()
  indices <- list()
  effect_size_vec <- rep(NA, length(groups))
  for (i in 1:length(groups)) {
    subgroup_index <- which(df[[grouping_col]] == groups[i])
    indices[[i]] <- subgroup_index
    subgroup <- df[subgroup_index, ]
    effect_size_vec[i] <- effect_size * sd(subgroup$Y)
    r <- data.frame(def = groups[i], 'MSE' = sum((subgroup$preds - subgroup$Y)^2)/nrow(subgroup), 'num' = nrow(subgroup),
                    'variance' = var(subgroup$Y), "resid_var" = var(subgroup$preds-subgroup$Y))
    results <- rbind(results, r)
  }


  results$samp_size <- unlist(lapply(as.numeric(results$resid_var),samp_size, effect_size = effect_size, alpha = alpha, beta = beta))
  results$samp_size_without <- unlist(lapply(as.numeric(results$variance), samp_size, effect_size = effect_size, alpha = alpha, beta = beta))

  # results$samp_size <- unlist(purrr::map2(as.numeric(results$resid_var), effect_size_vec, samp_size, alpha = alpha, beta = beta))
  # results$samp_size_without <- unlist(purrr::map2(as.numeric(results$variance), effect_size_vec, samp_size, alpha = alpha, beta = beta))

  #the sorting is going to be tricky here, we have to resort indices in the same way?
  sort_order <- order(results$def)

  #got to be a more efficient way to do this? TO - DO
  indices_new <- list()
  for (j in 1:length(indices)) {
    indices_new[[j]] <- indices[[sort_order[j]]]
  }

  return(list(arrange(results, def), indices_new))
}
