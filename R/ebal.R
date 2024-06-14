# %% ####################################################
#' Difference in differences with entropy balancing weights
#' @param Y the observation matrix.
#' @param N0 the number of control units (N_co in the paper). Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps (T_pre in the paper). Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param demean Remove unit FEs before balancing? TRUE allows for extrapolation
#' @return vector of treatment effects
#' @export
eb_estimate = function(Y, N0, T0, demean = TRUE) {
  N = nrow(Y); T = ncol(Y); N1 = N - N0; T1 = T - T0;
  # pre-treatment outcome matrix
  X0 = Y[1:N0, 1:T0]
  X1 = Y[(N0 + 1):N, 1:T0, drop = FALSE]
  if (demean) {
    # store unit intercepts
    X0fe = rowMeans(X0)
    X1fe = rowMeans(X1)
    X0 = X0 - X0fe
    X1 = X1 - X1fe
  }
  # balance on demeaned data
  eb_wts = as.numeric(ebal::ebal_torch(X0, colMeans(X1))$Weights.ebal)
  Y1 = colMeans(Y[(N0 + 1):N, (T0 + 1):T, drop = FALSE])
  Y0 = apply(Y[1:N0, (T0 + 1):T], 2, \(x) weighted.mean(x, eb_wts))
  # differences in each period
  es_coefs = (colMeans(Y[(N0 + 1):N, , drop = FALSE]) -
                 apply(Y[1:N0, , drop = FALSE], 2, \(x) weighted.mean(x, eb_wts))
              )
  # att only
  att_t = Y1 - Y0 # ATT_t
  estimate = mean(att_t)
  attr(estimate, "weights") = eb_wts
  attr(estimate, "event_study") = es_coefs
  return(estimate)
}
