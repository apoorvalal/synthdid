# %% ####################################################
#' Difference in differences with entropy balancing weights
#' @param Y the observation matrix.
#' @param N0 the number of control units (N_co in the paper). Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps (T_pre in the paper). Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param demean Remove unit FEs before balancing? TRUE allows for extrapolation
#' @return vector of treatment effects
#' @export
eb_estimate = function(Y, N0, T0, demean = FALSE) {
  N = nrow(Y); T = ncol(Y); N1 = N - N0; T1 = T - T0;
  # pre-treatment outcome matrix
  X0 = Y[1:N0, 1:T0]
  X1 = Y[(N0 + 1):N, 1:T0, drop = FALSE]
  Yvec = colMeans(Y[(N0 + 1):N, , drop = FALSE])
  if (demean) {
		# store unit intercepts
    X0fe = rowMeans(X0)
    X1fe = rowMeans(X1)
    X0 = X0 - X0fe
    X1 = X1 - X1fe
  }
  # fit weights
  eb_wts = as.numeric(ebal::ebal_torch(X0, colMeans(X1))$Weights.ebal)
  intercept = ifelse(demean, mean(X1fe) - weighted.mean(X0fe, eb_wts), 0)
  imputed_Y0 = apply(Y[1:N0, ], 2, \(x) weighted.mean(x, eb_wts))
  event_study = colMeans(Y[(N0 + 1):N, , drop = FALSE]) - imputed_Y0
  # ATT
  estimate = mean(event_study[(T0 + 1):length(event_study)]) - intercept
  if(demean) attr(estimate, "intercept") = intercept
  attr(estimate, "weights") = eb_wts
  attr(estimate, "event_study") = event_study
  return(estimate)
}

