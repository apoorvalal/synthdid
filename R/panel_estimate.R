#' Implement many panel estimators at once. Implements DiD, Synth, Matrix Completion, SDID with and without Intercept and time weights
#' @param df data frame (must be a balanced panel)
#' @param unit_id identifier for units (column name)
#' @param time_id identifier for time periods (column name)
#' @param treatment identifier for binary treatment (column name)
#' @param outcome identifier for outcome (column name)
#' @param mccores number of cores to parallelize inference (whenever possible)
#' @param reps number of reps for bootstrap
#' @param infmethod Inference method (must be in c("bootstrap", "jackknife", "placebo"))
#' @return matrix with estimates and standard errors
#' @references Arkhangelsky et al (2021), "Synthetic Difference in Differences", AER
#' @export
#' @importFrom MCPanel mcnnm_cv
panel_estimate = function(df, unit_id, time_id, treatment, outcome,
    mccores = 8, reps = 200,
    infmethod = c("jackknife", "placebo", "bootstrap")
    ) {
  infmethod = match.arg(infmethod)
  # call internal function to reshape data to N X T matrix with treated obs at the bottom
  setup = panelMatrices(df, unit_id, time_id, treatment, outcome)
  ######################################################################
  # nonstandard estimators : MC is self-contained, others involve existing functions
  # (feed in different weighting schemes)
  ######################################################################
  mc_estimate = function(Y, N0, T0) {
    N1 = nrow(Y) - N0
    T1 = ncol(Y) - T0
    W = outer(c(rep(0, N0), rep(1, N1)), c(rep(0, T0), rep(1, T1)))
    mc_pred = MCPanel::mcnnm_cv(Y, 1 - W, num_lam_L = 20)
    mc_fit = mc_pred$L + outer(mc_pred$u, mc_pred$v, '+')
    mc_est = sum(W * (Y - mc_fit)) / sum(W)
    mc_est
  }
  # standard error constructor for matrix comletion
  mc_placebo_se = function(Y, N0, T0) {
    N1 = nrow(Y) - N0
    theta = \(ind)  mc_estimate(Y[ind, ], length(ind) - N1, T0)
    sqrt((reps - 1) / reps) * sd(
      mcreplicate::mc_replicate(reps, theta(sample(1:N0)), mc.cores = mccores),
    )
  }
  difp_estimate = \(Y, N0, T0){
    synthdid_estimate(Y, N0, T0, weights = list(lambda = rep(1 / T0, T0)), eta.omega = 1e-6)
  }
  difp_estimate_reg = \(Y, N0, T0){
    synthdid_estimate(Y, N0, T0, weights = list(lambda = rep(1 / T0, T0)))
  }
  sc_estimate_reg = \(Y, N0, T0){
    sc_estimate(Y, N0, T0, eta.omega = ((nrow(Y) - N0) * (ncol(Y) - T0))^(1 / 4))
  }
  sdid_estimate_noint = \(Y, N0, T0){
    synthdid_estimate(Y, N0, T0, omega.intercept = FALSE)
  }
  # list of estimation functions
  estimators = list(
    "DID" = did_estimate,
    "Synthetic Control (SC)" = sc_estimate,
    "Synthetic DID (SDID)" = synthdid_estimate,
    "SDID (No Intercept)" = sdid_estimate_noint,
    "SC with FEs (DIFP)" = difp_estimate,
    "Matrix Completion" = mc_estimate,
    "SC (Regularized)" = sc_estimate_reg,
    "DIFP (Regularized)" = difp_estimate_reg
  )
  # apply estimator functions to setup matrix
  estimates = lapply(estimators, \(e) e(setup$Y, setup$N0, setup$T0))
  # inference
  standard_errors = mapply(function(e, n) {
    switch(n,
      "Matrix Completion" = mc_placebo_se(setup$Y, setup$N0, setup$T0),
      sqrt(vcov(e, method = infmethod))
    )
  }, estimates, names(estimators))
  # return table as a matrix
  rbind(unlist(estimates), unlist(standard_errors))
}
