#' Implement many panel estimators at once. Implements DiD, Synth, Matrix Completion, SDID with and without Intercept and time weights
#' @param df data frame (must be a balanced panel)
#' @param unit_id identifier for units (column name)
#' @param time_id identifier for time periods (column name)
#' @param treatment identifier for binary treatment (column name)
#' @param outcome identifier for outcome (column name)
#' @param mccores number of cores to parallelize inference (whenever possible)
#' @param reps number of reps for bootstrap
#' @param methods character vector of methods to implement. Must be subset of c("DID", "Synthetic Control (SC)", "Synthetic DID (SDID)", "Time Weighted DID", "SDID (No Intercept)", "SC with FEs (DIFP)", "Matrix Completion", "SC (Regularized)", "DIFP (Regularized)")
#' @param infmethod Inference method (must be in c("bootstrap", "jackknife", "placebo"))
#' @return list with matrix with estimates and standard errors and underlying `synthdid_estimate` objects, which can individually be plotted / summarized.
#' @references Arkhangelsky et al (2021), "Synthetic Difference in Differences", AER
#' @export
#' @importFrom MCPanel mcnnm_cv
#' @importFrom mcreplicate mc_replicate
panel_estimate = function(df, unit_id, time_id, treatment, outcome,
    mccores = 8, reps = 200,
    methods = c("DID", "Synthetic Control (SC)", "Synthetic DID (SDID)", "Matrix Completion"),
    infmethod = c("jackknife", "placebo", "bootstrap")
    ) {
  # check
  if (requireNamespace("MCPanel", quietly = TRUE)) {
    .ignore = tryCatch(attachNamespace("MCPanel"), error = function(e) e)
  } else {
    stop("Please install the MCPanel Package for this function.")
  }
  infmethod = match.arg(infmethod)
  # call internal function to reshape data to N X T matrix with treated obs at the bottom
  setup = panelMatrices(df, unit_id, time_id, treatment, outcome)
  ######################################################################
  # nonstandard estimators : MC is self-contained, others involve existing functions
  # (feed in different weighting schemes)
  ######################################################################
  mc_imputation = function(Y, N0, T0, N1, T1){
    W = outer(c(rep(0, N0), rep(1, N1)), c(rep(0, T0), rep(1, T1)))
    mc_pred = MCPanel::mcnnm_cv(Y, 1 - W, num_lam_L = 20)
    mc_fit = mc_pred$L + outer(mc_pred$u, mc_pred$v, '+')
  }
  mc_estimate = function(Y, N0, T0) {
    N1 = nrow(Y) - N0
    T1 = ncol(Y) - T0
    W = outer(c(rep(0, N0), rep(1, N1)), c(rep(0, T0), rep(1, T1)))
    mc_fit = mc_imputation(Y, N0, T0, N1, T1)
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
  # special cases of existing fns
  difp_estimate = \(Y, N0, T0) {
    synthdid_estimate(Y, N0, T0, weights = list(lambda = rep(1 / T0, T0)), eta.omega = 1e-6)
  }
  difp_estimate_reg = \(Y, N0, T0) {
    synthdid_estimate(Y, N0, T0, weights = list(lambda = rep(1 / T0, T0)))
  }
  twdid_estimate_reg = \(Y, N0, T0) {
    synthdid_estimate(Y, N0, T0, weights = list(omega = rep(1 / N0, N0)))
  }
  sc_estimate_reg = \(Y, N0, T0) {
    sc_estimate(Y, N0, T0, eta.omega = ((nrow(Y) - N0) * (ncol(Y) - T0))^(1 / 4))
  }
  sdid_estimate_noint = \(Y, N0, T0) {
    synthdid_estimate(Y, N0, T0, omega.intercept = FALSE)
  }
  # list of estimation functions
  all_estimators = list(
    "DID" = did_estimate,
    "Synthetic Control (SC)" = sc_estimate,
    "Synthetic DID (SDID)" = synthdid_estimate,
    "Time Weighted DID"   = twdid_estimate_reg,
    "SDID (No Intercept)" = sdid_estimate_noint,
    "SC with FEs (DIFP)" = difp_estimate,
    "Matrix Completion" = mc_estimate,
    "SC (Regularized)" = sc_estimate_reg,
    "DIFP (Regularized)" = difp_estimate_reg
  )
  # select subset
  selected = which(!is.na(match(names(all_estimators), methods)))
  estimators = all_estimators[selected]
  # apply estimator functions to setup matrix
  estimates = lapply(estimators, \(e) e(setup$Y, setup$N0, setup$T0))
  # inference
  standard_errors = mapply(\(e, n) {
    switch(n,
      "Matrix Completion" = mc_placebo_se(setup$Y, setup$N0, setup$T0),
      sqrt(vcov(e, method = infmethod, ncores = mccores, replications = reps))
    )
  }, estimates, names(estimators))
  # return table as a matrix
  sumtab = rbind(unlist(estimates), unlist(standard_errors))
  rownames(sumtab) = c("Est", "Std. Error")
  retobj = list(
      ests = estimates,
      summary_table = sumtab
    )
  class(retobj) = "panelEstimate"
  return(retobj)
}

# %% ####################################################
#' print method (prints summary table)
#' @param data [data.frame, matrix] data table (n x p)
#' @export
print.panelEstimate = \(x) {
  print(x$summary_table)
}

# construct event study estimates from class synthdid_estimate
event_study_coefs = \(est) {
  setup = attr(est, 'setup')
  weights = attr(est, 'weights')
  Y = setup$Y
  N0 = setup$N0; N1 = nrow(Y) - N0
  T0 = setup$T0; T1 = ncol(Y) - T0
  lambda.synth = c(weights$lambda, rep(0, T1))
  lambda.target = c(rep(0, T0), rep(1 / T1, T1))
  omega.synth = c(weights$omega, rep(0, N1))
  omega.target = c(rep(0, N0), rep(1 / N1, N1))
  offset = c((omega.target - omega.synth) %*% Y %*% lambda.synth)
  obs.trajectory = as.numeric(omega.target %*% Y)
  syn.trajectory = as.numeric(omega.synth %*% Y) + offset
  obs.trajectory - syn.trajectory
}

# evstudy_fig = \(ests, title, yl = 3) {
#   # compute ATTs (added to legend)
#   atts = round(ests[8,], 3)
#   colours = RColorBrewer::brewer.pal(ncol(ests), "Paired")
#   matplot(ests,
#     type = 'l', col = colours, lwd = 2, lty = 2,
#     ylab = "Coef. Estimate", xlab = "Time", axes = F,
#     ylim = c(-1*yl, yl),
#     main = title
#   )
#   matpoints(ests, pch = 16, col = colours, lwd = 2, lty = 2)
#   abline(h = 0, lty = 3)
#   abline(v = 7.5, lty = 4)
#   legend("bottomleft", colnames(ests),
#     col = colours, lty = 1, ncol = 2, cex = 0.9,
#     x.intersp = 0.8, text.width = 2, lwd = 5
#   )
#   axis(2)
#   axis(side = 1, at = 1:nrow(ests), labels = seq(1992, 2020, 4))
# }

