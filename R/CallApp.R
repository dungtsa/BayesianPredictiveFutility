#' @name Bayesian_Predictive_App
#' @title Bayesian predictive design in Shinny application
#' @description This function will run Shiny application for the Bayesian predictive design in single arm early  phase II clinical trial.  The function will automatically generate a statistical plan.
#' @export

Bayesian_Predictive_App <- function() {
  appDir <- system.file("shiny-examples", "app", package = "BayesianPredictiveFutility")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `BayesianPredictiveFutility`.", call. = FALSE)
  }

  shiny::runApp(paste(appDir,'/app.R',sep = ''), launch.browser = T,host = getOption( "127.0.0.1"))
}


