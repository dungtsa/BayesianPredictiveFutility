# BayesianPredictiveFutility
The tool provides futility interim analysis plan using the Bayesian predictive design in single arm early  phase II clinical trial. It also generates  statistical plan so clinicians could easily incorporate it into the clinical trial protocol.
(reference: Application of Bayesian predictive probability for interim analysis in single-arm early phase II trial. Chen et al; submitted).


## Features

* The shiny applictaion provides futility interim analysis plan for the Bayesian predictive design design in single arm early  phase II clinical trial and generates a statistical plan to be easily incorporated  into the clinical trial protocol. 

## Installation

Simply run the following from an R console:

```r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("dungtsa/BayesianPredictiveFutility",force = TRUE)
```

## Getting Started

```r
require("BayesianPredictiveFutility")
Bayesian_predictive_Shinny_app()
```
-------------------------------
Snapshot of shiny app: initial 
![snapshot of shiny app: initial](inst/img/shiny1.png)

-------------------------------
Snapshot of shiny app: output
![snapshot of shiny app: output](inst/img/shiny2.png)

-------------------------------
Snapshot of shiny app: output in Word format (through "download" button)
![snapshot of shiny app: output](inst/img/shiny3.png)

output
![output](Supplementary_Simulation_prior_distribution.pdf)

