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
library("BayesianPredictiveFutility")
Bayesian_Predictive_App()
```


-------------------------------
Procedure of BayesianPredictiveFutility R Shiny App  
![snapshot of shiny app: initial](inst/img/ProcedureOfBayesianPredictiveFutilityRShinyApp.png)

-------------------------------
Example of two-stage case
![Example of two-stage case](Example/Demonstration_Two_Stage.pdf)

-------------------------------
Example of three-stage case
![Example of three-stage case](Example/Demonstration_Three_Stage.pdf)

-------------------------------
Example of multi-stage case
![Example of multi-stage case](Example/Demonstration_Multi_Stage.pdf)
