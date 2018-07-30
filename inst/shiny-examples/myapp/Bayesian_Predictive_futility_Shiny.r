library(shiny)
library(knitr)
library(rmarkdown)
library(extraDistr)

ui <- fluidPage(
  headerPanel('Interim Analysis for Futility Based on Bayesian Predictive Probability'),
  sidebarLayout(
    sidebarPanel(
      textInput('nlist', 'Enter sample size for each interim analysis (comma delimited; e.g., 25,25)', "25,25"),
      numericInput("pTarget",label = "Probability of  unfavorable (p0)",min = 0, max = 1, value = 0.3),
      numericInput("pH1",label = "Probability of  favorable (p1: H1 with p1>p0)",min = 0, max = 1, value = 0.5),
      numericInput("theta",label = "Cutoff for greater p0 (prob(x>p0|data)>theta",min = 0, max = 1,value = 0.95),
      numericInput("cutoffPredictive",label = "Cutoff for predictive probability",min = 0, max = 1,value = 0.05),
      numericInput("betaA",label = "beta A parameter",min = 0, max = 1,value = 1),
      numericInput("betaB",label = "beta B parameter",min = 0, max = 1,value = 1),
      numericInput("simN",label = "number of simulation",value = 10000),
      textInput('outcomeLabel', label = 'Name of outcome', value = "response"),
      textInput('armLabel', label = 'Name of the Arm', value = "treatment"),
      actionButton("Submit","Calculate"),
      radioButtons('format', 'Document format', c('Word','PDF', 'HTML' ),
                   inline = TRUE),
      downloadButton('downloadReport')
    ),
    
    mainPanel(
      h2("Summary of Interim Analysis for Futility"),
      textOutput("summary"),
      tags$hr(),
      
      h3("Details"),
      textOutput("futility"),
      tags$hr(),
      
      h4("Table 1: Stopping Boundary for Futility"),
      tableOutput("stoppingBoundary"),
      tags$hr(),
      
      
      h4("Table 2: Bayesian Predictive Probability"),
      uiOutput("predictiveData"),
      tags$hr(),
      
      h4("Table 3: Sensitivity Analysis"),
      tableOutput("sensitivityData"),
      tags$hr(),
      
      # This is the dynamic UI for the plots
      h4("Figure 1: Bayesian Predictive Probability"),
      uiOutput("plotPredictive"),
      
      h4("Figure 2: Sensitivity"),
      plotOutput("plotSensitivityOverall"),
      
      h4("Figure 3: Probability of Stopping Early by Each Interim Analysis"),
      uiOutput("plotSensitivityInd"),
      
      tags$hr()
    )
    
  )
)


server <- function(input,output){
  
  
  get.result.futility <- eventReactive(input$Submit, {
    n.list <- as.numeric(unlist(strsplit(input$nlist,",")))
    pTarget <- input$pTarget
    pH1 <- input$pH1
    if (pH1 > pTarget) 
      increase <- TRUE else 
        increase <- FALSE
    theta <- input$theta
    cutoffPredictive <- input$cutoffPredictive
    outcomeLabel <- input$outcomeLabel
    armLabel <- input$armLabel
    betaA <- input$betaA
    betaB <- input$betaB
    simN <- input$simN
    if (increase) {
      # Increase Selected
      tmp99 <- Bayesian.predictive.futility.increase.fun(n.list = n.list,
                                                         p.target = pTarget,
                                                         p.h1 = pH1,
                                                         cutoff.n.for.greater.p0 = theta,
                                                         predictive.cutoff = cutoffPredictive,
                                                         beta.a = betaA,
                                                         beta.b = betaB,
                                                         sim.n = simN,
                                                         outcome.tmp = outcomeLabel,
                                                         arm.name = armLabel) 
      
    } else {
      # Decrease Selected
      tmp99 <- Bayesian.predictive.futility.reduction.fun(n.list = n.list,
                                                          p.target = pTarget,
                                                          p.h1 = pH1,
                                                          cutoff.n.for.greater.p0 = theta,
                                                          predictive.cutoff = cutoffPredictive,
                                                          beta.a = betaA,
                                                          beta.b = betaB,
                                                          sim.n = simN,
                                                          outcome.tmp = outcomeLabel,
                                                          arm.name = armLabel) 
      
    }
    
    n.needed.for.greater.p0 <- tmp99$n.needed.for.greater.p0
    n.predictive.cutoff <- tmp99$n.predictive.cutoff
    all.data <- tmp99$all.data
    predictive.prob.data <- tmp99$predictive.prob.data
    sensitivity.data <- tmp99$sensitivity.data
    sensitivity.data.ind <- tmp99$sensitivity.data.ind
    kk2.len.list <- tmp99$kk2.len.list
    n.additional.needed <- n.needed.for.greater.p0 - n.predictive.cutoff
    
    len.n.list <- length(n.list) - 1
    stage.index <- ifelse(len.n.list == 1,' 1st ',ifelse(len.n.list == 2,' 2nd ',ifelse(len.n.list == 3,' 3rd ',paste(' ',len.n.list,'th ',sep = ''))))
    
    
    a1 <- paste0('A Bayesian approach for futility analysis is used to calculate posterior probability and predictive probability for the rate of ',outcomeLabel,' with a non-informative beta prior, beta(',betaA,',',betaB,'). We consider a ',round(pTarget*100),'% rate or ',ifelse(increase, 'lower', 'higher'),' of ',outcomeLabel, ' as ineffective for the treatment. Thus, we expect the ',armLabel,' arm is not ',ifelse(increase, 'worse', 'better'),' than the historical control if the posterior probability of the rate (',outcomeLabel,') ',ifelse(increase, 'greater', 'less'),' than ',round(pTarget*100),'% is higher than ',theta,' (i.e., prob(rate of ',outcomeLabel,ifelse(increase, '>', '<'),round(pTarget*100),'% |data)>',theta,') ). \n \n With a total ',sum(n.list),' patients in ',armLabel,' arm, the number of patients with ',outcomeLabel,' needs to be ',n.needed.for.greater.p0,' or ',ifelse(increase, 'more', 'less'),' in order to meet the criteria. Therefore, we use the number of ',n.needed.for.greater.p0,' patients to guide the predictive probability. Specifically, for the first ',n.list[1],' patients in the 1st interim analysis, there are ',n.list[1] + 1,' ways for number of patients with ',outcomeLabel,' from 0, 1, to, ',n.list[1],'. In each case, given the number of patients with ',outcomeLabel,', $s$, in the first ',n.list[1],' patients, we calculate predictive probability of $',n.needed.for.greater.p0,'-s$ or ',ifelse(increase, 'more', 'less'),' patients with ',outcomeLabel,' in the future remaining ',sum(n.list) - n.list[1])
    nn1 <- n.list[1]
    nn2 <- sum(n.list) - n.list[1]
    n.predictive.cutoff.text <- n.predictive.cutoff
    n.predictive.cutoff.text[n.predictive.cutoff.text == (-99)] <- 'Pass'
    
    index.example <- n.predictive.cutoff == (-99)
    nn1.example <- sum((n.list[-length(n.list)])[index.example]) + (n.list[-length(n.list)])[!index.example][1]
    nn2.example <- sum(n.list) - nn1.example
    n.predictive.cutoff.example <- n.predictive.cutoff[!index.example][1]
    
    a1 <- paste0(a1,' patients, i.e., $\\sum_{i=',n.needed.for.greater.p0,'-s}^{',nn2,'} {',nn2,' \\choose i } \\frac{beta(',betaA,'+s+i,',betaB,'+(',nn1,'-s)+(',nn2,'-i))}{beta(',betaA,'+s,',betaB,'+(',nn1,'-s))}$.
              Calculation of predictive probability is based on beta binominal distribution for the number of patients with ',outcomeLabel,' in the future remaining ',nn2,' patients given a beta distribution for the  rate of ',outcomeLabel,', $beta(',betaA,'+s,',betaB,'+',nn1,'-s)$. 
              For example, if there are ',n.predictive.cutoff.example[1],' patients with ',outcomeLabel,' in the first ',nn1.example,' patients, the predictive probability of ',n.needed.for.greater.p0 - n.predictive.cutoff.example[1],' or ',ifelse(increase, 'more', 'less'),' patients with ',outcomeLabel,' in the future remaining ',nn2.example,' patients would be 
              $\\sum_{i=',n.needed.for.greater.p0 - n.predictive.cutoff.example[1],'}^{',nn2.example,'}{',nn2.example,' \\choose i } \\frac{beta(',betaA,'+',n.predictive.cutoff.example[1],'+i,',betaB,'+(',nn1,'-',n.predictive.cutoff.example[1],')+(',nn2.example,'-i))}{beta(',betaA,'+',n.predictive.cutoff.example[1],',',betaB,'+(',nn1.example,'-',n.predictive.cutoff.example[1],'))} =',round(predictive.prob.data[[1]][predictive.prob.data[[1]][,3] == n.predictive.cutoff.example[1],5],3),'$.
              \n \n The predictive probability is also calculated for each of the remaining interim analyses to evaluate the chance of ',n.needed.for.greater.p0,'-s or ',ifelse(increase, 'more', 'less'),' patients with ',outcomeLabel,' in the future remaining patients given s patients with ',outcomeLabel,' in the current stage of interim analysis. ',
                ' Figure 1 and Table 2 lists predictive probability for all scenarios of number of patients with ',outcomeLabel,' in each interim analysis and the associated largest number of patients with ',outcomeLabel,' needed in the future remaining patients to have at ',ifelse(increase, 'least', 'most'),' a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,'. ', 
                ' \n \n We consider that a ',round(cutoffPredictive*100),'% cutoff of the predictive probability will give little chance to have ',n.needed.for.greater.p0,' patients or ',ifelse(increase, 'more', 'less'),' with ',outcomeLabel,' at the end of study. Thus with this cutoff, the stopping rule (Table 1) will be: the trial will be stopped if there are ',ifelse(length(n.list) == 2,paste(n.predictive.cutoff,' or ',ifelse(increase, 'less', 'more'),' patients with ',outcomeLabel,' in the 1st interim analysis. ',sep = ''),paste(paste(paste(n.predictive.cutoff.text[-length(n.predictive.cutoff)],collapse = ', '),n.predictive.cutoff.text[length(n.predictive.cutoff)],sep = ', and '),' or ',ifelse(increase, 'less', 'more'),' patients with ',outcomeLabel,' in the 1st to ',stage.index,' interim analysis, respectively. ',sep = '')),
                ' Sensitivity analysis (Figure 2 and Table 3) for this stopping rule shows that if the true rate of ',outcomeLabel,' is ',round(pTarget*100),'%, the chance to reach at ',ifelse(increase, 'least', 'most'),' a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,2] * 100),'% (Type I error), however the probability to stop the trial early is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,1],2)*100 ,'%. \n \n When the true rate of ',outcomeLabel,' is ',as.numeric(pH1)*100,'%, then the chance to reach at ',ifelse(increase, 'least', 'most'),' ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[paste(pH1),2] * 100),'%, and the corresponding probability to stop the treatment early is ',round(sensitivity.data[paste(pH1),1],2)*100,'%.',ifelse(ncol(sensitivity.data.ind) > 3,' \n \n Figure 3 shows the probability of stopping at each interim analysis.',''))
    
    
    a1_summary <- paste0('Futility evaluation is implemented in ',length(n.predictive.cutoff),ifelse(length(n.predictive.cutoff) == 1,' interim analysis',' interim analyses'),', with stops at the first ',paste(cumsum(n.list[-length(n.list)]),collapse = ', '),' patients, with ', n.list[length(n.list)],' patients in the last stage, for a total of ',sum(n.list),' patients . With an unfavorable rate set at ',round(pTarget*100),'% and posterior probability of ', theta,', a total of at least ', n.needed.for.greater.p0, ' of the ', sum(n.list),' patients must have ', outcomeLabel, ' to be able to detect a difference. Given a ',round(cutoffPredictive*100),'% cutoff of the predictive probability (i.e. chance needed to ultimately see a difference, given the results from the ', ifelse(length(n.predictive.cutoff) == 1,' interim analysis',' interim analyses'), ') the stopping rule (Table 1) will be: the trial will be stopped if there are ',ifelse(length(n.list) == 2,paste(n.predictive.cutoff,' or ',ifelse(increase, 'less', 'more'),' patients with ',outcomeLabel,' in the 1st interim analysis. ',sep = ''),paste(paste(paste(n.predictive.cutoff.text[-length(n.predictive.cutoff)],collapse = ', '),n.predictive.cutoff.text[length(n.predictive.cutoff)],sep = ', and '),' or ',ifelse(increase, 'less', 'more'),' patients with ',outcomeLabel,' in the 1st to ',stage.index,' interim analysis, respectively. ',sep = '')), ' \n \n Sensitivity analysis (Figure 2 and Table 3) shows that if the true rate of ',outcomeLabel,' is ',round(pTarget*100),'%, the chance to reach at ',ifelse(increase, 'least', 'most'),' a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,2] * 100),'% (Type I error), however the probability to stop the trial early is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,1],2)*100 ,'%. \n \n If true rate of ',outcomeLabel,' is indeed ',as.numeric(pH1)*100,'%, then the chance to reach at ',ifelse(increase, 'least', 'most'),' ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[paste(pH1),2] * 100),'%, and the corresponding probability to stop the treatment early is ',round(sensitivity.data[paste(pH1),1],2)*100,'%.')
                         
                         
                         
    a2 <- paste('Application of Bayesian predictive probability for interim analysis in single-arm early phase II trial. Chen et al; submitted.')
    
    tmp98 <- list(a1_summary = a1_summary,
                  a1 = a1,
                  a2 = a2,
                  n.needed.for.greater.p0 = n.needed.for.greater.p0,
                  n.predictive.cutoff = n.predictive.cutoff,
                  all.data = all.data,
                  predictive.prob.data = predictive.prob.data,
                  sensitivity.data = unique(round(sensitivity.data,4)),
                  sensitivity.data.ind = unique(round(sensitivity.data.ind,4)),
                  n.list = n.list,
                  pH1 = pH1,
                  pTarget = pTarget,
                  increase = increase,
                  theta = theta,
                  cutoffPredictive = cutoffPredictive,
                  outcomeLabel = outcomeLabel,
                  armLabel = armLabel,
                  betaA = betaA,
                  betaB = betaB,
                  kk2.len.list = kk2.len.list
    )
    return(tmp98)
  })
  
  
  output$summary <- renderText({
    get.result.futility()$a1_summary
  })
  
  output$futility <- renderText({
    get.result.futility()$a1
  })
  
  output$stoppingBoundary <- renderTable({
    aa0 <- get.result.futility()
    aa <- aa0$n.list
    cutoff.tmp <- c(aa0$n.predictive.cutoff,aa0$n.needed.for.greater.p0 - 1)
    cutoff.tmp[cutoff.tmp == (-99)] <- 'PASS'
    aa <- rbind(c(1:(length(aa) - 1),'Final'),cumsum(aa),aa,cutoff.tmp)
    dimnames(aa)[[1]] <- c('Stage of interim analysis','Sample size up to the current stage','Sample size at each stage','Stopping boundary')
    aa
  },rownames = T,colnames = F, digits = 0)
  
  
  output$predictiveData <- renderUI({
    n.list <- get.result.futility()$n.list
    table_output_list <- lapply(1:(length(n.list) - 1), function(i) {
      tablename <- paste("PredictiveTable", i, sep = "")
      tableOutput(tablename)
    })
    do.call(tagList, table_output_list)
  })
  
  observe({
    aa0 <- get.result.futility()
    aa.list <- aa0$predictive.prob.data
    
    for (j in 1:length(aa.list))
    {
      local({
        my_j <- j
        tablename <- paste("PredictiveTable", my_j, sep = "")
        stage.index <- ifelse(my_j == 1,' 1st ',ifelse(my_j == 2,' 2nd ',ifelse(my_j == 3,' 3rd ',paste(' ',my_j,'th ',sep = ''))))
        aa <- aa.list[[my_j]][,3:5,drop = F]
        outcomeLabel <- aa0$outcomeLabel
        dimnames(aa)[[2]] <- c(paste('number of patients with ', outcomeLabel,' in the ',stage.index,' interim analysis',sep = ''),paste(ifelse(aa0$increase, 'minimum', 'maximum'),' number of patients with ',outcomeLabel,' ',ifelse(aa0$increase, 'needed', 'allowed'),' in the future remaining patients',sep = ''),'predictive probability')
        output[[tablename]] <- renderTable({
          aa
        })
      })
    }
  })
  
  
  output$sensitivityData <- renderTable({
    aa0 <- get.result.futility()
    aa <- aa0$sensitivity.data
    aa <- cbind(as.numeric(dimnames(aa)[[1]]),aa)
    dimnames(aa)[[2]] <- c('true rate',' overall probability of early stopping the trial',paste('probability to have at ',ifelse(aa0$increase, 'least', 'most'),' ',aa0$n.needed.for.greater.p0,' patients with ',aa0$outcomeLabel,sep = ''))
    aa
    
  })
  
  output$plotPredictive <- renderPlot({
    res <-  get.result.futility()
    n.list <- res$n.list
    pTarget <- res$pTarget
    theta <- res$theta
    predictive.cutoff <- cutoffPredictive <- res$cutoffPredictive
    outcomeLabel <- res$outcomeLabel
    armLabel <- res$armLabel
    betaA <- res$betaA
    betaB <- res$betaB
    kk2.len.list <- res$kk2.len.list
    
    n.predictive.cutoff.list <- res$n.predictive.cutoff
    
    tmp3.list <- res$predictive.prob.data
    
    
    par(mfrow = c(2,ceiling(length(tmp3.list)/2)))
    for (j in 1:length(tmp3.list)) 
    {
      tmp3 <- tmp3.list[[j]]
      n1 <- sum(n.list[1:j])
      n2 <- sum(n.list[(j + 1):length(n.list)])
      kk2.len <- kk2.len.list[[j]]
      n.predictive.cutoff <- n.predictive.cutoff.list[j]
      stage.index <- ifelse(j == 1,' 1st ',ifelse(j == 2,' 2nd ',ifelse(j == 3,' 3rd ',paste(' ',j,'th ',sep = ''))))
      plot(tmp3[,5],type = 'h',axes = F,xlab = paste('number of patients with ', outcomeLabel,' in the',stage.index,' interim analysis with ',n1,' patients',sep = ''),ylab = 'predictive probability',lwd = 3,cex.lab = 1.5,ylim = c(0,1))
      box();axis(2,cex.axis = 1.2)
      axis(2,predictive.cutoff,cex.axis = 1.2)
      axis(1,1:dim(tmp3)[1],tmp3[,3],cex.axis = 1.2)
      axis(3,1:dim(tmp3)[1],tmp3[,4],cex.axis = 1)
      abline(h = predictive.cutoff,col = 2,lty = 2)
      mtext(paste(ifelse(res$increase, 'minimum', 'maximum'),' number of patients with ',outcomeLabel,' ',ifelse(res$increase, 'needed', 'allowed'),' in the future remaining ', n2,' patients',sep = ''),side = 3,line = 1.8)
      title(paste(stage.index,' Interim Analysis of for Futility \n\n',sep = ''))
      arrows(kk2.len, 1, kk2.len, y1 = 0.1,col = 3,lwd = 3)
      text(1:dim(tmp3)[1],tmp3[,5],round(tmp3[,5],2),col = 2,cex = 2)
    }
    
  })
  
  
  # Insert the right number of plot output objects into the web page
  output$plotPredictive <- renderUI({
    n.list <- get.result.futility()$n.list
    plot_output_list <- lapply(1:(length(n.list) - 1), function(j) {
      plotname <- paste("plot", j, sep = "")
      plotOutput(plotname)
      #plotOutput(plotname, height = 280, width = 250)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  observe({
    res <-  get.result.futility()
    n.list <- res$n.list
    pTarget <- res$pTarget
    theta <- res$theta
    predictive.cutoff <- cutoffPredictive <- res$cutoffPredictive
    outcomeLabel <- res$outcomeLabel
    armLabel <- res$armLabel
    betaA <- res$betaA
    betaB <- res$betaB
    kk2.len.list <- res$kk2.len.list
    
    n.predictive.cutoff.list <- res$n.predictive.cutoff
    
    tmp3.list <- res$predictive.prob.data
    
    for (j in 1:length(tmp3.list)) {
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderPlot() will be the same across all instances, because
      # of when the expression is evaluated.
      local({
        my_j <- j
        plotname <- paste("plot", my_j, sep = "")
        tmp3 <- tmp3.list[[my_j]]
        n1 <- sum(n.list[1:my_j])
        n2 <- sum(n.list[(my_j + 1):length(n.list)])
        kk2.len <- kk2.len.list[[my_j]]
        n.predictive.cutoff <- n.predictive.cutoff.list[my_j]
        stage.index <- ifelse(my_j == 1,' 1st ',ifelse(my_j == 2,' 2nd ',ifelse(my_j == 3,' 3rd ',paste(' ',my_j,'th ',sep = ''))))
        
        output[[plotname]] <- renderPlot({
          plot(tmp3[,5],type = 'h',axes = F,xlab = paste('number of patients with ', outcomeLabel,' in the',stage.index,' interim analysis with ',n1,' patients',sep = ''),ylab = 'predictive probability',lwd = 3,cex.lab = 1.5,ylim = c(0,1))
          box();axis(2,cex.axis = 1.2)
          axis(2,predictive.cutoff,cex.axis = 1.2)
          axis(1,1:dim(tmp3)[1],tmp3[,3],cex.axis = 1.2)
          axis(3,1:dim(tmp3)[1],tmp3[,4],cex.axis = 1)
          abline(h = predictive.cutoff,col = 2,lty = 2)
          mtext(paste(ifelse(res$increase, 'minimum', 'maximum'),' number of patients with ',outcomeLabel,' ',ifelse(res$increase, 'needed', 'allowed'),' in the future remaining ', n2,' patients',sep = ''),side = 3,line = 1.8)
          title(paste(stage.index,' Interim Analysis of for Futility \n\n',sep = ''))
          arrows(kk2.len, 1, kk2.len, y1 = 0.1,col = 3,lwd = 3)
          text(1:dim(tmp3)[1],tmp3[,5],round(tmp3[,5],2),col = 2,cex = 2)
        })
      })
    }
  })  
  
  output$plotSensitivityOverall <- renderPlot({
    res <-  get.result.futility()
    tmp5 <- res$sensitivity.data
    tmp5.ind <- res$sensitivity.data.ind
    n.list <- res$n.list
    pTarget <- res$pTarget
    theta <- res$theta
    cutoffPredictive <- res$cutoffPredictive
    outcomeLabel <- res$outcomeLabel
    armLabel <- res$armLabel
    betaA <- res$betaA
    betaB <- res$betaB
    n.predictive.cutoff <- res$n.predictive.cutoff
    n.needed.for.greater.p0 <- res$n.needed.for.greater.p0
    
    par(mfrow = c(2,1))
    aa1 <- barplot(tmp5[,1],xlab = paste('Rate of ',outcomeLabel,sep = ''),ylab = 'Probability of early stopping the trial overall',main = 'Sensitivity Analysis \n Overall probability of early stopping the trial',lwd = 3,cex.lab = 1.5,ylim = c(0,1),cex.main = 2,cex.lab = 1.5,cex.axis = 1.5,cex = 1.5)
    text(aa1,tmp5[,1],round(tmp5[,1],2),col = 2,cex = 2)
    
    aa1 <- barplot(tmp5[,2],xlab = paste('Rate of ',outcomeLabel,sep = ''),ylab = paste('Probability of at ',ifelse(res$increase, 'least', 'most'),' ', n.needed.for.greater.p0,' patients with ', outcomeLabel,sep = ''),main = paste('Probability of at ',ifelse(res$increase, 'least', 'most'),' ', n.needed.for.greater.p0,' patients with ',outcomeLabel,sep = ''),lwd = 3,cex.lab = 1.5,ylim = c(0,1),cex.main = 2,cex.lab = 1.5,cex.axis = 1.5,cex = 1.5)
    text(aa1,tmp5[,2] + .1,round(tmp5[,2],2),col = 2,cex = 2)
    par(mfrow = c(1,1))
    
  })
  
  output$plotSensitivityInd <- renderUI({
    n.list <- get.result.futility()$n.list
    plot_output_list <- lapply(1:(length(n.list) - 1), function(j) {
      plotname <- paste("SenIndplot", j, sep = "")
      plotOutput(plotname)
      #plotOutput(plotname, height = 280, width = 250)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  observe({
    res <-  get.result.futility()
    tmp5 <- res$sensitivity.data
    tmp5.ind <- res$sensitivity.data.ind
    n.list <- res$n.list
    pTarget <- res$pTarget
    theta <- res$theta
    cutoffPredictive <- res$cutoffPredictive
    outcomeLabel <- res$outcomeLabel
    armLabel <- res$armLabel
    betaA <- res$betaA
    betaB <- res$betaB
    n.predictive.cutoff <- res$n.predictive.cutoff
    n.needed.for.greater.p0 <- res$n.needed.for.greater.p0
    
    dim.n <- dim(tmp5.ind)[2]
    
    #If just one interim analysis this plot is redundant, so excluding
    if (dim.n > 3) {
      for (j in 1:(dim.n - 1))
      {
        local({
          my_j <- j
          plotname <- paste("SenIndplot", my_j, sep = "")
          output[[plotname]] <- renderPlot({
            aa1 <- barplot(tmp5.ind[,my_j],xlab = paste('Rate of ',outcomeLabel,sep = ''),ylab = paste('Probability of early stopping the trial at interim',j),main = paste('Interim Analysis: ',my_j,'\n Probability of early stopping the trial',sep = ''),lwd = 3,cex.lab = 1.5,ylim = c(0,1),cex.main = 2,cex.lab = 1.5,cex.axis = 1.5,cex = 1.5)
            text(aa1,tmp5.ind[,my_j],round(tmp5.ind[,my_j],2),col = 2,cex = 2)
          })
        })
      }
    }
  })
  
  
  
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste('bayesian_futility_by_predictive_probability', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
    }, 
    content = function(file) {
      out <- render('bayesian_futility_by_predictive_probability.Rmd', switch(
        input$format,
        PDF = pdf_document(), HTML = html_document(), Word = word_document()
      ))
      file.rename(out, file)
    }
  )
  
}

shinyApp(ui = ui,server = server)