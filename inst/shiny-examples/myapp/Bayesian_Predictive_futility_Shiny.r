library(shiny)
library(knitr)
library(rmarkdown)
library(extraDistr)




ui <- fluidPage(
  headerPanel('Interim Analysis for Futility Based on Bayesian Predictive Probability'),
  tabsetPanel(
    tabPanel("Interim Analysis for Futility Based on Bayesian Predictive Probability  (P1>P0)",
             sidebarLayout(
               sidebarPanel(
                 textInput('nlistA', 'Enter sample size for each interim analysis (comma delimited; e.g., 25,25)', "25,25"),
                 numericInput("pTargetA",label = "Probability of  unfavorable (p0)",min = 0, max = 1, value = 0.3),
                 numericInput("pH1A",label = "Probability of  favorable (p1: H1 with p1>p0)",min = 0, max = 1, value = 0.5),
                 numericInput("thetaA",label = "Cutoff for greater p0 (prob(x>p0|data)>theta",min = 0, max = 1,value = 0.95),
                 numericInput("cutoffPredictiveA",label = "Cutoff for predictive probability",min = 0, max = 1,value = 0.05),
                 numericInput("betaAA",label = "beta A parameter",min = 0, max = 1,value = 1),
                 numericInput("betaBA",label = "beta B parameter",min = 0, max = 1,value = 1),
                 numericInput("simNA",label = "number of simulation",value = 10000),
                 textInput('outcomeLabelA', label = 'Name of outcome', value = "response"),
                 textInput('armLabelA', label = 'Name of the Arm', value = "treatment"),
                 actionButton("SubmitA","Calculate"),
                 radioButtons('formatA', 'Document format', c('Word','PDF', 'HTML' ),
                              inline = TRUE),
                 downloadButton('downloadReportIncrease')
               ),
               
               mainPanel(
                 h2("Interim Analysis for Futility"),
                 textOutput("futilityA"),
                 tags$hr(),
                 
                 h4("Table 0: Stopping Boundary for Futility"),
                 tableOutput("stoppingBoundaryA"),
                 tags$hr(),
                 
                 
                 h4("Table 1: Bayesian Predictive Probability"),
                 uiOutput("predictiveDataA"),
                 tags$hr(),
                 
                 h4("Table 2: Sensitivity Analysis"),
                 tableOutput("sensitivityDataA"),
                 tags$hr(),
                 
                 #h4("Figure 11: Bayesian Predictive Probability"),
                 #plotOutput("plotPredictiveB"),
                 
                 # This is the dynamic UI for the plots
                 h4("Figure 1: Bayesian Predictive Probability"),
                 uiOutput("plotPredictiveA"),
                 
                 h4("Figure 2: Sensitivity"),
                 plotOutput("plotSensitivityOverallA"),
                 
                 h4("Figure 3: Probability of Stopping Early by Each Interim Analysis"),
                 uiOutput("plotSensitivityIndA"),
                 
                 tags$hr()
               ))),
    #---this is for (P1<p0)-------------------    
    tabPanel("Interim Analysis for Futility Based on Bayesian Predictive Probability  (P1<P0)",
             sidebarLayout(
               sidebarPanel(
                 textInput('nlistB', 'Enter sample size for each interim analysis (comma delimited)', "26,26,23"),
                 numericInput("pTargetB",label = "Probability of  unfavorable (p0)",min = 0, max = 1, value = 0.27),
                 numericInput("pH1B",label = "Probability of  favorable (p1: H1 with p1<p0)",min = 0, max = 1, value = 0.17),
                 numericInput("thetaB",label = "Cutoff for greater p0 (prob(x<p0|data)>theta",min = 0, max = 1,value = 0.9),
                 numericInput("cutoffPredictiveB",label = "Cutoff for predictive probability",min = 0, max = 1,value = 0.05),
                 numericInput("betaAB",label = "beta A parameter",min = 0, max = 1,value = 1),
                 numericInput("betaBB",label = "beta B parameter",min = 0, max = 1,value = 1),
                 numericInput("simNB",label = "number of simulation",value = 10000),
                 textInput('outcomeLabelB', label = 'Name of outcome', value = "event"),
                 textInput('armLabelB', label = 'Name of the Arm', value = "treatment"),
                 actionButton("SubmitB","Calculate"),
                 radioButtons('formatB', 'Document format', c('Word','PDF', 'HTML' ),
                              inline = TRUE),
                 downloadButton('downloadReportReduction')
               ),

               mainPanel(
                 h2("Interim Analysis for Futility"),
                 textOutput("futilityB"),
                 tags$hr(),

                 h4("Table 1: Bayesian Predictive Probability"),
                 uiOutput("predictiveDataB"),
                 tags$hr(),

                 h4("Table 2: Sensitivity Analysis"),
                 tableOutput("sensitivityDataB"),
                 tags$hr(),

                 #h4("Figure 11: Bayesian Predictive Probability"),
                 #plotOutput("plotPredictiveB"),

                 # This is the dynamic UI for the plots
                 h4("Figure 1: Bayesian Predictive Probability"),
                 uiOutput("plotPredictiveB"),

                 h4("Figure 2: Sensitivity"),
                 plotOutput("plotSensitivityOverallB"),

                 h4("Figure 3: Probability of Stopping Early by Each Interim Analysis"),
                 uiOutput("plotSensitivityIndB"),
                 tags$hr()
               )))
    
  )
)


server <- function(input,output){
  
  
  #---Futility: p1>p0--------
  
  get.result.futilityA <- eventReactive(input$SubmitA, {
    n.list <- as.numeric(unlist(strsplit(input$nlistA,",")))
    pTarget <- input$pTargetA
    pH1 <- input$pH1A
    theta <- input$thetaA
    cutoffPredictive <- input$cutoffPredictiveA
    outcomeLabel <- input$outcomeLabelA
    armLabel <- input$armLabelA
    betaA <- input$betaAA
    betaB <- input$betaBA
    simN <- input$simNA
    if (pH1 > pTarget)
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
    
    
    a1 <- paste('UPUPUPUPFutility evaluation will be implementded in ',length(n.predictive.cutoff),ifelse(length(n.predictive.cutoff) == 1,' interim analysis ',' interim analyses '),' for the first ',paste(cumsum(n.list[-length(n.list)]),collapse = ', '),' patients. ',' A Bayesian approach for futility analysis is used to calculate posterior probability and predictive probability for the rate of ',outcomeLabel,' with a non-informative beta prior, beta(',betaA,',',betaB,'). We consider a ',round(pTarget*100),'% rate or higher of ',outcomeLabel, ' as ineffective for the treatment. Thus, we expect the treatment ',armLabel,' is not worse than the historical control if the posterior probability of the rate (',outcomeLabel,') greater than ',round(pTarget*100),'% is higher than ',theta,' (i.e., prob(rate of ',outcomeLabel,'>',round(pTarget*100),'% |data)>',theta,') ). With a total ',sum(n.list),' patients in treatment ',armLabel,', the number of patients with ',outcomeLabel,' needs to be ',n.needed.for.greater.p0,' or more in order to meet the criteria. Therefore, we use the number of ',n.needed.for.greater.p0,' patients to guide the predictive probability. Specifically, for the first ',n.list[1],' patients in the 1st interim analysis, there are ',n.list[1] + 1,' ways for number of patients with ',outcomeLabel,' from 0, 1, to, ',n.list[1],'. In each case, given the number of patients with ',outcomeLabel,', $s$, in the first ',n.list[1],' patients, we calculate predictive probability of $',n.needed.for.greater.p0,'-s$ or more patients with ',outcomeLabel,' in the future remaining ',sum(n.list) - n.list[1],sep = '')
    nn1 <- n.list[1]
    nn2 <- sum(n.list) - n.list[1]
    n.predictive.cutoff.text <- n.predictive.cutoff
    n.predictive.cutoff.text[n.predictive.cutoff.text == (-99)] <- 'Pass'
    
    index.example <- n.predictive.cutoff == (-99)
    nn1.example <- sum((n.list[-length(n.list)])[index.example]) + (n.list[-length(n.list)])[!index.example][1]
    nn2.example <- sum(n.list) - nn1.example
    n.predictive.cutoff.example <- n.predictive.cutoff[!index.example][1]
    
    a1 <- paste(a1,' patients, i.e., $\\sum_{i=',n.needed.for.greater.p0,'-s}^{',nn2,'} {',nn2,' \\choose i } \\frac{beta(',betaA,'+s+i,',betaB,'+(',nn1,'-s)+(',nn2,'-i))}{beta(',betaA,'+s,',betaB,'+(',nn1,'-s))}$.
              Calculation of predictive probability is based on beta binominal distribution for the number of patients with ',outcomeLabel,' in the future remaining ',nn2,' patients given a beta distribution for the  rate of ',outcomeLabel,', $beta(',betaA,'+s,',betaB,'+',nn1,'-s)$. 
              For example, if there are ',n.predictive.cutoff.example[1],' patients with ',outcomeLabel,' in the first ',nn1.example,' patients, the predictive probability of ',n.needed.for.greater.p0 - n.predictive.cutoff.example[1],' or more patients with ',outcomeLabel,' in the future remaining ',nn2.example,' patients would be 
              $\\sum_{i=',n.needed.for.greater.p0 - n.predictive.cutoff.example[1],'}^{',nn2.example,'}{',nn2.example,' \\choose i } \\frac{beta(',betaA,'+',n.predictive.cutoff.example[1],'+i,',betaB,'+(',nn1,'-',n.predictive.cutoff.example[1],')+(',nn2.example,'-i))}{beta(',betaA,'+',n.predictive.cutoff.example[1],',',betaB,'+(',nn1.example,'-',n.predictive.cutoff.example[1],'))} =',round(predictive.prob.data[[1]][predictive.prob.data[[1]][,3] == n.predictive.cutoff.example[1],5],3),'$.
              The predictive probability is also calculated for each of the remaining interim analyses to evaluate the chance of ',n.needed.for.greater.p0,'-s or more patients with ',outcomeLabel,' in the future remaining patients given s patients with ',outcomeLabel,' in the current stage of interim analysis. ',
              ' Figure 1 and Table 1 lists predictive probability for all scenarios of number of patients with ',outcomeLabel,' in each interim analysis and the associated largest number of patients with ',outcomeLabel,' needed in the future remaining patients to have at least a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,'. ', 
              ' We consider that a ',round(cutoffPredictive*100),'% cutoff of the predictive probability will give little chance to have ',n.needed.for.greater.p0,' patients or more with ',outcomeLabel,' at the end of study. Thus with this cutoff, the stopping rule (Table 0) will be: the trial will be stopped if there are ',ifelse(length(n.list) == 2,paste(n.predictive.cutoff,' or less patients with ',outcomeLabel,' in the 1st interim analysis. ',sep = ''),paste(paste(paste(n.predictive.cutoff.text[-length(n.predictive.cutoff)],collapse = ', '),n.predictive.cutoff.text[length(n.predictive.cutoff)],sep = ', and '),' or less patients with ',outcomeLabel,' in the 1st to ',stage.index,' interim analysis, respectively. ',sep = '')),
              ' Sensitivity analysis (Figure 2 and Table 2) for this stopping rule shows that if the true rate of ',outcomeLabel,' is ',round(pTarget*100),'%, the chance to reach a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,2],2),' (Type I error). On the other hand, the probability to stop the trial early is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,1],2)*100 ,'%. When the true rate of ',outcomeLabel,' is ',as.numeric(pH1)*100,'%, the chance to reach at least ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[paste(pH1),2],2),'. The corresponding probability to stop the treatment early is ',round(sensitivity.data[paste(pH1),1],2)*100,'%.',sep = '')
    
    
    a2 <- paste('Application of Bayesian predictive probability for interim analysis in single-arm early phase II trial. Chen et al; submitted.')
    
    tmp98 <- list(a1 = a1,
                a2 = a2,
                n.needed.for.greater.p0 = n.needed.for.greater.p0,
                n.predictive.cutoff = n.predictive.cutoff,
                all.data = all.data,
                predictive.prob.data = predictive.prob.data,
                sensitivity.data = round(sensitivity.data,4),
                sensitivity.data.ind = round(sensitivity.data.ind,4),
                n.list = n.list,
                pTarget = pTarget,
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
  
  
  output$futilityA <- renderText({
    get.result.futilityA()$a1
  })
  
  output$stoppingBoundaryA <- renderTable({
    aa0 <- get.result.futilityA()
    aa <- aa0$n.list
    cutoff.tmp <- c(aa0$n.predictive.cutoff,aa0$n.needed.for.greater.p0 - 1)
    cutoff.tmp[cutoff.tmp == (-99)] <- 'PASS'
    aa <- rbind(c(1:(length(aa) - 1),'Final'),cumsum(aa),aa,cutoff.tmp)
    dimnames(aa)[[1]] <- c('Stage of interim analysis','Sample size up to the current stage','Sample size at each stage','Stopping boundary')
    aa
  },rownames = T,colnames = F, digits = 0)
  
  
  output$predictiveDataA <- renderUI({
    n.list <- get.result.futilityA()$n.list
    table_output_list <- lapply(1:(length(n.list) - 1), function(i) {
      tablename <- paste("PredictiveTable", i, sep = "")
      tableOutput(tablename)
    })
    do.call(tagList, table_output_list)
  })
  
  observe({
    aa0 <- get.result.futilityA()
    aa.list <- aa0$predictive.prob.data
    
    for (j in 1:length(aa.list))
    {
      local({
        my_j <- j
        tablename <- paste("PredictiveTable", my_j, sep = "")
        stage.index <- ifelse(my_j == 1,' 1st ',ifelse(my_j == 2,' 2nd ',ifelse(my_j == 3,' 3rd ',paste(' ',my_j,'th ',sep = ''))))
        aa <- aa.list[[my_j]][,3:5,drop = F]
        outcomeLabel <- aa0$outcomeLabel
        dimnames(aa)[[2]] <- c(paste('number of patients with ', outcomeLabel,' in the ',stage.index,' interim analysis',sep = ''),paste('minimum number of patients with ',outcomeLabel,' needed in the furture remaining patients',sep = ''),'predictive probability')
        output[[tablename]] <- renderTable({
          aa
        })
      })
    }
  })
  
  
  output$sensitivityDataA <- renderTable({
    aa0 <- get.result.futilityA()
    aa <- aa0$sensitivity.data
    aa <- cbind(as.numeric(dimnames(aa)[[1]]),aa)
    dimnames(aa)[[2]] <- c('true rate',' overall probability of early stopping the trial',paste('probability to have at least ',aa0$n.needed.for.greater.p0,' patients with ',aa0$outcomeLabel,sep = ''))
    aa
    
  })
  
  output$plotPredictiveA <- renderPlot({
    res <-  get.result.futilityA()
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
      mtext(paste('minimum number of patients with ',outcomeLabel,' needed in the furture remaining ', n2,' patients',sep = ''),side = 3,line = 1.8)
      title(paste(stage.index,' Interim Analysis of for Futility \n\n',sep = ''))
      arrows(kk2.len, 1, kk2.len, y1 = 0.1,col = 3,lwd = 3)
      text(1:dim(tmp3)[1],tmp3[,5],round(tmp3[,5],2),col = 2,cex = 2)
    }
    
  })
  
  
  # Insert the right number of plot output objects into the web page
  output$plotPredictiveA <- renderUI({
    n.list <- get.result.futilityA()$n.list
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
    res <-  get.result.futilityA()
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
          mtext(paste('minimum number of patients with ',outcomeLabel,' needed in the furture remaining ', n2,' patients',sep = ''),side = 3,line = 1.8)
          title(paste(stage.index,' Interim Analysis of for Futility \n\n',sep = ''))
          arrows(kk2.len, 1, kk2.len, y1 = 0.1,col = 3,lwd = 3)
          text(1:dim(tmp3)[1],tmp3[,5],round(tmp3[,5],2),col = 2,cex = 2)
        })
      })
    }
  })  
  
  output$plotSensitivityOverallA <- renderPlot({
    res <-  get.result.futilityA()
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
    
    aa1 <- barplot(tmp5[,2],xlab = paste('Rate of ',outcomeLabel,sep = ''),ylab = paste('Probability of at least ', n.needed.for.greater.p0,' patients with ', outcomeLabel,sep = ''),main = paste('Probability of at least ', n.needed.for.greater.p0,' patients with ',outcomeLabel,sep = ''),lwd = 3,cex.lab = 1.5,ylim = c(0,1),cex.main = 2,cex.lab = 1.5,cex.axis = 1.5,cex = 1.5)
    text(aa1,tmp5[,2] + .1,round(tmp5[,2],2),col = 2,cex = 2)
    par(mfrow = c(1,1))
    
  })
  
  output$plotSensitivityIndA <- renderUI({
    n.list <- get.result.futilityA()$n.list
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
    res <-  get.result.futilityA()
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
  })
  
  
  
  output$downloadReportIncrease <- downloadHandler(
    filename = function() {
      paste('bayesian_futility_by_predictive_probability', sep = '.', switch(
        input$formatA, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
    }, 
    content = function(file) {
      out <- render('bayesian_futility_by_predictive_probability.Rmd', switch(
        input$formatA,
        PDF = pdf_document(), HTML = html_document(), Word = word_document()
      ))
      file.rename(out, file)
    }
  )
  
  
  
  
  #---Futility: p1<p0--------
  
  get.result.futilityB <- eventReactive(input$SubmitB, {
    n.list <- as.numeric(unlist(strsplit(input$nlistB,",")))
    pTarget <- input$pTargetB
    pH1 <- input$pH1B
    theta <- input$thetaB
    cutoffPredictive <- input$cutoffPredictiveB
    outcomeLabel <- input$outcomeLabelB
    armLabel <- input$armLabelB
    betaA <- input$betaAB
    betaB <- input$betaBB
    simN <- input$simNB
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
    
    #    a0 <- paste('i.e., $\\sum_{i=',nn2,'-s}^{',nn2,'} {',nn2,' \\choose i \\frac{beta(',betaA,'+s + i,',betaB,'+(',nn1,'-s)+(',nn2,'-i))}{beta(',betaA,'+s,',betaB,'+(',nn1,'-s))}} = 2%$.',sep = '') 
    
    
    a1 <- paste('DOWNDOWNDOWNFutility evaluation will be implementded in ',length(n.predictive.cutoff),ifelse(length(n.predictive.cutoff) == 1,' interim analysis ',' interim analyses '),' for the first ',paste(cumsum(n.list[-length(n.list)]),collapse = ', '),' patients. ',' A Bayesian approach for futility analysis is used to calculate posterior probability and predictive probability for the rate of ',outcomeLabel,' with a non-informative beta prior, beta(',betaA,',',betaB,'). We consider a ',round(pTarget*100),'% rate or higher of ',outcomeLabel, ' as ineffective for the treatment. Thus, we expect the treatment ',armLabel,' is not worse than the historical control if the posterior probability of the rate (',outcomeLabel,') less than ',round(pTarget*100),'% is higher than ',theta,' (i.e., prob(rate of ',outcomeLabel,'<',round(pTarget*100),'% |data)>',theta,') ). With a total ',sum(n.list),' patients in treatment ',armLabel,', the number of patients with ',outcomeLabel,' needs to be ',n.needed.for.greater.p0,' or less in order to meet the criteria. Therefore, we use the number of ',n.needed.for.greater.p0,' patients to guide the predictive probability. Specifically, for the first ',n.list[1],' patients in the 1st interim analysis, there are ',n.list[1] + 1,' ways for number of patients with ',outcomeLabel,' from 0, 1, to, ',n.list[1],'. In each case, given the number of patients with ',outcomeLabel,', $s$, in the first ',n.list[1],' patients, we calculate predictive probability of $',n.needed.for.greater.p0,'-s$ or less patients with ',outcomeLabel,' in the future remaining ',sum(n.list) - n.list[1],sep = '')
    nn1 <- n.list[1]
    nn2 <- sum(n.list) - n.list[1]
    
    
    a1 <- paste(a1,' patients, i.e., $\\sum_{i=',n.needed.for.greater.p0,'-s}^{',0,'} {',nn2,' \\choose i } \\frac{beta(',betaA,'+s+i,',betaB,'+(',nn1,'-s)+(',nn2,'-i))}{beta(',betaA,'+s,',betaB,'+(',nn1,'-s))}$. 
              Calculation of predictive probability is based on beta binominal distribution for the number of patients with ',outcomeLabel,' in the future remaining ',nn2,' patients given a beta distribution for the  rate of ',outcomeLabel,', $beta(',betaA,'+s,',betaB,'+',nn1,'-s)$. 
              For example, if there are ',n.predictive.cutoff[1],' patients with ',outcomeLabel,' in the first ',nn1,' patients, the predictive probability of ',n.needed.for.greater.p0 - n.predictive.cutoff[1],' or less patients with ',outcomeLabel,' in the future remaining ',nn2,' patients would be 
              $\\sum_{i=',n.needed.for.greater.p0 - n.predictive.cutoff[1],'}^{',0,'} {',nn2,' \\choose i } \\frac{beta(',betaA,'+',n.predictive.cutoff[1],'+i,',betaB,'+(',nn1,'-',n.predictive.cutoff[1],')+(',nn2,'-i))}{beta(',betaA,'+',n.predictive.cutoff[1],',',betaB,'+(',nn1,'-',n.predictive.cutoff[1],'))} =',round(predictive.prob.data[[1]][predictive.prob.data[[1]][,3] == n.predictive.cutoff[1],5],3),'$.
              The predictive probability is also calculated for the 2nd interim analysis to evaluate the chance of ',n.needed.for.greater.p0,'-s or less patients with event in the future remaining patients given s patients with event in the current stage of interim analysis. ',
              ' Figure 1 and Table 1 lists predictive probability for all scenarios of number of patients with event in each interim analysis and the associated largest number of patients with event needed in the future remaining patients to have at most a total of ',n.needed.for.greater.p0,' patients with event. ', 
              ' We consider that a ',round(cutoffPredictive*100),'% cutoff of the predictive probability will give little chance to have ',n.needed.for.greater.p0,' patients or less with event at the end of study. Thus with this cutoff, the stopping rule will be: the trial will be stopped if there are ',ifelse(length(n.list) == 2,paste(n.predictive.cutoff,' or more patients with event in the 1st interim analysis. ',sep = ''),paste(paste(n.predictive.cutoff,collapse = ' and '),' or more patients with event in the 1st to ',stage.index,' interim analysis, respectively. ',sep = '')),
              ' Sensitivity analysis (Figure 2 and Table 2) for this stopping rule shows that if the true rate of ',outcomeLabel,' is ',round(pTarget*100),'%, the chance to reach a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,2],2),' (Type I error). On the other hand, the probability to stop the trial early is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,1],2)*100 ,'%. When the true rate of ',outcomeLabel,' is ',as.numeric(pH1)*100,'%, the chance to reach at most ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[paste(pH1),2],2),'. The corresponding probability to stop the treatment early is ',round(sensitivity.data[paste(pH1),1],2)*100,'%.',sep = '')
    
    
    a2 <- paste('Saville BR, Connor JT, Ayers GD, et al: The utility of Bayesian predictive probabilities for interim monitoring of clinical trials. Clin Trials 11:485-493, 2014')
    
    tmp98 <- list(a1 = a1,
                a2 = a2,
                n.needed.for.greater.p0 = n.needed.for.greater.p0,
                n.predictive.cutoff = n.predictive.cutoff,
                all.data = all.data,
                predictive.prob.data = predictive.prob.data,
                sensitivity.data = round(sensitivity.data,4),
                sensitivity.data.ind = round(sensitivity.data.ind,4),
                n.list = n.list,
                pTarget = pTarget,
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
  
  
  output$futilityB <- renderText({
    get.result.futilityB()$a1
  })
  
  
  output$predictiveDataB <- renderUI({
    n.list <- get.result.futilityB()$n.list
    table_output_list <- lapply(1:(length(n.list) - 1), function(i) {
      tablename <- paste("PredictiveTable", i, sep = "")
      tableOutput(tablename)
    })
    do.call(tagList, table_output_list)
  })
  
  observe({
    aa0 <- get.result.futilityB()
    aa.list <- aa0$predictive.prob.data
    
    for (j in 1:length(aa.list))
    {
      local({
        my_j <- j
        tablename <- paste("PredictiveTable", my_j, sep = "")
        stage.index <- ifelse(my_j == 1,' 1st ',ifelse(my_j == 2,' 2nd ',ifelse(my_j == 3,' 3rd ',paste(' ',my_j,'th ',sep = ''))))
        aa <- aa.list[[my_j]][,3:5]
        outcomeLabel <- aa0$outcomeLabel
        dimnames(aa)[[2]] <- c(paste('number of patients with ', outcomeLabel,' in the ',stage.index,' interim analysis',sep = ''),paste('largest number of patients with ',outcomeLabel,' needed in the future remaining patients',sep = ''),'predictive probability')
        output[[tablename]] <- renderTable({
          aa
        })
      })
    }
  })
  
  
  output$sensitivityDataB <- renderTable({
    aa0 <- get.result.futilityB()
    aa <- aa0$sensitivity.data
    aa <- cbind(as.numeric(dimnames(aa)[[1]]),aa)
    dimnames(aa)[[2]] <- c('true rate',' overall probability of early stopping the trial',paste('probability to have at most ',aa0$n.needed.for.greater.p0,' patients with ',aa0$outcomeLabel,sep = ''))
    aa
    
  })
  
  output$plotPredictiveB <- renderPlot({
    res <-  get.result.futilityB()
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
    if (1 > 2)
    {
      tmp3 <- res$predictive.prob.data
      kk2 <- tmp3[tmp3[,5] < cutoffPredictive,1]
      kk2.len <- length(kk2)
      
      plot(tmp3[,5],type = 'h',axes = F,xlab = paste('number of patients with ', outcomeLabel,' in the 1st ',nn1,' patients',sep = ''),ylab = 'predictive probability',lwd = 3,cex.lab = 1.5)
      box();axis(2,cex.axis = 1.2)
      axis(2,cutoffPredictive,cex.axis = 1.2)
      axis(1,1:dim(tmp3)[1],tmp3[,3],cex.axis = 1.2)
      axis(3,1:dim(tmp3)[1],tmp3[,4],cex.axis = 1)
      abline(h = cutoffPredictive,col = 2,lty = 2)
      mtext(paste('largest number of patients with ',outcomeLabel,' needed in the furture remaining ', nn2,' patients',sep = ''),side = 3,line = 1.8)
      title(paste('Interim Analysis of Futility for Treatment ',armLabel,'\n\n',sep = ''))
      arrows(kk2.len, 1, kk2.len, y1 = 0.1,col = 3,lwd = 3)
      text(1:dim(tmp3)[1],tmp3[,5],round(tmp3[,5],2),col = 2,cex = 2)
    }
    
    par(mfrow = c(2,ceiling(length(tmp3.list)/2)))
    for (j in 1:length(tmp3.list)) 
    {
      tmp3 <- tmp3.list[[j]]
      n1 <- sum(n.list[1:j])
      n2 <- sum(n.list[(j + 1):length(n.list)])
      kk2.len <- kk2.len.list[[j]]
      n.predictive.cutoff <- n.predictive.cutoff.list[j]
      stage.index <- ifelse(j == 1,' 1st ',ifelse(j == 2,' 2nd ',ifelse(j == 3,' 3rd ',paste(' ',j,'th ',sep = ''))))
      plot(tmp3[,5],type = 'h',axes = F,xlab = paste('number of patients with ', outcomeLabel,' in the',stage.index,' interim analysis with ',n1,' patients',sep = ''),ylab = 'predictive probability',lwd = 3,cex.lab = 1.5)
      box();axis(2,cex.axis = 1.2)
      axis(2,predictive.cutoff,cex.axis = 1.2)
      axis(1,1:dim(tmp3)[1],tmp3[,3],cex.axis = 1.2)
      axis(3,1:dim(tmp3)[1],tmp3[,4],cex.axis = 1)
      abline(h = predictive.cutoff,col = 2,lty = 2)
      mtext(paste('largest number of patients with ',outcomeLabel,' needed in the furture remaining ', n2,' patients',sep = ''),side = 3,line = 1.8)
      title(paste(stage.index,' Interim Analysis of for Futility \n\n',sep = ''))
      arrows(kk2.len, 1, kk2.len, y1 = 0.1,col = 3,lwd = 3)
      text(1:dim(tmp3)[1],tmp3[,5],round(tmp3[,5],2),col = 2,cex = 2)
    }
    
  })
  
  
  
  # Insert the right number of plot output objects into the web page
  output$plotPredictiveB <- renderUI({
    n.list <- get.result.futilityB()$n.list
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
    res <-  get.result.futilityB()
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
          plot(tmp3[,5],type = 'h',axes = F,xlab = paste('number of patients with ', outcomeLabel,' in the',stage.index,' interim analysis with ',n1,' patients',sep = ''),ylab = 'predictive probability',lwd = 3,cex.lab = 1.5)
          box();axis(2,cex.axis = 1.2)
          axis(2,predictive.cutoff,cex.axis = 1.2)
          axis(1,1:dim(tmp3)[1],tmp3[,3],cex.axis = 1.2)
          axis(3,1:dim(tmp3)[1],tmp3[,4],cex.axis = 1)
          abline(h = predictive.cutoff,col = 2,lty = 2)
          mtext(paste('largest number of patients with ',outcomeLabel,' needed in the furture remaining ', n2,' patients',sep = ''),side = 3,line = 1.8)
          title(paste(stage.index,' Interim Analysis of for Futility \n\n',sep = ''))
          arrows(kk2.len, 1, kk2.len, y1 = 0.1,col = 3,lwd = 3)
          text(1:dim(tmp3)[1],tmp3[,5],round(tmp3[,5],2),col = 2,cex = 2)
        })
      })
    }
  })  
  
  output$plotSensitivityOverallB <- renderPlot({
    res <- get.result.futilityB()
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
    
    aa1 <- barplot(tmp5[,2],xlab = paste('Rate of ',outcomeLabel,sep = ''),ylab = paste('Probability of at most ', n.needed.for.greater.p0,' patients with ', outcomeLabel,sep = ''),main = paste('Probability of at most ', n.needed.for.greater.p0,' patients with ',outcomeLabel,sep = ''),lwd = 3,cex.lab = 1.5,ylim = c(0,1),cex.main = 2,cex.lab = 1.5,cex.axis = 1.5,cex = 1.5)
    text(aa1,tmp5[,2] + .1,round(tmp5[,2],2),col = 2,cex = 2)
    par(mfrow = c(1,1))
    
  })
  
  output$plotSensitivityIndB <- renderUI({
    n.list <- get.result.futilityB()$n.list
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
    res <-  get.result.futilityB()
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
  })
  
  
  
  output$downloadReportReduction <- downloadHandler(
    filename = function() {
      paste('bayesian_futility_by_predictive_probability', sep = '.', switch(
        input$formatA, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
    }, 
    
    content = function(file) {
      out <- render('bayesian_futility_by_predictive_probability.Rmd', switch(
        input$formatA,
        PDF = pdf_document(), HTML = html_document(), Word = word_document()
      ))
      file.rename(out, file)
    }
  )
  
  
}

shinyApp(ui = ui,server = server)