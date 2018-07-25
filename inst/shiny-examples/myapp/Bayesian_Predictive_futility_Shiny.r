library(shiny)
library(knitr)
library(rmarkdown)
library(extraDistr)



ui <- fluidPage(
  headerPanel('Interim Analysis for Futility Based on Bayesian Predictive Probability'),

  tabsetPanel(

    tabPanel("Interim Analysis for Futility Based on Bayesian Predictive Probability",
             sidebarLayout(
               sidebarPanel(
                 numericInput("num1", label ="Sample size for 1st interim analysis", value = 30),
                 numericInput("num2", label ="Sample size for the 2nd interim analysis", value = 38),
                 numericInput("pTarget",label="Probability of least favorable (p0)",min = 0, max = 1, value = 0.45),
                 numericInput("pH1",label="Probability of most favorable (p1: H1)",min = 0, max = 1, value = 0.65),
                 numericInput("theta",label="Cutoff for greater p0 (prob(x>p0|data)>theta",min = 0, max = 1,value = 0.95),
                 numericInput("cutoffPredictive",label="Cutoff for predictive probability",min = 0, max = 1,value = 0.05),
                 numericInput("betaA",label="beta A parameter",min = 0, max = 1,value = 1),
                 numericInput("betaB",label="beta B parameter",min = 0, max = 1,value = 1),
                 numericInput("simN",label="number of simulation",value = 10000),
                 textInput('outcomeLabel', label='Name of outcome', value = "response"),
                 textInput('armLabel', label='Name of the Arm', value = "treatment"),
                 actionButton("Submit1","Calculate"),
                 radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                              inline = TRUE),
                 downloadButton('downloadReport')
               ),

               mainPanel(
                 h2("Interim Analysis for Futility"),
                 textOutput("futility"),
                 tags$hr(),

                 h4("Table 1: Bayesian Predictive Probability"),
                 tableOutput("predictiveData"),
                 tags$hr(),

                 h4("Table 2: Sensitivity Analysis"),
                 tableOutput("sensitivityData"),
                 tags$hr(),

                 h4("Figure 1: Bayesian Predictive Probability"),
                 plotOutput("plotPredictive"),

                 h4("Figure 2: Sensitivity"),
                 plotOutput("plotSensitivity"),

                 h4("Reference"),
                 verbatimTextOutput("reference"),
                 tags$hr()

               )))
  )
)


server <- function(input,output){

  #---Futility--------

  get.result.futility <- eventReactive(input$Submit1, {
    nn1<-input$num1
    nn2<-input$num2
    pTarget<-input$pTarget
    pH1<-input$pH1
    theta<-input$theta
    cutoffPredictive<-input$cutoffPredictive
    outcomeLabel<-input$outcomeLabel
    armLabel<-input$armLabel
    betaA<-input$betaA
    betaB<-input$betaB
    simN<-input$simN
    print(outcomeLabel)
    tmp99<-Bayesian.predictive.futility.fun(n1=nn1,n2=nn2,
                                                p.target=pTarget,
                                                p.h1=pH1,
                                                cutoff.n.for.greater.p0=theta,
                                                predictive.cutoff=cutoffPredictive,
                                                beta.a=betaA,
                                                beta.b=betaB,
                                                sim.n=simN,
                                                outcome.tmp=outcomeLabel,
                                                arm.name=armLabel)
    n.needed.for.greater.p0<-tmp99$n.needed.for.greater.p0
    n.predictive.cutoff<-tmp99$n.predictive.cutoff
    all.data<-tmp99$all.data
    predictive.prob.data<-tmp99$predictive.prob.data
    sensitivity.data<-tmp99$sensitivity.data
    n.additional.needed<-n.needed.for.greater.p0-n.predictive.cutoff

    a0<-paste('i.e., $\\sum_{i=',nn2,'-s}^{',nn2,'} {',nn2,' \\choose i \\frac{beta(',betaA,'+s+i,',betaB,'+(',nn1,'-s)+(',nn2,'-i))}{beta(',betaA,'+s,',betaB,'+(',nn1,'-s))}} = 2%$.',sep='')


    a1<-paste('One interim analysis for futility will be done for the arm ',armLabel,' after the first ',nn1, ' patients of the arm ',armLabel,'.  A Bayesian approach for futility analysis is used to calculate posterior probability and predictive probability for the rate of ',outcomeLabel,' with a non-informative beta prior, beta(',betaA,',',betaB,'). We consider a ',round(pTarget*100),'% rate or less of ',outcomeLabel, ' as ineffective. Thus, we expect arm ',armLabel,' is not worse than the historical control if the posterior probability of the rate (',outcomeLabel,') greater than ',round(pTarget*100),'% is higher than ',theta,' (i.e., prob(rate of ',outcomeLabel,'>',round(pTarget*100),'% |data)>',theta,') ). With a total ',nn1+nn2,' patients in arm ',armLabel,', it will need at least ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' to meet the criteria. Therefore, we use the number of ',n.needed.for.greater.p0,' patients to guide the predictive probability. Specifically, for the first ',nn1,' patients in the interim analysis, there are ',nn1+1,' ways for number of patients with ',outcomeLabel,' from 0, 1, to, ',nn1,'. In each case, given the number of patients with ',outcomeLabel,', $s$, in the first ',nn1,' patients, we calculate predictive probability of $',n.needed.for.greater.p0,'-s$ or more patients with ',outcomeLabel,' in the future remaining ',nn2,sep='')
a1<-paste(a1,' patients, i.e., $\\sum_{i=',n.needed.for.greater.p0,'-s}^{',nn2,'} {',nn2,' \\choose i } \\frac{beta(',betaA,'+s+i,',betaB,'+(',nn1,'-s)+(',nn2,'-i))}{beta(',betaA,'+s,',betaB,'+(',nn1,'-s))}$.
Calculation of predictive probability is based on beta binominal distribution for the number of patients with ',outcomeLabel,' in the future remaining ',nn2,' patients given a beta distribution for the  rate of ',outcomeLabel,', $beta(',betaA,'+s,',betaB,'+',nn1,'-s)$.
For example, if there are ',n.predictive.cutoff,' patients with ',outcomeLabel,' in the first ',nn1,' patients, the predictive probability of ',n.needed.for.greater.p0-n.predictive.cutoff,' or more patients with ',outcomeLabel,' in the future remaining ',nn2,' patients would be
$\\sum_{i=',n.needed.for.greater.p0-n.predictive.cutoff,'}^{',nn2,'} {',nn2,' \\choose i } \\frac{beta(',betaA,'+',n.predictive.cutoff,'+i,',betaB,'+(',nn1,'-',n.predictive.cutoff,')+(',nn2,'-i))}{beta(',betaA,'+',n.predictive.cutoff,',',betaB,'+(',nn1,'-',n.predictive.cutoff,'))} =',round(predictive.prob.data[predictive.prob.data[,3]==n.predictive.cutoff,5],3),'$.
Figure 1 lists predictive probability for all scenarios of number of patients with ',outcomeLabel,' in the first ',nn1,' patients (interim analysis) and number of patients with ',outcomeLabel,'  needed in the future remaining ',nn2,' patients to have at least a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,'. With a ',round(cutoffPredictive*100),'% cutoff for the predictive probability, the stopping rule is if there are ',         n.predictive.cutoff,' or less patients with ',outcomeLabel,' for the first ',nn1,' patients in the interim analysis (i.e., little chance to reach a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at the end of study), we consider arm ',armLabel,' is ineffective and will be stopped.   Sensitivity analysis (Figure 2) for this stopping rule shows that if the true rate of ',outcomeLabel,' is ',round(pTarget*100),'%, the chance to reach a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]])==pTarget,2],2),'. The probability to stop the arm early is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]])==pTarget,1],2)*100 ,'%. When the true rate of ',outcomeLabel,' is ',as.numeric(dimnames(sensitivity.data)[[1]][1])*100,'%, the chance to reach a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[1,2],2),'. The corresponding probability to stop the arm early increases to ',round(sensitivity.data[1,1],2)*100,'%.',sep='')

    a2<-paste('Chen, et al: Application of Bayesian predictive probability for interim analysis in single-arm early phase II trial. submitted')


    tmp98<-list(a1=a1,
                a2=a2,
                n.needed.for.greater.p0=n.needed.for.greater.p0,
                n.predictive.cutoff=n.predictive.cutoff,
                all.data=round(all.data,4),
                predictive.prob.data=round(predictive.prob.data,4),
                sensitivity.data=round(sensitivity.data,4),
                nn1=nn1,
                nn2=nn2,
                pTarget=pTarget,
                theta=theta,
                cutoffPredictive=cutoffPredictive,
                outcomeLabel=outcomeLabel,
                armLabel=armLabel,
                betaA=betaA,
                betaB=betaB
                )
    print(outcomeLabel)
    return(tmp98)
  })


  output$futility <- renderText({
    get.result.futility()$a1
  })


  output$predictiveData <- renderTable({
    aa0<-get.result.futility()
    aa<-aa0$predictive.prob.data[,3:5]
    outcomeLabel<-aa0$outcomeLabel
    dimnames(aa)[[2]]<-c(paste('number of patients with ', outcomeLabel,' in the 1st interim analysis',sep=''),paste('number of patients with ',outcomeLabel,' needed in the furture remaining patients',sep=''),'predictive probability')
    aa
  })

  output$sensitivityData <- renderTable({
     #get.result.futility()$sensitivity.data
    aa0<-get.result.futility()
    aa<-aa0$sensitivity.data
    aa<-cbind(as.numeric(dimnames(aa)[[1]]),aa)
    dimnames(aa)[[2]]<-c('true rate','probability of early stopping the arm',paste('probability to have at least ',aa0$n.needed.for.greater.p0,' patients with ',aa0$outcomeLabel,sep=''))
    aa

  })

  output$plotPredictive<-renderPlot({
    res <-  get.result.futility()
    nn1<-res$nn1
    nn2<-res$nn2
    pTarget<-res$pTarget
    theta<-res$theta
    cutoffPredictive<-res$cutoffPredictive
    outcomeLabel<-res$outcomeLabel
    armLabel<-res$armLabel
    betaA<-res$betaA
    betaB<-res$betaB

    n.predictive.cutoff<-res$n.predictive.cutoff

    tmp3<-res$predictive.prob.data
    kk2<-tmp3[tmp3[,5]<cutoffPredictive,1]
    kk2.len<-length(kk2)

      plot(tmp3[,5],type='h',axes=F,xlab=paste('number of patients with ', outcomeLabel,' in the 1st ',nn1,' patients',sep=''),ylab='predictive probability',lwd=3,cex.lab=1.5)
      box();axis(2,cex.axis=1.2)
      axis(2,cutoffPredictive,cex.axis=1.2)
      axis(1,1:dim(tmp3)[1],tmp3[,3],cex.axis=1.2)
      axis(3,1:dim(tmp3)[1],tmp3[,4],cex.axis=1)
      abline(h=cutoffPredictive,col=2,lty=2)
      mtext(paste('number of patients with ',outcomeLabel,' needed in the furture remaining ', nn2,' patients',sep=''),side=3,line=1.8)
      title(paste('Interim Analysis of Futility for Arm ',armLabel,'\n\n',sep=''))
      arrows(kk2.len, 1, kk2.len, y1 = 0.1,col=3,lwd=3)
      text(1:dim(tmp3)[1],tmp3[,5],round(tmp3[,5],2),col=2,cex=2)

  })

  output$plotSensitivity<-renderPlot({
    res <-  get.result.futility()
    tmp5<-res$sensitivity.data
    nn1<-res$nn1
    nn2<-res$nn2
    pTarget<-res$pTarget
    theta<-res$theta
    cutoffPredictive<-res$cutoffPredictive
    outcomeLabel<-res$outcomeLabel
    armLabel<-res$armLabel
    betaA<-res$betaA
    betaB<-res$betaB
    n.predictive.cutoff<-res$n.predictive.cutoff
    n.needed.for.greater.p0<-res$n.needed.for.greater.p0

    par(mfrow=c(2,1))
    aa1<-barplot(tmp5[,1],xlab=paste('Rate of ',outcomeLabel,sep=''),ylab='Probability of early stopping the arm',main='Sensitivity Analysis \n Probability of early stopping the arm',lwd=3,cex.lab=1.5,ylim=c(0,1),cex.main=2,cex.lab=1.5,cex.axis=1.5,cex=1.5)
    text(aa1,tmp5[,1],round(tmp5[,1],2),col=2,cex=2)

    aa1<-barplot(tmp5[,2],xlab=paste('Rate of ',outcomeLabel,sep=''),ylab=paste('Probability of at least ', n.needed.for.greater.p0,' patients with ', outcomeLabel,sep=''),main=paste('Probability of at least ', n.needed.for.greater.p0,' patients with ',outcomeLabel,sep=''),lwd=3,cex.lab=1.5,ylim=c(0,1),cex.main=2,cex.lab=1.5,cex.axis=1.5,cex=1.5)
    text(aa1,tmp5[,2]+.1,round(tmp5[,2],2),col=2,cex=2)
    par(mfrow=c(1,1))

  })

  output$reference <- renderText({
    get.result.futility()$a2
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

shinyApp(ui=ui,server=server)
