#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



# Define UI for application that draws a histogram
ui <- fluidPage(
  headerPanel('Interim Analysis for Futility Based on Bayesian Predictive Probability'),
  sidebarLayout(
    sidebarPanel(
      textInput('project_title', 'Project Title', "My Awesome Project"),
      textInput('author_input', 'Author(s)', "Dung-Tsa Chen"),
      textInput('nlist', 'Enter sample size for each interim analysis (comma delimited; e.g., 25,25 for 1st and 2nd stage)', "25,25"),
      numericInput("pTarget", label = "Probability of  unfavorable under null hypothesis (p0) (e.g., p0=30% response rate at H0)", min = 0, max = 1, value = 0.3, step = 0.05),
      numericInput("pH1", label = "Probability of  favorable under alternative hypothesis (p1) (e.g., p1=50% response rate at H1)", min = 0, max = 1, value = 0.5, step = 0.05),
      numericInput("theta", label = "Threshold of posterior probability to define treatment efficacy (0-100%; suggested value: 80-99%; a higher threshold requires a larger number of responders to claim efficacy and leads to a lower power)", min = 0, max = 1, value = 0.95, step = 0.05),
      numericInput("cutoffPredictive", label = "Cutoff of predictive probability to define stopping rule (0-100%; suggested value: 0.01-0.3; a higher threshold requires a larger number of responders to advance to the next stage and leads to terminate the trial early)", min = 0, max = 1, value = 0.05, step = 0.05),
      numericInput("betaA", label = "beta A parameter (degree of response in treatment)", min = 0, max = 1, value = 1, step = 0.05),
      numericInput("betaB", label = "beta B parameter (degree of non-response in treatment)", min = 0, max = 1, value = 1, step = 0.05),
      radioButtons('analysis_type', 'Analysis Type', c('Analytical', 'Simulation'),
                   inline = TRUE),
      conditionalPanel(
        condition = "input.analysis_type == 'Simulation'",
        numericInput("simN", label = "number of simulations", min = 1, value = 10000, step = 1000),
        numericInput("seed", label = "seed", min = 1, max = 2147483647, value = 42534253, step = 1)
      ),
      textInput('outcomeLabel', label = 'Name of outcome', value = "response"),
      textInput('armLabel', label = 'Name of the Arm', value = "treatment"),
      actionButton("CalculateButton", "Calculate"),
      radioButtons('format', 'Document format', c('Word', 'PDF'),
                   inline = TRUE),
      downloadButton('downloadReport')
    ),
    mainPanel(
              h2("Summary of Interim Analysis for Futility"),
              span(textOutput("summary"), style="color:red"),
#              textOutput("summary"),
              tags$hr(),

              h3("Details"),
              textOutput("futility"),
              tags$hr(),

              h4("Table 1: Stopping Boundary for Futility"),
              tableOutput("stoppingBoundary"),
              tags$hr(),

              h4("Table 2: Bayesian Predictive Probability for Stopping Rule"),
              uiOutput("predictiveData"),
              tags$hr(),

              h4("Table 3: Performance (Probability of Early Termination (PET), Type I error, and Power)"),
              tableOutput("sensitivityData"),
              tags$hr(),

              h4("Table 4: Sensitivity Analysis: Predictive Probability"),
              tableOutput("SensitivityAnalysisPredictive"),
              tags$hr(),

              h4("Table 5: Sensitivity Analysis: Posterior Probability"),
              tableOutput("SensitivityAnalysisPosterior"),
              tags$hr(),

              h4("Table 6: Sensitivity Analysis: Sample Size"),
              tableOutput("SensitivityAnalysisSampleSize"),
              tags$hr(),

              h4("Table 7: Sensitivity Analysis: Beta Prior Distribution"),
              tableOutput("SensitivityAnalysisPrior"),
              tags$hr(),


              # This is the dynamic UI for the plots
              h4("Figure 1: Bayesian Predictive Probability for Stopping Rule"),
              uiOutput("plotPredictive"),

              h4("Figure 2: Performance (Probability of Early Termination (PET), Type I error, and Power)"),
              plotOutput("plotSensitivityOverall"),

              h4("Figure 3: Probability of Stopping Early by Each Interim Analysis"),
              uiOutput("plotSensitivityInd"),

              h4("Figure 4: Sensitivity Analysis: Predictive Probability"),
              plotOutput("plotSensitivityAnalysisPredictive"),

              h4("Figure 5: Sensitivity Analysis: Posterior Probability"),
              plotOutput("plotSensitivityAnalysisPosterior"),

              h4("Figure 6: Sensitivity Analysis: Sample Size"),
              plotOutput("plotSensitivityAnalysisSampleSize"),

              h4("Figure 7: Sensitivity Analysis: Beta Prior Distribution"),
              plotOutput("plotSensitivityAnalysisPrior"),

              tags$hr()
    ))
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  #---Futility--------
  get.result.futility <- eventReactive(input$CalculateButton, {
    # logging::logdebug("testing")

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
    seed <- input$seed
    analysis_type <- input$analysis_type
    simN <- input$simN
    if (input$project_title == '') project_title = ' ' else project_title = input$project_title
    if (input$author_input == '') author_input = ' ' else author_input = input$author_input

    tmp99 <- evaluateInterim(
      ns = n.list,
      p.target = pTarget,
      p.h1 = pH1,
      cutoff.n.for.greater.p0 = theta,
      predictive.cutoff = cutoffPredictive,
      beta.a = betaA,
      beta.b = betaB,
      analysis.type = analysis_type,
      sim.n = simN,
      seed = seed,
      outcome.tmp = outcomeLabel,
      arm.name = armLabel
    )

    if(length(tmp99)==1)
    {
      a1_summary<-tmp99
      tmp98 <- list(a1_summary = a1_summary)
    }
    else
    {
    #print(tmp99$n.predictive.cutoff)
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

    #print('#--sensitivity analysis----')
  #--sensitivity analysis----
    p0<-pTarget
    p1<-pH1
    cutoff.n.for.greater.p0.tmp<-theta
    predictive.cutoff.tmp<-cutoffPredictive
    ns.list<-n.list
    beta.prior.list<-c(betaA,betaB)

    sensitivity.analysis.ans<-list()

    #---sample size----
    #print('Sample size')
    parameter.value<-c(-5:5)
    var.name.tmp<-paste('n',1:length(ns.list),sep='')
    tmp0<-numeric()

    tmp0<-numeric()
    for(i in 1:length(parameter.value))
    {
      a1=evaluateInterim(ns=ns.list+parameter.value[i],p.target = p0,p.h1=p1,predictive.cutoff=predictive.cutoff.tmp,beta.a =beta.prior.list[1],beta.b = beta.prior.list[2],cutoff.n.for.greater.p0=cutoff.n.for.greater.p0.tmp)

      if(length(a1)!=1)
      {
      tmp1<-c(ns.list+parameter.value[i],a1$sensitivity.data[paste(p0),'prob.stop.overall'],
              a1$sensitivity.data[paste(p0),'prob.above.threshold.of.reject.Ho'],
              a1$sensitivity.data[paste(p1),'prob.above.threshold.of.reject.Ho'])
      names(tmp1)<-c(var.name.tmp,'PET','typeI','power')
      tmp0<-rbind(tmp0,tmp1)
      }
    }
    rownames(tmp0)<-NULL
    sensitivity.analysis.ans$sample.size<-list(ans=tmp0,var.name.tmp=var.name.tmp)

    #---predictive probability----
    #print('predictive probability')
    parameter.value<-unique(sort(c(0.01,0.05,0.1,0.15,0.2,0.25,0.3,predictive.cutoff.tmp)))
    var.name.tmp<-'Cutoff of predictive probability'
    tmp0<-numeric()
    for(i in 1:length(parameter.value))
    {
      a1=evaluateInterim(
        p.target = p0,p.h1=p1,
        predictive.cutoff=parameter.value[i],
        cutoff.n.for.greater.p0=cutoff.n.for.greater.p0.tmp,
        ns=ns.list,
        beta.a =beta.prior.list[1],
        beta.b = beta.prior.list[2]
      )

      if(length(a1)!=1)
      {
      tmp1<-c(parameter.value[i],a1$sensitivity.data[paste(p0),'prob.stop.overall'],
              a1$sensitivity.data[paste(p0),'prob.above.threshold.of.reject.Ho'],
              a1$sensitivity.data[paste(p1),'prob.above.threshold.of.reject.Ho'])
      names(tmp1)<-c(var.name.tmp,'PET','typeI','power')
      tmp0<-rbind(tmp0,tmp1)
      }
    }
    rownames(tmp0)<-NULL
    sensitivity.analysis.ans$predictive.prob<-list(ans=tmp0,var.name.tmp=var.name.tmp)

    #---threshold of posterior probability----
    #cutoff.n.for.greater.p0.tmp<-.95
    #print('posterior')
    parameter.value<-unique(sort(c((80:99)/100,cutoff.n.for.greater.p0.tmp)))
    var.name.tmp<-'threshold of posterior probability'
    tmp0<-numeric()
    for(i in 1:length(parameter.value))
    {
      a1=evaluateInterim(
        p.target = p0,p.h1=p1,
        predictive.cutoff=predictive.cutoff.tmp,
        cutoff.n.for.greater.p0=parameter.value[i],
        ns=ns.list,
        beta.a =beta.prior.list[1],
        beta.b = beta.prior.list[2]
      )

      if(length(a1)!=1)
      {
      tmp1<-c(parameter.value[i],a1$sensitivity.data[paste(p0),'prob.stop.overall'],
              a1$sensitivity.data[paste(p0),'prob.above.threshold.of.reject.Ho'],
              a1$sensitivity.data[paste(p1),'prob.above.threshold.of.reject.Ho'])
      names(tmp1)<-c(var.name.tmp,'PET','typeI','power')
      tmp0<-rbind(tmp0,tmp1)
      }
    }
    rownames(tmp0)<-NULL
    sensitivity.analysis.ans$posterior.prob<-list(ans=tmp0,var.name.tmp=var.name.tmp)

    #---prior distribution----
    #print('prior')
    #beta.a.tmp<-1
    #beta.b.tmp<-1
    estBetaParams <- function(mu, var) {
      # mu: 0-1
      #var: 0-0.5^2
      alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
      beta <- alpha * (1 / mu - 1)
      return(params = list(alpha = alpha, beta = beta))
    }

    sd1<-c(.05,.1,.2,.3)
    parameter.value<-rbind(c(0,0),c(0,1),c(1,0),c(1,1),t(apply(as.matrix(sd1^2),1,function(x) unlist(estBetaParams(mu=p0,var=x)))),t(apply(as.matrix(sd1^2),1,function(x) unlist(estBetaParams(mu=p1,var=x)))))
    pos.value.status<-apply(parameter.value,1,function(x) all(x>=0))
    parameter.value<-parameter.value[pos.value.status,]

    prior.redundant.status<-apply(parameter.value,1,function(x) all(x==beta.prior.list))
    if(!any(prior.redundant.status))
      parameter.value<-rbind(parameter.value,beta.prior.list)
#print('print(prior.redundant.status)')
#print(parameter.value)
    var.name.tmp<-c('beta.a','beta.b')
    tmp0<-ii<-numeric()
    for(i in 1:dim(parameter.value)[1])
    {
      a1=evaluateInterim(
        p.target = p0,p.h1=p1,
        predictive.cutoff=predictive.cutoff.tmp,
        cutoff.n.for.greater.p0=cutoff.n.for.greater.p0.tmp,
        ns=ns.list,
        beta.a =parameter.value[i,1],
        beta.b = parameter.value[i,2]
      )

      if(length(a1)!=1)
      {
        ii<-c(ii,i)
        #print(ii)
        #print(tmp0)
      tmp1<-c(parameter.value[i,1:2],a1$sensitivity.data[paste(p0),'prob.stop.overall'],
              a1$sensitivity.data[paste(p0),'prob.above.threshold.of.reject.Ho'],
              a1$sensitivity.data[paste(p1),'prob.above.threshold.of.reject.Ho'])
      names(tmp1)<-c(var.name.tmp,'PET','typeI','power')
      tmp0<-rbind(tmp0,tmp1)
      }
    }

    #print(tmp0)
    #print('--B--')
    if(any(prior.redundant.status)) dimnames(tmp0)[[1]]<-c(apply(parameter.value[1:4,],1,paste,collapse='/'),paste(p0,' (SD=',sd1,')',sep=''),paste(p1,' (SD=',sd1,')',sep=''))[ii] else
      dimnames(tmp0)[[1]]<-c(apply(parameter.value[1:4,],1,paste,collapse='/'),paste(p0,' (SD=',sd1,')',sep=''),paste(p1,' (SD=',sd1,')',sep=''),paste(beta.prior.list,collapse='/'))[ii]
#print(tmp0)
    #rownames(tmp0)<-NULL
    sensitivity.analysis.ans$prior<-list(ans=tmp0,var.name.tmp=var.name.tmp)

    sensitivity.analysis.summary<-sapply(sensitivity.analysis.ans,function(y)
      apply(y$ans,2,function(x) paste(range(round(x,2)),collapse='-')),simplify=F)


    a1 <- paste0('A Bayesian approach for futility analysis is used to calculate posterior probability and predictive probability for the rate of ',outcomeLabel,' with a non-informative beta prior, beta(',betaA,',',betaB,'), using ', ifelse(analysis_type == 'Analytical', 'the analytical form', paste0('a simulation, with ', simN,' replications (Seed:',seed,')')),'. We consider a ',round(pTarget*100),'% rate or ',ifelse(increase, 'lower', 'higher'),' of ',outcomeLabel, ' as ineffective for the treatment. Thus, we expect the ',armLabel,' arm is promising if the posterior probability of the rate (',outcomeLabel,') ',ifelse(increase, 'greater', 'less'),' than ',round(pTarget*100),'% is higher than ',theta,' (i.e., prob(rate of ',outcomeLabel,ifelse(increase, '>', '<'),round(pTarget*100),'% |data)>',theta,') ). \n \n With a total ',sum(n.list),' patients in ',armLabel,' arm, the number of patients with ',outcomeLabel,' needs to be ',n.needed.for.greater.p0,' or ',ifelse(increase, 'more', 'less'),' in order to meet the criteria. Therefore, we use the number of ',n.needed.for.greater.p0,' patients to guide the predictive probability. Specifically, given the number of patients with ',outcomeLabel,', $s$, in the first ',n.list[1],' patients, we calculate predictive probability of $',n.needed.for.greater.p0,'-s$ or ',ifelse(increase, 'more', 'less'),' patients with ',outcomeLabel,' in the future remaining ',sum(n.list) - n.list[1])

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
                 ' \n \n We consider that a ',round(cutoffPredictive*100),'% cutoff of the predictive probability will give unlikely chance to have ',n.needed.for.greater.p0,' patients or ',ifelse(increase, 'more', 'less'),' with ',outcomeLabel,' at the end of study. Thus with this cutoff, the stopping rule (Table 1) will be: the trial will be stopped if there are ',ifelse(length(n.list) == 2,paste(n.predictive.cutoff,' or ',ifelse(increase, 'less', 'more'),' patients with ',outcomeLabel,' in the 1st interim analysis. ',sep = ''),paste(paste(paste(n.predictive.cutoff.text[-length(n.predictive.cutoff)],collapse = ', '),n.predictive.cutoff.text[length(n.predictive.cutoff)],sep = ', and '),' or ',ifelse(increase, 'less', 'more'),' patients with ',outcomeLabel,' in the 1st to ',stage.index,' interim analysis, respectively. ',sep = '')),
                 ' Performance of this stopping rule (Figure 2 and Table 3) shows that if the true rate of ',outcomeLabel,' is ',round(pTarget*100),'%, the chance to reach at ',ifelse(increase, 'least', 'most'),' a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,2] * 100),'% (Type I error), however the probability of early termination (PET) is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,1],2)*100 ,'%.  When the true rate of ',outcomeLabel,' is ',as.numeric(pH1)*100,'%, then the chance to reach at ',ifelse(increase, 'least', 'most'),' ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[paste(pH1),2] * 100),'% (power), and the corresponding probability to stop the treatment early is ',round(sensitivity.data[paste(pH1),1],2)*100,'%.',ifelse(ncol(sensitivity.data.ind) >= 3,'  Figure 3 shows the probability of stopping at each interim analysis.',''))


    #print(sensitivity.analysis.summary)
a1<-paste0(a1,' \n \n Sensitivity analysis (Table 4-7 and Figure 4-7) evaluates four parameters for their impact on performance (PET, type I error, and power): cutoff for the predictive probability, threshold for posterior probability of response rate, sample size, and beta prior distribtuion of the response rate.
Evaluation is conducted for each parameter when the values of other parameters are fixed.
When the cutoff of the predictive probability for the stopping rule is ',paste(range(sensitivity.analysis.ans$predictive.prob$ans[,1]),collapse='-'),', the range is ', sensitivity.analysis.summary$predictive.prob['PET'],' for PET, ',sensitivity.analysis.summary$predictive.prob['typeI'],' for type I error, and ',sensitivity.analysis.summary$predictive.prob['power'],' for power (Table 4 and Figure 4). ',sep='')

a1<-paste0(a1,' When the threshold for posterior probability is ',paste(range(sensitivity.analysis.ans$posterior.prob$ans[,1]),collapse='-'),', the range is ', sensitivity.analysis.summary$posterior.prob['PET'],' for PET, ',sensitivity.analysis.summary$posterior.prob['typeI'],' for type I error, and ',sensitivity.analysis.summary$posterior.prob['power'],' for power (Table 5 and Figure 5). ',sep='')

n.diff<-range(sensitivity.analysis.ans$sample.size$ans[,1]-ns.list[1])

a1<-paste0(a1,' When the sample size of each stage is in the magnitude from ', ifelse(prod(n.diff)<0,paste0('decrease by ',n.diff[1],' to increase by ',n.diff[2],sep=''),paste('increase by ',paste(n.diff,collapse='-'),sep='')),', the range is ', sensitivity.analysis.summary$sample.size['PET'],' for PET, ',sensitivity.analysis.summary$sample.size['typeI'],' for type I error, and ',sensitivity.analysis.summary$sample.size['power'],' for power (Table 6 and Figure 6). ',sep='')

#a1<-paste0(a1,' When the sample size of each stage is in the magnitude from decrease of 5 to increase of 5, the range is ', sensitivity.analysis.summary$sample.size['PET'],' for PET, ',sensitivity.analysis.summary$sample.size['typeI'],' for type I error, and ',sensitivity.analysis.summary$sample.size['power'],' for power (Table 6 and Figure 6). ',sep='')

a1<-paste0(a1,' When the beta prior varies from non-informative prior to  the one with a response rate at the null or alternative hypothesis and  a series of standard deviation (SD), the range is ', sensitivity.analysis.summary$prior['PET'],' for PET, ',sensitivity.analysis.summary$prior['typeI'],' for type I error, and ',sensitivity.analysis.summary$prior['power'],' for power (Table 7 and Figure 7). ',sep='')

    a1_summary <- paste0('Futility evaluation is implemented in ',ifelse(length(n.predictive.cutoff) == 1, paste0(' one interim analysis with ',n.list[1],' patients',sep=''),paste0(length(n.predictive.cutoff),paste0(' interim analyses with ',paste((n.list[-length(n.list)]),collapse = ', '),' patients in stage ',paste(paste(1:length(n.predictive.cutoff),collapse=', '),', respectively',sep=''),sep=''),sep='')),', and ', n.list[length(n.list)],' patients in the last stage, for a total of ',sum(n.list),' patients. With an unfavorable rate set at ',round(pTarget*100),'% (null hypothesis) and posterior probability of ', theta,' as the threshold, a total of at least ', n.needed.for.greater.p0, ' of the ', sum(n.list),' patients must have ', outcomeLabel, ' to be able to claim treatment efficacy. Given a ',round(cutoffPredictive*100),'% cutoff of the predictive probability (i.e. chance to stop the trial in the ', ifelse(length(n.predictive.cutoff) == 1,' interim analysis',' interim analyses'), ') the stopping rule (Table 1) will be: the trial will be stopped if there are ',ifelse(length(n.list) == 2,paste(n.predictive.cutoff,' or ',ifelse(increase, 'less', 'more'),' patients with ',outcomeLabel,' in the 1st interim analysis. ',sep = ''),paste(paste(paste(n.predictive.cutoff.text[-length(n.predictive.cutoff)],collapse = ', '),n.predictive.cutoff.text[length(n.predictive.cutoff)],sep = ', and '),' or ',ifelse(increase, 'less', 'more'),' patients with ',outcomeLabel,' in the 1st to ',stage.index,' interim analysis, respectively. ',sep = '')), ' \n \n Performance of the design (Figure 2 and Table 3) shows that if the true rate of ',outcomeLabel,' is ',round(pTarget*100),'%, the chance to reach at ',ifelse(increase, 'least', 'most'),' a total of ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,2] * 100),'% (Type I error), however the probability to stop the trial early is ',round(sensitivity.data[as.numeric(dimnames(sensitivity.data)[[1]]) == pTarget,1],2)*100 ,'%.  If true rate of ',outcomeLabel,' is indeed ',as.numeric(pH1)*100,'%, then the chance to reach at ',ifelse(increase, 'least', 'most'),' ',n.needed.for.greater.p0,' patients with ',outcomeLabel,' at end of the study is ',round(sensitivity.data[paste(pH1),2] * 100),'% (power), and the corresponding probability to stop the treatment early is ',round(sensitivity.data[paste(pH1),1],2)*100,'%.')

      a2 <- c('Chen et al, Application of Bayesian predictive probability for interim analysis in single-arm early phase II trial. Submitted',
  'Lee JJ and Liu DD. A predictive probability design for phase II cancer clinical trials. Clinical trials. 2008; 5: 93-106.')

    tmp98 <- list(project_title = project_title,
                  author_input = author_input,
                  a1_summary = a1_summary,
                  a1 = a1,
                  a2 = a2,
                  n.needed.for.greater.p0 = n.needed.for.greater.p0,
                  n.predictive.cutoff = n.predictive.cutoff,
                  all.data = all.data,
                  predictive.prob.data = predictive.prob.data,
                  sensitivity.data = round(sensitivity.data,4),
                  sensitivity.data.ind = round(sensitivity.data.ind,4),
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
                  kk2.len.list = kk2.len.list,sensitivity.analysis.ans=sensitivity.analysis.ans
    )
    }

    #print('End function')
    return(tmp98)
  })

  output$summary <- renderText({
    get.result.futility()$a1_summary
  })

  output$futility <- renderText({
    aa0<-get.result.futility()
    if(length(aa0)==1) '' else
      aa0$a1
  })

  output$stoppingBoundary <- renderTable({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    aa <- aa0$n.list
    cutoff.tmp <- c(aa0$n.predictive.cutoff,aa0$n.needed.for.greater.p0 - 1)
    cutoff.tmp[cutoff.tmp == (-99)] <- 'PASS'
    aa <- rbind(c(1:(length(aa) - 1),'Final'),cumsum(aa),aa,cutoff.tmp)
    dimnames(aa)[[1]] <- c('Stage of interim analysis','Sample size up to the current stage','Sample size at each stage','Stopping boundary')
    aa
    }
  },rownames = T,colnames = F, digits = 0)


  output$predictiveData <- renderUI({
    aa0<-get.result.futility()
    if(length(aa0)==1) '' else
    {
    n.list <- get.result.futility()$n.list
    table_output_list <- lapply(1:(length(n.list) - 1), function(i) {
      tablename <- paste("PredictiveTable", i, sep = "")
      tableOutput(tablename)
    })
    do.call(tagList, table_output_list)
    }
  })

  observe({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    aa.list <- aa0$predictive.prob.data

    for (j in 1:length(aa.list))
    {
      local({
        my_j <- j
        tablename <- paste("PredictiveTable", my_j, sep = "")
        stage.index <- ifelse(my_j == 1,' 1st ',ifelse(my_j == 2,' 2nd ',ifelse(my_j == 3,' 3rd ',paste(' ',my_j,'th ',sep = ''))))
        aa <- aa.list[[my_j]][,3:5,drop = (input$pH1 < input$pTarget)]
        aa[,dim(aa)[2]] <- round(aa[,dim(aa)[2]],3)
        outcomeLabel <- aa0$outcomeLabel
        dimnames(aa)[[2]] <- c(paste('number of patients with ', outcomeLabel,' in the ',stage.index,' interim analysis',sep = ''),paste(ifelse(aa0$increase, 'minimum', 'maximum'),' number of patients with ',outcomeLabel,' ',ifelse(aa0$increase, 'needed', 'allowed'),' in the future remaining patients',sep = ''),'predictive probability')
        output[[tablename]] <- renderTable({
          aa
        })
      })
    }
    }
  })


  output$sensitivityData <- renderTable({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    aa <- aa0$sensitivity.data
    aa <- cbind(as.numeric(dimnames(aa)[[1]]),aa)
    dimnames(aa)[[2]] <- c('true rate',' overall probability of early stopping the trial',paste('probability to have at ',ifelse(aa0$increase, 'least', 'most'),' ',aa0$n.needed.for.greater.p0,' patients with ',aa0$outcomeLabel,sep = ''))
    aa
    }

  })

  output$plotPredictive <- renderPlot({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
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
      text(1:dim(tmp3)[1],pmin(tmp3[,5] + .1, .9),round(tmp3[,5],2),col = 2,cex = 1.75)
    }
    }

  })


  # Insert the right number of plot output objects into the web page
  output$plotPredictive <- renderUI({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    n.list <- get.result.futility()$n.list
    plot_output_list <- lapply(1:(length(n.list) - 1), function(j) {
      plotname <- paste("plot", j, sep = "")
      plotOutput(plotname)
      #plotOutput(plotname, height = 280, width = 250)
    })

    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
    }
  })

  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  observe({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
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
          text(1:dim(tmp3)[1],pmin(tmp3[,5] + .1, .8),round(tmp3[,5],2),col = 2,cex = 2)
        })
      })
    }
    }
  })

  output$plotSensitivityOverall <- renderPlot({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
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
    aa1 <- barplot(tmp5[,1],xlab = paste('Rate of ',outcomeLabel,sep = ''),ylab = 'Probability of early stopping the trial overall',main = 'Overall probability of early stopping the trial (PET)',lwd = 3,cex.lab = 1.5,ylim = c(0,1),cex.main = 2,cex.lab = 1.5,cex.axis = 1.5,cex = 1.5)
    text(aa1,pmin(tmp5[,1], .8),round(tmp5[,1],2),col = 2,cex = 2)

    aa1 <- barplot(tmp5[,2],xlab = paste('Rate of ',outcomeLabel,sep = ''),ylab = paste('Probability of at ',ifelse(res$increase, 'least', 'most'),' ', n.needed.for.greater.p0,' patients with ', outcomeLabel,sep = ''),main = paste('Probability of at ',ifelse(res$increase, 'least', 'most'),' ', n.needed.for.greater.p0,' patients with ',outcomeLabel,sep = ''),lwd = 3,cex.lab = 1.5,ylim = c(0,1),cex.main = 2,cex.lab = 1.5,cex.axis = 1.5,cex = 1.5)
    text(aa1,pmin(tmp5[,2] + .1, .8),round(tmp5[,2],2),col = 2,cex = 2)
    par(mfrow = c(1,1))
    }

  })

  output$plotSensitivityInd <- renderUI({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    n.list <- get.result.futility()$n.list
    plot_output_list <- lapply(1:(length(n.list) - 1), function(j) {
      plotname <- paste("SenIndplot", j, sep = "")
      plotOutput(plotname)
      #plotOutput(plotname, height = 280, width = 250)
    })

    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
    }
  })

  observe({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
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
    if (dim.n > 2) {
      for (j in 1:(dim.n - 1))
      {
        local({
          my_j <- j
          plotname <- paste("SenIndplot", my_j, sep = "")
          output[[plotname]] <- renderPlot({
            aa1 <- barplot(tmp5.ind[,my_j],xlab = paste('Rate of ',outcomeLabel,sep = ''),ylab = paste('Probability of early stopping the trial at interim',j),main = paste('Interim Analysis: ',my_j,'\n Probability of early stopping the trial',sep = ''),lwd = 3,cex.lab = 1.5,ylim = c(0,1),cex.main = 2,cex.lab = 1.5,cex.axis = 1.5,cex = 1.5)
            text(aa1,pmin(tmp5.ind[,my_j] + .1,.8),round(tmp5.ind[,my_j],2),col = 2,cex = 2)
          })
        })
      }
    }
    }
  })

#--sensitivity analysis-----
  #sensitivity.analysis.ans$predictive.prob
  #sensitivity.analysis.ans$posterior.prob
  #sensitivity.analysis.ans$prior
  #sensitivity.analysis.ans$sample.size
  #tableOutput("SensitivityAnalysisPredictive")
  #tableOutput("SensitivityAnalysisPosterior")
  #tableOutput("SensitivityAnalysisSampleSize")
  #tableOutput("SensitivityAnalysisPrior")


  output$SensitivityAnalysisPredictive <- renderTable({
    #print('SensitivityAnalysisPredictive')
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    aa <- aa0$sensitivity.analysis.ans$predictive.prob
    aa$ans
    }
  })

  output$SensitivityAnalysisPosterior <- renderTable({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    aa <- aa0$sensitivity.analysis.ans$posterior.prob
    aa$ans
    }
  })

  output$SensitivityAnalysisSampleSize <- renderTable({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    aa <- aa0$sensitivity.analysis.ans$sample.size
    aa$ans
    }
  })

  output$SensitivityAnalysisPrior <- renderTable({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    aa <- aa0$sensitivity.analysis.ans$prior
    ans<-data.frame(prior=dimnames(aa$ans)[[1]],aa$ans)
    ans
    }
  })

  output$plotSensitivityAnalysisPredictive <- renderPlot({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    res <-  get.result.futility()
    aa0<-res$sensitivity.analysis.ans$predictive.prob
    tmp0<-aa0$ans
    var.name.tmp<-aa0$var.name.tmp
    #print('plotSensitivityAnalysisPredictive')

    par(mfrow=c(3,1))
    name2<-c('PET','typeI','power')
    name21<-c('PET','type I error','power')
    for(i in 1:length(name2))
    {
      plot(1:dim(tmp0)[1],tmp0[,name2[i]],xlab=var.name.tmp,ylab=name21[i],main=name21[i],axes=F)
      box()
      axis(2)
      axis(1,1:dim(tmp0)[1],tmp0[,var.name.tmp])
      index1<-tmp0[,1]==res$cutoffPredictive
      points((1:dim(tmp0)[1])[index1],tmp0[index1,name2[i]],col=2,cex=4,pch=2)

    }
    }
  })

  output$plotSensitivityAnalysisPosterior <- renderPlot({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    res <-  get.result.futility()
    aa0<-res$sensitivity.analysis.ans$posterior.prob
    tmp0<-aa0$ans
    var.name.tmp<-aa0$var.name.tmp

    par(mfrow=c(3,1))
    name2<-c('PET','typeI','power')
    name21<-c('PET','type I error','power')
    for(i in 1:length(name2))
    {
      plot(1:dim(tmp0)[1],tmp0[,name2[i]],xlab=var.name.tmp,ylab=name21[i],main=name21[i],axes=F)
      box()
      axis(2)
      axis(1,1:dim(tmp0)[1],tmp0[,var.name.tmp])
      index1<-tmp0[,1]==res$theta
      points((1:dim(tmp0)[1])[index1],tmp0[index1,name2[i]],col=2,cex=4,pch=2)

    }
    }
  })

  output$plotSensitivityAnalysisSampleSize <- renderPlot({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    res <-  get.result.futility()
    aa0<-res$sensitivity.analysis.ans$sample.size
    tmp0<-aa0$ans
    var.name.tmp<-aa0$var.name.tmp
    par(mfrow=c(3,1))
    name2<-c('PET','typeI','power')
    name21<-c('PET','type I error','power')
    for(i in 1:length(name2))
    {
      plot(1:dim(tmp0)[1],tmp0[,name2[i]],xlab='sample size',ylab=name21[i],main=name21[i],axes=F)
      box()
      axis(2)
      nn1<-apply(tmp0[,var.name.tmp],1,paste,collapse='/')
      axis(1,1:dim(tmp0)[1],nn1)
      index1<-nn1==paste(res$n.list,collapse ='/')
      points((1:dim(tmp0)[1])[index1],tmp0[index1,name2[i]],col=2,cex=4,pch=2)
    }
    }
    })

  output$plotSensitivityAnalysisPrior <- renderPlot({
    aa0 <- get.result.futility()
    if(length(aa0)==1) '' else
    {
    res <-  get.result.futility()
    aa0<-res$sensitivity.analysis.ans$prior
    tmp0<-aa0$ans
    var.name.tmp<-aa0$var.name.tmp
    beta.prior.list<-c(res$betaA,res$betaB)

    par(mfrow=c(3,1))
    name2<-c('PET','typeI','power')
    name21<-c('PET','type I error','power')
    for(i in 1:length(name2))
    {
      plot(1:dim(tmp0)[1],tmp0[,name2[i]],xlab='prior distribution: beta.a/beta.b or response rate (SD)',ylab=name21[i],main=name21[i],axes=F)
      box()
      axis(2)
      nn1<-apply(round(tmp0[,1:2],2),1,paste,collapse='/')
      axis(1,1:dim(tmp0)[1],dimnames(tmp0)[[1]])
      index1<-nn1==paste(beta.prior.list,collapse ='/')
      points((1:dim(tmp0)[1])[index1],tmp0[index1,name2[i]],col=2,cex=4,pch=2)
    }
    }
  })





  #---output report----
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste('bayesian_futility_by_predictive_probability', sep = '.',
            switch(input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'))
    },
    content = function(file) {
      out <- rmarkdown::render(input = 'bayesian_futility_by_predictive_probability_multi_stage.Rmd',
                    output_format =
                      switch(input$format,
                             PDF = rmarkdown::pdf_document(),
                             HTML = rmarkdown::html_document(),
                             Word = rmarkdown::word_document()
                      ),
                      params = list(set_title = input$project_title, set_author = input$author_input)
                    )
      file.rename(out, file)
    }
  )

}

# Run the application
shinyApp(ui = ui, server = server)

