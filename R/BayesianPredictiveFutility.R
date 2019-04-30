#’ A packeage for calculating BayesianPredictiveFutility
#’
#’ @docType package
#’ @name BayesianPredictiveFutility
NULL

#' Utility function to create list of lists
#'
#' @param lst, parent  list to add specified vector to
#' @param ..., vector to be added to list
#'
#' @return parent list with vector added as last element
#' @examples
#' newList = lappend(list(), c(1,2,3,4,5))
#' @export
lappend <- function(lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

#' Creates different combinations of out trial outcomes for a list stage sizes
#'
#' @param sizes an array of stage sizes
#' @param nreqCum an array with the required positive outcomes for all stages
#'
#' @return void, results are stored in a global variable
#' @export
#'
#' @examples
generateCombinationsRecursively <- function(sizes,nreqCum){

  for (stage in 1:length(sizes)) {
    #initiate an inner list to store by stage, accessed within recursive loop as a global variable
    combinations[[stage]] <<-  list()

    if (stage <=  length(nreqCum)) {
      #initiate value vector to store cumulative success counts in
      values = rep(0,stage)

      #set number of successes for last stage to the number of successes that needs to be exceeded
      values[stage] = nreqCum[stage]

      res <- combGen(stages = stage, currStage = 1,sizes, cumsum(sizes),nreqCum,values)
    }
  }

}

#' A function for generating outcome combinations for calculations using analytic form
#'
#' @param stages total number of stages to generate combinations for
#' @param currStage value for handling depth of recursion
#' @param sizes a vector of sample sizes for the stages
#' @param sizesCum  a vector of cumulative sample sizes for the stages
#' @param reqCum  a vector with the required cumulative counts for successful trial outcome
#' @param values  a list for recursively passing on combinations vector for storing in global list
#'
#' @return
#' @export
#'
#' @examples
combGen <- function(stages, currStage, sizes, sizesCum, reqCum, values) {

  #when last stage has been reached
  if (currStage == stages) {
    combinations[[stages]] <<- lappend(combinations[[stages]],values)
    return()
  }

  prevCount = values[currStage - 1]

  #conditions for minimum number of patients passing current stage
  # begin from maximum of
  # 1) count from previous stage
  # 2) required count from current stage + 1
  minSuccesses = max(prevCount,reqCum[currStage] + 1)

  #conditions for maximum number of patients passing current stage
  # end from minimum of
  # 1) cumulative number of patients at current stage
  # 2) required cumulative successes for final stage
  # 3) counts from previous stage plus number of patients in current stage
  maxSuccesses = min(sizesCum[currStage],reqCum[stages],prevCount+sizes[currStage])

  if (maxSuccesses >= minSuccesses)
    for (successes in minSuccesses:maxSuccesses) {
      values[currStage] = successes
      combGen(stages, currStage + 1,sizes, sizesCum, reqCum, values)
    }
}


#' Calculates probability for a sequence of outcomes fulfilling the criteria for analytic form power/sensitivity calculations
#'
#' @param passCombinations an array with outcomes for all stages
#' @param totals an array with the total number of patients in the study design
#' @param p probability to test sensitivity/power for
#'
#' @return probability for outcome combination
#' @export
#'
#' @examples
prob.fun <- function(passCombinationCum,stageSamples,p, INCREASE) {
  tmp_length <- length(passCombinationCum)

  # format from total counts to utilize counts per stage
  pcDiff = diff(passCombinationCum)
  pcDiff  = c(passCombinationCum[1],pcDiff)

  prod(dbinom(pcDiff[-tmp_length],stageSamples[-tmp_length],prob = p)) *
    pbinom(pcDiff[tmp_length],stageSamples[tmp_length],prob = p)
}

#' Calculates probabilities for power/sensitivity calculations
#'
#' @param pp probability to evaluate
#' @param nn.list array of stage sizes
#' @param passCombinations a list of lists with generated combinations for all stages
#'
#' @return a data object with probabilities and labels
#' @export
#'
#' @examples
analyticForm <- function(single_pp,nn.list,passCombinations, INCREASE) {

  #Probability of Early Termination
  PETs = numeric(0)

  for (stage in 1:length(nn.list)) {

    stageCombinations = passCombinations[[stage]]

    pet = 0

    # check that there is valid combinations by which a stage may fail
    if (length(stageCombinations) > 0) {
      pet = sum(sapply(stageCombinations,prob.fun,stageSamples = nn.list[1:stage],p = single_pp, INCREASE = INCREASE))
    }
    PETs = c(PETs,pet)
  }

  sumProb = sum(PETs)

  names(PETs) <- paste('PET.',1:length(PETs),'stage',sep = '')
  # todo confirm this bug has been corrected
  # Don't want the last PET anymore (not needed to Overall prob of stopping early)
  PETs <- c(PETs[-length(PETs)],1 - sumProb)
  names(PETs)[length(PETs)] <- 'prob.success'

  return(PETs)
}


#' Calculate beta binomial probabilities for increase scenario
#'
#' @param n1
#' @param n2
#' @param beta.a
#' @param beta.b
#' @param k
#'
#' @return
#' @export
#'
#' @examples
calcPbb <- function(n1,n2,beta.a,beta.b,k,INCREASE = TRUE){
  successes <- beta.a + k
  fails <- beta.b + (n1 - k)

  if (INCREASE) {
    return(1 - pbbinom(0:(n2 - 1),n2,successes,fails))
  }
  else{
    return(pbbinom(0:(n2),n2,successes,fails))
  }
}

#' Calculate beta binomial probabilities for reduction scenario
#'
#' @param n1
#' @param n2
#' @param beta.a
#' @param beta.b
#' @param k
#'
#' @return
#' @export
#'
#' @examples
calcPbbForReduction <- function(n1,n2,beta.a,beta.b,k){
  successes <- beta.a + k
  fails <- beta.b + (n1 - k)

  return(pbbinom(0:(n2),n2,successes,fails))
}

#' Calculate probability of successful trial based on current and remaining tested patients
#'
#' @param n1 number of tested patients
#' @param n2 number of remaining patients
#' @param betabinProbs
#' @param n.needed.for.greater.p0
#' @param increase boolean value specifying study type, TRUE = increase FALSE = reduce
#'
#' @return
#' @export
#'
#' @examples
getSuccessProbability <- function(n1,n2,betabinProbs,n.needed.for.greater.p0,increase = TRUE){
  tmp3 <- numeric(0)

  if (increase) {
    neededResponderRange <- max(0, n.needed.for.greater.p0 - n2):min(n1,(n.needed.for.greater.p0-1))
  }
  else{
    neededResponderRange <- min(n1,n2,n.needed.for.greater.p0):max(0)
  }

  for (i in neededResponderRange)
  {
    remaining <- min(n2, n.needed.for.greater.p0 - i)

    tmp <- betabinProbs[paste(i),paste(remaining),drop = F]
    tmp3 <- rbind(tmp3, c(i,remaining,as.numeric(dimnames(tmp)[[1]]),as.numeric(dimnames(tmp)[[2]]), as.vector(tmp)))
  }

  return(tmp3)
}

#' Simulation based sensitivity analysis
#'
#' @param pp
#' @param n1
#' @param n2
#' @param r1
#' @param r
#'
#' @return
#' @export
#'
#' @examples
sensAnalysis <- function(pp,n1,n2,r1,r)
{
  # r1: x<= r1|1st stage: stop trial at 1st stage
  # r: x<= r|at end of study: no success
  prob.stage1 <- pbinom(r1,n1,pp)
  # case which pass the 1st stage, but fail in the 2nd stage
  r.tmp <- seq(r1 + 1,min(n1,r))
  # prob(passing 1st stage, but fail in 2nd stage)
  prob.stage2 <- sum(dbinom(r.tmp,n1,pp) * pbinom(r - r.tmp,n2,pp))
  c(prob.stage1,1 - (prob.stage1 + prob.stage2))
}


#' Run simulations for sensitivity/power calculations
#'
#' @param senrates
#' @param sim.n
#' @param n.list
#' @param n.predictive.cutoff.list
#' @param n.needed.for.greater.p0
#' @param incrMode
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
simulateInterims <- function(senrates,sim.n,n.list,n.predictive.cutoff.list,n.needed.for.greater.p0,incrMode, seed = 42534253){

  stageSensitivities <- numeric(0)

  for (i in 1:length(senrates)) {
    p <- senrates[i]
    data.n <- as.matrix(n.list)
    dimnames(data.n)[[1]] <- paste('Interim',1:length(n.list))

    set.seed(seed)
    binomSim <- apply(data.n,1,function(x) rbinom(sim.n,x,p))

    binomSim.cum <- t(apply(binomSim[,-dim(binomSim)[2],drop = F],1,cumsum))
    tmp40 <- rbind(unlist(n.predictive.cutoff.list),binomSim.cum)
    if (incrMode) {
      stagesFailure <- cbind(apply(tmp40,2,function(x) x[-1] <= x[1]),apply(binomSim,1,sum) >= n.needed.for.greater.p0)
    }
    else{
      stagesFailure <- cbind(apply(tmp40,2,function(x) x[-1] >= x[1]),apply(binomSim,1,sum) <= n.needed.for.greater.p0)
    }


    #------------------------------------------------
    tmp50 <- apply(stagesFailure[,-dim(stagesFailure)[2],drop = F],1,function(x) {index1 <- cumsum(x)>0;ifelse(sum(index1) == 0,0,(1:length(index1))[index1][1])})
    tmp50 <- factor(tmp50,levels = 0:(length(n.list) - 1))
    tmp501 <- table(tmp50)/length(tmp50)
    #prob.stop <- tmp501[-grep('0',names(tmp501))]
    prob.stop <- tmp501[names(tmp501) != '0'] #--remove the success one---
    #------------------------------------------------
    prob.stop.name <- paste('prob.stop.interim.',1:(dim(stagesFailure)[2] - 1),sep = '')
    # ---status of no stopping in each simulation---
    predFailureList <- apply(stagesFailure[,-dim(stagesFailure)[2],drop = F],1,sum) == 0
    pwr <- mean(predFailureList*stagesFailure[,dim(stagesFailure)[2]])
    stageSensitivities <- rbind(stageSensitivities, c(prob.stop,pwr))
    #---1st to last 2 columns are prob of early stopping for each interim analysis
    #----last column is power or type I error: prob(pass all stages and # of event <= the cutoff)
  }

  dimnames(stageSensitivities) <- list(senrates,c(prob.stop.name,'prob.above.threshold.of.reject.Ho'))
  rejectProbs <- t(apply(stageSensitivities,1,function(x) c(sum(x[-length(x)]),x[length(x)])))
  dimnames(rejectProbs) <- list(senrates,c('prob.stop.overall','prob.above.threshold.of.reject.Ho'))

  return(list(stageSensitivities = stageSensitivities, sensitivity.data = rejectProbs))
}

#todo: utilize if performance is an issue
#' More efficient implementation for calculating beta binomial distributions, overcoming performance drawback with r-style pushing to arrays
#'
#' @param n1
#' @param n2
#' @param beta.a
#' @param beta.b
#'
#' @return
#' @export
#'
#' @examples
calcBetabinProbabilitiesInPlace <- function(n1,n2,beta.a,beta.b) {
  nRow <- n1
  nCol <- n2

  r1 <- as.list(calcPbb(n1,n2,beta.a,beta.b,0))
  d <- data.frame(r1,
                  stringsAsFactors = FALSE)

  if (nRow > 1) {
    d <- d[rep.int(1,nRow + 1),]
      # lose data.frame class for a while
      # changes what S3 methods implement
      # assignment.
      d <- as.list(d)
    for (i in seq.int(1,nRow,1)) {
      ri <- calcPbb(n1,n2,beta.a,beta.b,i)
      for (j in seq_len(nCol)) {
        print(j)
        d[[j]][i + 1] <- ri[[j]]
      }
    }
  }
    d <- data.frame(d,stringsAsFactors = FALSE)

  dimnames(d) <- list(0:n1,1:n2)

  return(d)
}

#' Calculate beta binomial probabilities for both increase and reduce scenario
#'
#' @param n1
#' @param n2
#' @param beta.a
#' @param beta.b
#' @param increase
#'
#' @return
#' @export
#'
#' @examples
calcBetabinProbabilities <- function(n1,n2,beta.a,beta.b,INCREASE = TRUE) {
  return(bbProb(n1,n2,beta.a,beta.b,INCREASE))
}


#' Calculate beta bionomial probabilities
#'
#' @param n1
#' @param n2
#' @param beta.a
#' @param beta.b
#' @param INCREASE
#'
#' @return
#' @export
#'
#' @examples
bbProb <- function(n1,n2,beta.a,beta.b, INCREASE = TRUE){

  tmp <- numeric(0)

  testRange <- 0:n1

  for (k in testRange)
  {
    tmp <- rbind(tmp,calcPbb(n1,n2,beta.a,beta.b,k,INCREASE))
  }

  dimnames(tmp) <- list(testRange,ifelse(INCREASE,1,0):n2)

  return(tmp)
}

bbProbIncr <- function(n1,n2,beta.a,beta.b){

  tmp <- numeric(0)

  testRange <- 0:n1

  for (k in testRange)
  {
    tmp <- rbind(tmp,calcPbb(n1,n2,beta.a,beta.b,k))
  }

  dimnames(tmp) <- list(testRange,1:n2)

  return(tmp)
}

bbProbRed <- function(n1,n2,beta.a,beta.b) {

  tmp <- numeric(0)

  # Is it necessary to go through all combinations of n1|n2? Seems like it is only used for GUI :TODO
  # Make it possible to loop through list of stage sample sizes :TODO

  testRange <- 0:n1
  for (k in testRange)
  {
    #pbbinom(x,n2,beta.a+k,beta.b+(n1-k)): prob(X<= x|beta(a,b,k))
    # so prob(X>x) =  1- pbbinom(x,n2,beta.a+k,beta.b+(n1-k))
    # therefore, it becomes prob(X>= 1), prob(X>= 2),...,prob(X>= n2)
    # for x = 0,1,...,n2-1
    #tmp2[k] <- 1-pbbinom(0:(n2-1),n2,beta.a+k,beta.b+(n1-k))
    #list[k] <- 1-pbbinom(0:(n2-1),n2,beta.a+k,beta.b+(n1-k))
    tmp <- rbind(tmp,calcPbbForReduction(n1,n2,beta.a,beta.b,k))
  }

  dimnames(tmp) <- list(testRange,0:n2)

  return(tmp)
}

#' Calculates number of positive outcomes from total number of patients that is required for the full study to be successful
#'
#' @param ns array of patients per stage
#' @param p.target target probability for successful study
#' @param beta.a beta a
#' @param beta.b beta b
#' @param cutoff.n.for.greater.p0 cutoff.n.for.greater.p0
#' @param increase boolean for specifying if study evaluation is aimed capture significant increase or decrease. TRUE = increase FALSE = reduce
#'
#' @return numeric value of outcomes that needs to be exceded for positive outcome of study
#' @export
#'
#' @examples
nNeededFromTotal <- function(ns, p.target, beta.a, beta.b, cutoff.n.for.greater.p0,increase){

  total_n <- sum(ns)

  responderCountRange <- 0:(total_n)

  betaProbs <- (1 - pbeta(p.target, beta.a + responderCountRange, beta.b + (total_n) - responderCountRange))

  if (increase) {
    nNeeded <- responderCountRange[betaProbs > cutoff.n.for.greater.p0][1]
  }
  else{
    nNeeded <- tail(responderCountRange[betaProbs < cutoff.n.for.greater.p0],1)
  }

  return(nNeeded)
}

nNeeded <- function(n1,n2, p.target, beta.a, beta.b, cutoff.n.for.greater.p0){
  responderCountRange <- 0:(n1 + n2)
  betaProbs <- (1 - pbeta(p.target, beta.a+responderCountRange, beta.b + (n1 + n2) - responderCountRange))

  return(responderCountRange[betaProbs > cutoff.n.for.greater.p0][1])
}


#' Main function for calculating statistics for the interim study design
#'
#' @param ns
#' @param p.target
#' @param cutoff.n.for.greater.p0
#' @param predictive.cutoff
#' @param beta.a
#' @param beta.b
#' @param outcome.tmp
#' @param arm.name
#' @param p.h1
#' @param analysis.type
#' @param sim.n
#' @param seed
#'
#' @return
#'
#' @examples
#'
#' @import shiny
#' @import knitr
#' @import extraDistr
#' @import rmarkdown
#' @export
#'

evaluateInterim <-
  function(ns,p.target = 0.2, cutoff.n.for.greater.p0 = 0.95, predictive.cutoff = 0.05, beta.a = 1, beta.b = 1, outcome.tmp = 'response', arm.name = 'B', p.h1 = .4, analysis.type = c('Analytical', 'Simulation'), sim.n = 10000, seed = 42534253)
  {
    analysis.type <- match.arg(analysis.type)
    if (p.h1 > p.target)
      INCREASE <- TRUE else
        INCREASE <- FALSE

    failureCountLengthList <- list()

    successProbList <- list()


    beta.a <- beta.a + 10^(-100)
    beta.b <- beta.b + 10^(-100)

    n.needed.for.greater.p0 <- nNeededFromTotal(ns, p.target, beta.a, beta.b, cutoff.n.for.greater.p0,INCREASE)
if(!(n.needed.for.greater.p0>0))
{
  warning('All PASS, no need for stopping')
  #return('no_stopping')
  return('All PASS, no need for stopping. You may want to lower beta prior values, increase  sample size or the threshold of the posterior probability for efficacy.  ')
}
    nrOfStages <- length(ns)

    #------------------------------------------------------------------------
    #section for determining the required cumulative number of successes for n-1 stages
    #------------------------------------------------------------------------
    n.predictive.cutoff.list <- list()
    betaBinProbList <- list()

    for (i in 1:(nrOfStages - 1)) {

      testedPatients <- sum(ns[1:i])
      remainingPatients <- sum(ns[(i + 1):nrOfStages])

      localBetabinProbs <- bbProb(testedPatients,remainingPatients,beta.a,beta.b,INCREASE)

      betaBinProbList[[i]] <- localBetabinProbs

      successProb <- getSuccessProbability(testedPatients,remainingPatients,localBetabinProbs,n.needed.for.greater.p0, INCREASE)
      successProbList[[i]] = successProb

      failureCounts <- successProb[successProb[,5] < predictive.cutoff,1]
      failureCountLength <- length(failureCounts)

      if (failureCountLength > 0) {
        failureCountLengthList[[i]] <- failureCountLength
        n.predictive.cutoff.list[[i]] <- failureCounts[failureCountLength]
      } else {
        failureCountLengthList[[i]] <- (-99)

        if (successProb[1,1] == 0) {
          n.predictive.cutoff.list[[i]] <- failureCounts[failureCountLength]
        } else {
          n.predictive.cutoff.list[[i]] <- successProb[1,1] - 1
        }
      }



    if (length(n.predictive.cutoff.list[[1]]) ==  0 || n.predictive.cutoff.list[[1]] ==  0)
    {
    #a1<-n.predictive.cutoff.list
    #index99<-sapply(a1,function(x) x==0|length(x)==0)
    #a11<-paste(paste(unlist(a1),' in stage',1:length(a1)),collapse=', ')
    #a12<-paste('No stopping at',paste(paste('stage',(1:length(index99))[index99]),collapse=', '))
    warning('With such a low unfavorable p value no events needed for stopping in stage ' ,i, ', hence this method is not appropriate')
    #a13<-paste(a12,'. Stopping boundary is ',a11, sep='')
    #return(a13)
    return(paste('No stopping in the 1st stage. Please increase the sample size in the 1st stage  and rerun the program.'))
    }
    }

    #------------------------------------------------------------------------

    pp <- sort(unique(round(c(p.target - c(0.1,0.05,.15,.2,.25,.3), p.target, p.target + c(0.1,0.05,.15,.2,.25,.3),p.h1), 5)))
    senRates <- pp[(pp > 0) & (pp < 1)]


    #------------------------------------------------------------------------
    #Calculate power and sensitivity using analytic form or simulation
    #------------------------------------------------------------------------
      required = c(unlist(n.predictive.cutoff.list),n.needed.for.greater.p0 - 1)

      combinations <<-  list()
      generateCombinationsRecursively(ns,required)

      #if less than specified number of stages use analytic form
      if (analysis.type ==  'Analytical')
      {
        #---1st to last 2 columns are prob of early stopping for each interim analysis
        #----last column is power or type I error: prob(pass all stages and # of event <= the curoff)
        #todo: [important] make sure that the values below are correct
        # stageSensitivities <- t(apply(t(senRates),2,analyticForm,nn.list = cumsum(ns), combinations))
        stageSensitivities <- t(apply(t(senRates),2,analyticForm,nn.list = ns, passCombinations = combinations, INCREASE = INCREASE))


        prob.stop.name <- paste('prob.stop.interim.',1:(dim(stageSensitivities)[2] - 1),sep = '')
        dimnames(stageSensitivities) <- list(senRates,c(prob.stop.name,'prob.above.threshold.of.reject.Ho'))
        sensProbs <- t(apply(stageSensitivities,1,function(x) c(sum(x[-length(x)]),x[length(x)])))
        dimnames(sensProbs) <- list(senRates,c('prob.stop.overall','prob.above.threshold.of.reject.Ho'))

     } else{
       # Simulation Method
        simulateInterims_output <- simulateInterims(senRates,sim.n,ns,n.predictive.cutoff.list, n.needed.for.greater.p0,INCREASE, seed = seed)
        stageSensitivities <- simulateInterims_output$stageSensitivities
        sensProbs <- simulateInterims_output$sensitivity.data
      }


      list(n.needed.for.greater.p0 = n.needed.for.greater.p0,
           n.predictive.cutoff = unlist(n.predictive.cutoff.list),
           all.data = betaBinProbList,
           predictive.prob.data = successProbList,
           sensitivity.data = sensProbs,
           sensitivity.data.ind = stageSensitivities,
           kk2.len.list = failureCountLengthList
      )

  }
