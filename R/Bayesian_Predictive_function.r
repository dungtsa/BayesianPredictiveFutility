#' @description The function generates report and statistical plan for futility interim analysis by the Bayesian predictive design for single arm early  phase II clinical trial
#' @param n1: sample size in the 1st stage
#' @param n2: sample size in the 2nd stage
#' @param p.target: response rate in null hypothesis
#' @param cutoff.n.for.greater.p0: the threshold for posterior probability. That is, prob(response rate>p.target|data) >cutoff.n.for.greater.p0
#' @param predictive.cutoff: the cutoff for the predictive probability
#' @param beta.a: beta prior for the 1st parameter
#' @param beta.b: beta prior for the 2nd parameter
#' @param outcome.tmp:  label for the outcome
#' @param arm.name:  label for the experimental arm
#' @param p.h1: response rate in alternative hypothesis
#' @param plot.status: indictor whether to generate plot
#' @name  Bayesian.predictive.futility.fun
#' @title Bayesian predictive design
#' @examples
#' #--- 15 subjects for the 1st and 2nd stages with H0:prob(response rate>0.2|data) >0.95 and the cutoff for the predictive probability=0.05
#' Bayesian.predictive.futility.fun(n1=15,n2=15,p.target=0.20,cutoff.n.for.greater.p0=0.95,predictive.cutoff=0.05,beta.a=1,beta.b=1,outcome.tmp='response',arm.name='A')
#' @export
#'

#-------
Bayesian.predictive.futility.fun<-
  function(n1,n2,p.target=0.2,cutoff.n.for.greater.p0=0.95,predictive.cutoff=0.05,beta.a=1,beta.b=1,outcome.tmp='response',arm.name='B',p.h1=.4,plot.status=F)
  {
    require(extraDistr)
    beta.a<-beta.a+10^(-100)
    beta.b<-beta.b+10^(-100)
    tmp2<-numeric(0)
    for(k in 0:n1)
    {
      #pbbinom(x,n2,beta.a+k,beta.b+(n1-k)): prob(X<=x|beta(a,b,k))
      # so prob(X>x)= 1- pbbinom(x,n2,beta.a+k,beta.b+(n1-k))
      # therefore, it becomes prob(X>=1), prob(X>=2),...,prob(X>=n2)
      # for x=0,1,...,n2-1
      tmp2<-rbind(tmp2,1-pbbinom(0:(n2-1),n2,beta.a+k,beta.b+(n1-k)))
    }
    dimnames(tmp2)<-list(0:n1,1:n2)

    kk<-0:(n1+n2)
    kk1<-kk[(1-pbeta(p.target,beta.a+kk,beta.b+(n1+n2)-kk))>cutoff.n.for.greater.p0]
    n.needed.for.greater.p0<-kk1[1]
    #if(n.needed.for.greater.p0==0) stop('The rule always reject H0')
    tmp3<-numeric(0)
    for(i in max(0,n.needed.for.greater.p0-n2):min(n1,(n.needed.for.greater.p0-1)))
    {
      #print(c(i,n.needed.for.greater.p0-i))
      tmp<-tmp2[paste(i),paste(min(n2,n.needed.for.greater.p0-i)),drop=F]
      tmp3<-rbind(tmp3, c(i,min(n2,n.needed.for.greater.p0-i),as.numeric(dimnames(tmp)[[1]]),as.numeric(dimnames(tmp)[[2]]), as.vector(tmp)))
    }

    kk2<-tmp3[tmp3[,5]<predictive.cutoff,1]
    kk2.len<-length(kk2)

    pp<-sort(unique(c(p.target-c(0.1,0.05),p.target,p.target+c(0.1,0.05,.15,.2,.25,.3),p.h1)))
    pp<-pp[(pp>0)&(pp<1)]


    if(length(kk2)==0)
    {
      sen.null<-matrix(NA,length(pp),2)
      dimnames(sen.null)<-list(pp,c('prob.stop','prob.above.threshold.of.reject.Ho'))
      list(n.needed.for.greater.p0=n.needed.for.greater.p0,
           n.predictive.cutoff=NA,
           all.data=tmp2,
           predictive.prob.data=tmp3,
           sensitivity.data=sen.null)
    } else
    {
      n.predictive.cutoff<-kk2[length(kk2)]
      if(plot.status)
      {

        plot(tmp3[,5],type='h',axes=F,xlab=paste('number of patients with ', outcome.tmp,' in the 1st ',n1,' patients',sep=''),ylab='predictive probability',lwd=3,cex.lab=1.5)
        box();axis(2,cex.axis=1.2)
        axis(2,predictive.cutoff,cex.axis=1.2)
        axis(1,1:dim(tmp3)[1],tmp3[,3],cex.axis=1.2)
        axis(3,1:dim(tmp3)[1],tmp3[,4],cex.axis=1)
        abline(h=predictive.cutoff,col=2,lty=2)
        mtext(paste('number of patients with ',outcome.tmp,' needed in the furture remaining ', n2,' patients',sep=''),side=3,line=1.8)
        title(paste('Interim Analysis of Futility for Arm ',arm.name,'\n\n',sep=''))
        arrows(kk2.len, 1, kk2.len, y1 = 0.1,col=3,lwd=3)
        text(1:dim(tmp3)[1],tmp3[,5],round(tmp3[,5],2),col=2,cex=2)
      }


      fun1<-function(pp,n1,n2,r1,r)
      {
        # r1: x<=r1|1st stage: stop trial at 1st stage
        # r: x<=r|at end of study: no success
        porb.stage1<-pbinom(r1,n1,pp)
        # case which pass the 1st stage, but fail in the 2nd stage
        r.tmp<-seq(r1+1,min(n1,r))
        # prob(passing 1st stage, but fail in 2nd stage)
        porb.stage2<-sum(dbinom(r.tmp,n1,pp)*pbinom(r-r.tmp,n2,pp))
        c(porb.stage1,1-(porb.stage1+porb.stage2))
      }
      tmp5<-t(apply(t(pp),2,fun1,n1=n1,n2=n2,r1=n.predictive.cutoff,r=n.needed.for.greater.p0-1))
      dimnames(tmp5)<-list(pp,c('prob.stop','prob.above.threshold.of.reject.Ho'))


      if(plot.status)
      {
        par(mfrow=c(2,1))
        aa1<-barplot(tmp5[,1],xlab=paste('Rate of ',outcome.tmp,sep=''),ylab='Probability of early stopping the arm',main='Sensitivity Analysis \n Probability of early stopping the treatment arm',lwd=3,cex.lab=1.5,ylim=c(0,1),cex.main=2,cex.lab=1.5,cex.axis=1.5,cex=1.5)
        text(aa1,tmp5[,1],round(tmp5[,1],2),col=2,cex=2)

        aa1<-barplot(tmp5[,2],xlab=paste('Rate of ',outcome.tmp,sep=''),ylab=paste('Probability of at least ', n.needed.for.greater.p0,' patients with ', outcome.tmp,sep=''),main=paste('Probability of at least ', n.needed.for.greater.p0,' patients with ',outcome.tmp,sep=''),lwd=3,cex.lab=1.5,ylim=c(0,1),cex.main=2,cex.lab=1.5,cex.axis=1.5,cex=1.5)
        text(aa1,tmp5[,2]+.1,round(tmp5[,2],2),col=2,cex=2)
        par(mfrow=c(1,1))
      }
      list(n.needed.for.greater.p0=n.needed.for.greater.p0,
           n.predictive.cutoff=n.predictive.cutoff,
           all.data=tmp2,
           predictive.prob.data=tmp3,
           sensitivity.data=tmp5)
    }
  }

