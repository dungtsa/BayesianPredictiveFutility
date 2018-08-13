#' @description The function generates report and statistical plan for futility interim analysis by the Bayesian predictive design for single arm early  phase II clinical trial (p1>p0)
#' @param n.list: sample sizes for each stage (vector length = 2)
#' @param p.target: response rate in null hypothesis
#' @param cutoff.n.for.greater.p0: the threshold for posterior probability. That is, prob(response rate>p.target|data) >cutoff.n.for.greater.p0
#' @param predictive.cutoff: the cutoff for the predictive probability
#' @param beta.a: beta prior for the 1st parameter
#' @param beta.b: beta prior for the 2nd parameter
#' @param outcome.tmp:  label for the outcome
#' @param arm.name:  label for the experimental arm
#' @param p.h1: response rate in alternative hypothesis
#' @param sim.n: number of simulations
#' @param plot.status: indictor whether to generate plot
#' @name  Bayesian.predictive.futility.increase.fun
#' @title Bayesian predictive design
#' @examples
#' #--- 15 subjects for the 1st and 2nd stages with H0:prob(response rate>0.2|data) >0.95 and the cutoff for the predictive probability=0.05
#' Bayesian.predictive.futility.increase.fun(n1=15,n2=15,p.target=0.20,cutoff.n.for.greater.p0=0.95,predictive.cutoff=0.05,beta.a=1,beta.b=1,outcome.tmp='response',arm.name='A')
#' @export


Bayesian.predictive.futility.increase.fun<-
  function(n.list,p.target=0.3,cutoff.n.for.greater.p0=0.9,predictive.cutoff=0.05,beta.a=1,beta.b=1,outcome.tmp='response',arm.name='B',p.h1=.5,sim.n=10000,plot.status=F)
  {
    if (p.target > p.h1) stop('"p.target" can not be greater than "p.h1". Considering using "Bayesian.predictive.futility.reduction.fun"')
    #pbbinom(x,n2,beta.a+k,beta.b+(n1-k)): prob(X<=x|beta(a,b,k))
    # so prob(X>x)= 1- pbbinom(x,n2,beta.a+k,beta.b+(n1-k))
    # therefore, it becomes prob(X>=1), prob(X>=2),...,prob(X>=n2)
    # for x=0,1,...,n2-1
    beta.a<-beta.a+10^(-100)
    beta.b<-beta.b+10^(-100)
    
    tmp2.list<-tmp3.list<-n.predictive.cutoff.list<-kk2.len.list<-list()
    n.sum<-sum(n.list)
    kk<-0:n.sum
    kk1<-kk[(1-pbeta(p.target,beta.a+kk,beta.b+n.sum-kk))>=cutoff.n.for.greater.p0]
    n.needed.for.greater.p0<-kk1[1]
    
    for(j in 1:(length(n.list)-1))
    {
      n1<-sum(n.list[1:j])
      n2<-sum(n.list[(j+1):length(n.list)])
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
      tmp2.list[[j]]<-tmp2
      tmp3<-numeric(0)
      for(i in max(0,n.needed.for.greater.p0-n2):min(n1,(n.needed.for.greater.p0-1)))
      {
        tmp<-tmp2[paste(i),paste(min(n2,n.needed.for.greater.p0-i)),drop=F]
        tmp3<-rbind(tmp3, c(i,min(n2,n.needed.for.greater.p0-i),as.numeric(dimnames(tmp)[[1]]),as.numeric(dimnames(tmp)[[2]]), as.vector(tmp)))
      }
      tmp3.list[[j]]<-tmp3
      kk2<-tmp3[tmp3[,5]<predictive.cutoff,1]
      if((length(kk2)>0))
      {
        kk2.len.list[[j]]<- length(kk2)
        n.predictive.cutoff.list[[j]]<-kk2[length(kk2)]
      }else 
      {if(tmp3[1,1]==0)
      {
        kk2.len.list[[j]]<-(-99)
        n.predictive.cutoff.list[[j]]<-(-99)
      }else 
      {
        kk2.len.list[[j]]<-(-99)
        n.predictive.cutoff.list[[j]]<-tmp3[1,1]-1 
      }}
      
    }
    
    if(plot.status)
    {
      for(j in 1:length(tmp3.list)) 
      {
        tmp3<-tmp3.list[[j]]
        n1<-sum(n.list[1:j])
        n2<-sum(n.list[(j+1):length(n.list)])
        kk2.len<-kk2.len.list[[j]]
        n.predictive.cutoff<-n.predictive.cutoff.list[[j]]
        stage.index<-ifelse(j==1,' 1st ',ifelse(j==2,' 2nd ',ifelse(j==3,' 3rd ',paste(' ',j,'th ',sep=''))))
        plot(tmp3[,5],type='h',axes=F,xlab=paste('number of patients with ', outcome.tmp,' in the',stage.index,' interim analysis with ',n1,' patients',sep=''),ylab='predictive probability',lwd=3,cex.lab=1.5,ylim=c(0,1))
        box();axis(2,cex.axis=1.2)
        axis(2,predictive.cutoff,cex.axis=1.2)
        axis(1,1:dim(tmp3)[1],tmp3[,3],cex.axis=1.2)
        axis(3,1:dim(tmp3)[1],tmp3[,4],cex.axis=1)
        abline(h=predictive.cutoff,col=2,lty=2)
        mtext(paste('largest number of patients with ',outcome.tmp,' needed in the furture remaining ', n2,' patients',sep=''),side=3,line=1.8)
        title(paste(stage.index,' Interim Analysis for Futility \n\n',sep=''))
        arrows(kk2.len, 1, kk2.len, y1 = 0.1,col=3,lwd=3)
        text(1:dim(tmp3)[1],tmp3[,5],round(tmp3[,5],2),col=2,cex=2)
      }
    }
    
    #pp<-sort(unique(c(p.target-c(0.1,0.05),p.target,p.target+c(0.1,0.05,.15,.2),p.h1)))
    pp<-sort(unique(c(p.target-c(0.1,0.05,.15,.2,.25,.3),p.target,p.target+c(0.1,0.05,.15,.2,.25,.3),p.h1)))
    pp<-pp[(pp>0)&(pp<1)]
    
    if(1>2)
    {
      #---for 2 stages only----
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
    }
    if(length(n.list)<5)
    {#--for 2-4 stages or more: use analytic form----
      fun1<-function(pp,nn.list,ri.list)
      {
        n.cum<-cumsum(nn.list)
        r.total<-ri.list[length(ri.list)]
        n.total<-n.cum[length(n.cum)]
        # r1: x<=r1|1st stage: stop trial at 1st stage
        # r: x<=r|at end of study: no success
        porb.1st.stop<-pbinom(ri.list[1],n.cum[1],pp)
        ans<-porb.1st.stop
        #--stop at 2nd stage---
        r.tmp<-seq(ri.list[1]+1,min(n.cum[1],ri.list[2]))
        porb.2nd.stop<-sum(dbinom(r.tmp,n.cum[1],pp)*pbinom(ri.list[2]-r.tmp,n.cum[2]-n.cum[1],pp))
        ans<-c(ans,porb.2nd.stop)
        #  print('----2---')
        if((length(ri.list)>=3)&(length(ri.list)<5))
        {
          #--stop at 3rd stage---
          #  print('----3rd---')
          a1=(ri.list[2]+1):ri.list[3]
          
          tmp3<-numeric(0)
          for(i in 1:length(a1))
          {
            tmp1<-(ri.list[1]+1):min(c(a1[i],n.cum[1]))
            tmp2<-sum(dbinom(tmp1,n.cum[1],pp)*dbinom(a1[i]-tmp1,n.cum[2]-n.cum[1],pp))*
              pbinom(ri.list[3]-a1[i],n.cum[3]-n.cum[2],pp)
            tmp3<-c(tmp3,tmp2)
          }
          porb.3rd.stop<-sum(tmp3)
          ans<-c(ans,porb.3rd.stop)
        }
        
        if(length(ri.list)==4)
        {
          #--stop at 4th stage---
          #  print('----4th---')
          a1=(ri.list[3]+1):ri.list[4]
          
          tmp4<-numeric(0)
          for(i in 1:length(a1))
          {
            tmp1<-(ri.list[1]+1):min(c(a1[i],n.cum[1]))
            tmp4.0<-numeric(0)
            for(j in 1:length(tmp1))
            {
              tmp2<-max(0,(ri.list[2]+1)-tmp1[j]):min(c(a1[i]-tmp1[j],n.cum[2]-n.cum[1]))
              tmp21<-sum(dbinom(tmp1[j],n.cum[1],pp)*dbinom(tmp2,n.cum[2]-n.cum[1],pp)*dbinom(a1[i]-tmp1[j]-tmp2,n.cum[3]-n.cum[2],pp))*
                pbinom(ri.list[4]-a1[i],n.cum[4]-n.cum[3],pp)
              tmp4.0<-c(tmp4.0,tmp21)
            }
            tmp4<-c(tmp4,sum(tmp4.0))
          }
          porb.4th.stop<-sum(tmp4)
          ans<-c(ans,porb.4th.stop)
        }
        
        names(ans)<-paste('PET.',1:length(ans),'stage',sep='')
        ans<-c(ans,1-sum(ans))
        names(ans)[length(ans)]<-'prob.success'
        ans
      }
      
      #---1st to last 2 columns are prob of early stopping for each interim analysis
      #----last column is power or type I error: prob(pass all stages and # of event <=the curoff)
      tmp5.ind<-t(apply(t(pp),2,fun1,nn.list=n.list,ri.list=c(unlist(n.predictive.cutoff.list),n.needed.for.greater.p0-1)))
      prob.stop.name<-paste('prob.stop.interim.',1:(dim(tmp5.ind)[2]-1),sep='')
      dimnames(tmp5.ind)<-list(pp,c(prob.stop.name,'prob.above.threshold.of.reject.Ho'))
      tmp5<-t(apply(tmp5.ind,1,function(x) c(sum(x[-(length(x)-c(1,0))]),x[length(x)])))
      dimnames(tmp5)<-list(pp,c('prob.early.stop.overall','prob.above.threshold.of.reject.Ho'))
    } else
    {  #--for 5 stages or more: use simulation----
      tmp5.ind<-numeric(0)
      for(i in 1:length(pp))
      {
        p1<-pp[i]
        data.n<-as.matrix(n.list)
        dimnames(data.n)[[1]]<-paste('Interim',1:length(n.list))
        tmp4<-apply(data.n,1,function(x) rbinom(sim.n,x,p1))
        tmp4.cum<-t(apply(tmp4[,-dim(tmp4)[2],drop=F],1,cumsum))
        tmp40<-rbind(unlist(n.predictive.cutoff.list),tmp4.cum)
        tmp41<-cbind(apply(tmp40,2,function(x) x[-1]<=x[1]),apply(tmp4,1,sum)>=n.needed.for.greater.p0) 
        tmp50<-apply(tmp41[,-dim(tmp41)[2],drop=F],1,function(x) {index1<-cumsum(x)>0;ifelse(sum(index1)==0,0,(1:length(index1))[index1][1])})
        tmp50<-factor(tmp50,level=0:(length(n.list)-1))
        tmp501<-table(tmp50)/length(tmp50)
        tmp501<-tmp501[names(tmp501)!='0'] #--remove the success one---
        prob.stop.name<-paste('prob.stop.interim.',1:(dim(tmp41)[2]-1),sep='')
        # ---status of no stopping in each simulation---
        tmp411<-apply(tmp41[,-dim(tmp41)[2],drop=F],1,sum)==0
        tmp5.ind<-rbind(tmp5.ind, c(tmp501,mean(tmp411*tmp41[,dim(tmp41)[2]])))
        #---1st to last 2 columns are prob of early stopping for each interim analysis
        #----last column is power or type I error: prob(pass all stages and # of event <=the curoff)
      }
      dimnames(tmp5.ind)<-list(pp,c(prob.stop.name,'prob.above.threshold.of.reject.Ho'))
      
      tmp5<-t(apply(tmp5.ind,1,function(x) c(sum(x[-length(x)]),x[length(x)])))
      dimnames(tmp5)<-list(pp,c('prob.early.stop.overall','prob.above.threshold.of.reject.Ho'))
    }
    
    if(plot.status)
    {
      par(mfrow=c(2,1))
      aa1<-barplot(tmp5[,1],xlab=paste('Rate of ',outcome.tmp,sep=''),ylab='Probability of early stopping the treatment overall',main='Sensitivity Analysis \n Overall probability of early stopping the treatment',lwd=3,cex.lab=1.5,ylim=c(0,1),cex.main=2,cex.lab=1.5,cex.axis=1.5,cex=1.5)
      text(aa1,tmp5[,1],round(tmp5[,1],2),col=2,cex=2)
      
      aa1<-barplot(tmp5[,2],xlab=paste('Rate of ',outcome.tmp,sep=''),ylab=paste('Probability of at most ', n.needed.for.greater.p0,' patients with ', outcome.tmp,sep=''),main=paste('Probability of at least ', n.needed.for.greater.p0,' patients with ',outcome.tmp,sep=''),lwd=3,cex.lab=1.5,ylim=c(0,1),cex.main=2,cex.lab=1.5,cex.axis=1.5,cex=1.5)
      text(aa1,tmp5[,2]+.1,round(tmp5[,2],2),col=2,cex=2)
      par(mfrow=c(1,1))
      dim.n<-dim(tmp5.ind)[2]
      if(dim.n>2)
      {
        par(mfrow=c(2,ceiling(dim.n-1)/2))
        for(j in 1:(dim.n-1))
        {
          aa1<-barplot(tmp5.ind[,j],xlab=paste('Rate of ',outcome.tmp,sep=''),ylab=paste('Probability of early stopping the treatment at interim',j),main=paste('Interim Analysis: ',j,'\n Probability of early stopping the treatment',sep=''),lwd=3,cex.lab=1.5,ylim=c(0,1),cex.main=2,cex.lab=1.5,cex.axis=1.5,cex=1.5)
          text(aa1,tmp5.ind[,j],round(tmp5.ind[,j],2),col=2,cex=2)
        }
        
      }
    }
    list(n.needed.for.greater.p0=n.needed.for.greater.p0,
         n.predictive.cutoff=unlist(n.predictive.cutoff.list),
         all.data=tmp2.list,
         predictive.prob.data=tmp3.list,
         sensitivity.data=tmp5,
         sensitivity.data.ind=tmp5.ind,
         kk2.len.list=kk2.len.list)
  }




#' @description The function generates report and statistical plan for futility interim analysis by the Bayesian predictive design for single arm early  phase II clinical trial (p1<p0)
#' @param n.list: sample sizes for each stage (vector length = 2)
#' @param p.target: response rate in null hypothesis
#' @param cutoff.n.for.greater.p0: the threshold for posterior probability. That is, prob(response rate>p.target|data) >cutoff.n.for.greater.p0
#' @param predictive.cutoff: the cutoff for the predictive probability
#' @param beta.a: beta prior for the 1st parameter
#' @param beta.b: beta prior for the 2nd parameter
#' @param outcome.tmp:  label for the outcome
#' @param arm.name:  label for the experimental arm
#' @param p.h1: response rate in alternative hypothesis
#' @param sim.n: number of simulations
#' @param plot.status: indictor whether to generate plot
#' @name  Bayesian.predictive.futility.reduction.fun
#' @title Bayesian predictive design
#' @examples
#' #--- 15 subjects for the 1st and 2nd stages with H0:prob(response rate>0.2|data) >0.95 and the cutoff for the predictive probability=0.05
#' Bayesian.predictive.futility.reduction.fun(n1=15,n2=15,p.target=0.20,cutoff.n.for.greater.p0=0.95,predictive.cutoff=0.05,beta.a=1,beta.b=1,outcome.tmp='response',arm.name='A')
#' @export


Bayesian.predictive.futility.reduction.fun<-
  function(n.list,p.target=0.27,cutoff.n.for.greater.p0=0.9,predictive.cutoff=0.05,beta.a=1,beta.b=1,outcome.tmp='event',arm.name='B',p.h1=.17,sim.n=10000,plot.status=F)
  {
    if (p.target < p.h1) stop('"p.target" can not be less than "p.h1". Considering using "Bayesian.predictive.futility.increase.fun"')
    
    beta.a<-beta.a+10^(-100)
    beta.b<-beta.b+10^(-100)
    
    tmp2.list<-tmp3.list<-n.predictive.cutoff.list<-kk2.len.list<-list()
    n.sum<-sum(n.list)
    kk<-0:n.sum
    kk1<-kk[pbeta(p.target,beta.a+kk,beta.b+n.sum-kk)>=cutoff.n.for.greater.p0]
    n.needed.for.greater.p0<-kk1[length(kk1)]
    
    for(j in 1:(length(n.list)-1))
    {
      n1<-sum(n.list[1:j])
      n2<-sum(n.list[(j+1):length(n.list)])
      tmp2<-numeric(0)
      for(k in 0:n1)
      {
        #pbbinom(x,n2,beta.a+k,beta.b+(n1-k)): prob(X<=x|beta(a,b,k))
        # so prob(X>x)= 1- pbbinom(x,n2,beta.a+k,beta.b+(n1-k))
        # therefore, it becomes prob(X>=1), prob(X>=2),...,prob(X>=n2)
        # for x=0,1,...,n2-1
        tmp2<-rbind(tmp2,pbbinom(0:(n2),n2,beta.a+k,beta.b+(n1-k)))
      }
      dimnames(tmp2)<-list(0:n1,0:n2)
      tmp2.list[[j]]<-tmp2
      
      tmp3<-numeric(0)
      for(i in min(n1,n2,n.needed.for.greater.p0):max(0))
      {
        tmp<-tmp2[paste(i),paste(min(n2,n.needed.for.greater.p0-i)),drop=F]
        tmp3<-rbind(tmp3, c(i,min(n2,n.needed.for.greater.p0-i),as.numeric(dimnames(tmp)[[1]]),as.numeric(dimnames(tmp)[[2]]), as.vector(tmp)))
      }
      tmp3.list[[j]]<-tmp3
      kk2<-tmp3[tmp3[,5]<predictive.cutoff,1]
      kk2.len<-length(kk2)
      kk2.len.list[[j]]<-kk2.len
      n.predictive.cutoff<-kk2[length(kk2)]
      n.predictive.cutoff.list[[j]]<-n.predictive.cutoff
    }
    
    if(plot.status)
    {
      for(j in 1:length(tmp3.list)) 
      {
        tmp3<-tmp3.list[[j]]
        n1<-sum(n.list[1:j])
        n2<-sum(n.list[(j+1):length(n.list)])
        kk2.len<-kk2.len.list[[j]]
        n.predictive.cutoff<-n.predictive.cutoff.list[[j]]
        stage.index<-ifelse(j==1,' 1st ',ifelse(j==2,' 2nd ',ifelse(j==3,' 3rd ',paste(' ',j,'th ',sep=''))))
        plot(tmp3[,5],type='h',axes=F,xlab=paste('number of patients with ', outcome.tmp,' in the',stage.index,' interim analysis with ',n1,' patients',sep=''),ylab='predictive probability',lwd=3,cex.lab=1.5)
        box();axis(2,cex.axis=1.2)
        axis(2,predictive.cutoff,cex.axis=1.2)
        axis(1,1:dim(tmp3)[1],tmp3[,3],cex.axis=1.2)
        axis(3,1:dim(tmp3)[1],tmp3[,4],cex.axis=1)
        abline(h=predictive.cutoff,col=2,lty=2)
        mtext(paste('largest number of patients with ',outcome.tmp,' needed in the furture remaining ', n2,' patients',sep=''),side=3,line=1.8)
        title(paste(stage.index,' Interim Analysis of for Futility \n\n',sep=''))
        arrows(kk2.len, 1, kk2.len, y1 = 0.1,col=3,lwd=3)
        text(1:dim(tmp3)[1],tmp3[,5],round(tmp3[,5],2),col=2,cex=2)
      }
    }
    
    #pp<-sort(unique(c(p.target-c(0.1,0.05),p.target,p.target+c(0.1,0.05,.15,.2),p.h1)))
    pp<-sort(unique(c(p.target-c(0.1,0.05,.15,.2,.25,.3),p.target,p.target+c(0.1,0.05,.15,.2,.25,.3),p.h1)))
    
    pp<-pp[(pp>0)&(pp<1)]
    tmp5.ind<-numeric(0)
    for(i in 1:length(pp))
    {
      p1<-pp[i]
      data.n<-as.matrix(n.list)
      dimnames(data.n)[[1]]<-paste('Interim',1:length(n.list))
      tmp4<-apply(data.n,1,function(x) rbinom(sim.n,x,p1))
      tmp4.cum<-t(apply(tmp4[,-dim(tmp4)[2],drop=F],1,cumsum))
      tmp40<-rbind(unlist(n.predictive.cutoff.list),tmp4.cum)
      tmp41<-cbind(apply(tmp40,2,function(x) x[-1]>=x[1]),apply(tmp4,1,sum)<=n.needed.for.greater.p0) 
      tmp50<-apply(tmp41[,-dim(tmp41)[2],drop=F],1,function(x) {index1<-cumsum(x)>0;ifelse(sum(index1)==0,0,(1:length(index1))[index1][1])})
      tmp50<-factor(tmp50,level=0:(length(n.list)-1))
      tmp501<-table(tmp50)/length(tmp50)
      tmp501<-tmp501[-grep('0',names(tmp501))]
      prob.stop.name<-paste('prob.stop.interim.',1:(dim(tmp41)[2]-1),sep='')
      # ---status of no stopping in each simulation---
      tmp411<-apply(tmp41[,-dim(tmp41)[2],drop=F],1,sum)==0
      tmp5.ind<-rbind(tmp5.ind, c(tmp501,mean(tmp411*tmp41[,dim(tmp41)[2]])))
      #---1st to last 2 columns are prob of early stopping for each interim analysis
      #----last column is power or type I error: prob(pass all stages and # of event <=the curoff)
    }
    
    dimnames(tmp5.ind)<-list(pp,c(prob.stop.name,'prob.above.threshold.of.reject.Ho'))
    tmp5<-t(apply(tmp5.ind,1,function(x) c(sum(x[-length(x)]),x[length(x)])))
    dimnames(tmp5)<-list(pp,c('prob.stop.overall','prob.above.threshold.of.reject.Ho'))
    
    if(plot.status)
    {
      par(mfrow=c(2,1))
      aa1<-barplot(tmp5[,1],xlab=paste('Rate of ',outcome.tmp,sep=''),ylab='Probability of early stopping the treatment overall',main='Sensitivity Analysis \n Overall probability of early stopping the treatment',lwd=3,cex.lab=1.5,ylim=c(0,1),cex.main=2,cex.lab=1.5,cex.axis=1.5,cex=1.5)
      text(aa1,tmp5[,1],round(tmp5[,1],2),col=2,cex=2)
      
      aa1<-barplot(tmp5[,2],xlab=paste('Rate of ',outcome.tmp,sep=''),ylab=paste('Probability of at most ', n.needed.for.greater.p0,' patients with ', outcome.tmp,sep=''),main=paste('Probability of at most ', n.needed.for.greater.p0,' patients with ',outcome.tmp,sep=''),lwd=3,cex.lab=1.5,ylim=c(0,1),cex.main=2,cex.lab=1.5,cex.axis=1.5,cex=1.5)
      text(aa1,tmp5[,2]+.1,round(tmp5[,2],2),col=2,cex=2)
      par(mfrow=c(1,1))
      dim.n<-dim(tmp5.ind)[2]
      if(dim.n>2)
      {
        par(mfrow=c(2,ceiling(dim.n-1)/2))
        for(j in 1:(dim.n-1))
        {
          aa1<-barplot(tmp5.ind[,j],xlab=paste('Rate of ',outcome.tmp,sep=''),ylab=paste('Probability of early stopping the treatment at interim',j),main=paste('Interim Analysis: ',j,'\n Probability of early stopping the treatment',sep=''),lwd=3,cex.lab=1.5,ylim=c(0,1),cex.main=2,cex.lab=1.5,cex.axis=1.5,cex=1.5)
          text(aa1,tmp5.ind[,j],round(tmp5.ind[,j],2),col=2,cex=2)
        }
        
      }
    }
    list(n.needed.for.greater.p0=n.needed.for.greater.p0,
         n.predictive.cutoff=unlist(n.predictive.cutoff.list),
         all.data=tmp2.list,
         predictive.prob.data=tmp3.list,
         sensitivity.data=tmp5,
         sensitivity.data.ind=tmp5.ind,
         kk2.len.list=kk2.len.list)
  }

