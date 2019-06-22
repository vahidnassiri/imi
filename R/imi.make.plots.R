imi.make.plots<- function(sim.out){
  Min.mvn=  apply(rbind(apply(sim.out$sim.mean$sim.cs3.mvn.mean[,1:5],2,min),
                        apply(sim.out$sim.mean$sim.cs7.mvn.mean[,1:5],2,min),
                        apply(sim.out$sim.mean$sim.ar13.mvn.mean[,1:5],2,min),
                        apply(sim.out$sim.mean$sim.ar13.mvn.mean[,1:5],2,min),
                        apply(sim.out$sim.mean$sim.log3.mean[,1:5],2,min),
                        apply(sim.out$sim.mean$sim.log7.mean[,1:5],2,min)),2,min)
  
  Max.mvn=  apply(rbind(apply(sim.out$sim.mean$sim.cs3.mvn.mean[,1:5],2,max),
                        apply(sim.out$sim.mean$sim.cs7.mvn.mean[,1:5],2,max),
                        apply(sim.out$sim.mean$sim.ar13.mvn.mean[,1:5],2,max),
                        apply(sim.out$sim.mean$sim.ar13.mvn.mean[,1:5],2,max),
                        apply(sim.out$sim.mean$sim.log7.mean[,1:5],2,max)),2,max)
  
  Main=c('Euclidean','Inf-norm','Mahalanobis-within','Mahalanobis-between',
         'Mahalanabos-comb')
  # add.mvn=paste(folder.address,paste(paste('\\mvn',scenario,sep=''),'.pdf',sep=''),sep='')
  # 
  # pdf(paste(paste('mvn',scenario,sep=''),'.pdf',sep=''))
  
  par(mfrow=c(2,3))
  for (i in 1:5){
    plot(sim.out$sim.mean$sim.cs3.mvn.mean[,i],ylim=c(Min.mvn[i],Max.mvn[i]),ylab='Distance (log-scale)',xlab='M',
         main=Main[i],log='y')
    #grid(nx=20)
    points(sim.out$sim.mean$sim.cs7.mvn.mean[,i],pch=2,col=2)
    points(sim.out$sim.mean$sim.ar13.mvn.mean[,i],pch=3,col=3)
    points(sim.out$sim.mean$sim.ar17.mvn.mean[,i],pch=4,col=4)
    points(sim.out$sim.mean$sim.log3.mean[,i],pch=5,col=6)
    points(sim.out$sim.mean$sim.log7.mean[,i],pch=6,col=8)
  }
  i=1
  plot(sim.out$sim.mean$sim.cs3.mvn.mean[,i],ylim=c(Min.mvn[i],Max.mvn[i]),ylab='',xlab='',
       main='',log='y',type="n",axes=FALSE,ann=FALSE)
  legend('center',c('CS-10%','CS-70%','AR1-10%','AR1-70%','Logreg-10%','Logreg-70%'),
         col=c(1,2,3,4,6,8),pch=c(1,2,3,4,5,6),cex=1.6,bty = "n")
  
  
  # dev.off()
  
  # Now plot for t-test
  
  dev.new()
  
  
  Min.t=  apply(rbind(apply(sim.out$sim.t$sim.cs3.t.mean,2,min),
                      apply(sim.out$sim.t$sim.cs7.t.mean,2,min),
                      apply(sim.out$sim.t$sim.ar13.t.mean,2,min),
                      apply(sim.out$sim.t$sim.ar13.t.mean,2,min)),2,min)
  
  Max.t=  apply(rbind(apply(sim.out$sim.t$sim.cs3.t.mean,2,max),
                      apply(sim.out$sim.t$sim.cs7.t.mean,2,max),
                      apply(sim.out$sim.t$sim.ar13.t.mean,2,max),
                      apply(sim.out$sim.t$sim.ar13.t.mean,2,max)),2,max)
  
  Main=c(expression(paste(mu[0],'=0')),expression(paste(mu[0],'=10')),
         expression(paste(mu[1]-mu[2],'=0')),
         expression(paste(mu[1]-mu[2],'=10')))
  # pdf(paste(paste('ttest',scenario,sep=''),'.pdf',sep=''))
  par(mfrow=c(2,2))
  for (i in 1:4){
    plot(sim.out$sim.t$sim.cs3.t.mean[,i],ylim=c(Min.t[i],Max.t[i]),ylab='Distance (log-scale)',xlab='M',
         main=Main[i],log='y')
    #grid(nx=20)
    points(sim.out$sim.t$sim.cs7.t.mean[,i],pch=2,col=2)
    points(sim.out$sim.t$sim.ar13.t.mean[,i],pch=3,col=3)
    points(sim.out$sim.t$sim.ar17.t.mean[,i],pch=4,col=4)
    if (i==3){
      legend('topright',c('CS-10%','CS-70%','AR1-10%','AR1-70%'),
             col=c(1,2,3,4),pch=c(1,2,3,4))
    }
    
    
  }
  # dev.off()
  
}
