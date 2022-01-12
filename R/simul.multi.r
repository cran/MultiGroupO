#' Simulation function of quantitative multigroup data under a
#' multivariate normal distribution
#'
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom rlist list.remove
#' @param g An vector of the size of each group
#' @param mean1 An vector of the population means structure
#' @param d distance d for the structure of population means
#' @param n.var 2x1 dimension vector whose first component is the number of
#' random variables to simulate and the second component number
#' of noise variables to simulate
#' @param sds2 An vector of the variances to simulate for each group noise variables
#' @param corr An vector of the correlation to simulate for each group and noise variables
#'
#' @return \strong{return an grafics}
#'
#' @examples
#'fun.sim(g=c(20,20),mean1=2,d=0,sds2=c(1,1,1),corr=c(0.5,0.5,0),n.var=c(50,1))
#' @export
#'
fun.sim<-function(g,mean1,d,n.var,sds2,corr)
{
 # set.seed(seed)
  means<-vector();
  k=1
  M<-length(g)
  x1<-NULL
  liscovar<-vector("list", M)
  mu<-NULL
  lisx<-NULL

  for(j in 1:length(g))
  {
    if(j==1)
    {
      means[j]<-mean1
    }
    else{
      means[j]<-mean1+(j-1)*d
    }
  }
  means<-c(means,0)


  for(i in 1:(length(g)+1))
  {
    if(i<=length(g)){
      covar=matrix(NA, nrow=n.var[k], ncol=n.var[k], byrow=FALSE)
      covar[lower.tri(covar,diag=F)]<-corr[i]*sqrt(sds2[i]*sds2[i])
      covar[upper.tri(covar)]<-corr[i]*sqrt(sds2[i]*sds2[i])
      diag(covar)<-sqrt(sds2[i]*sds2[i])
      liscovar[[i]]<-covar


      if(d==0)
      {
        mu[[i]]<- c(means[i],0,rep(0,n.var[k]-2))

      }
      else{
      if(i==1){
        mu[[i]]<- c(means[i],0,rep(0,n.var[k]-2))
      }
      else if(i==2){
        mu[[i]]<- c(0,means[i],seq(2*means[i],4*means[i],length.out=n.var[k]-2))
      }
      else if(i==3){
        mu[[i]]<- c(1,0,2,-means[i],seq(6*means[i],8*means[i],length.out=n.var[k]-4))
      }
      else if(i==4){
        mu[[i]]<- c(1,0,2,0,3,means[i],rep(0,n.var[k]-6))
      }
      else if(i==5){
        mu[[i]]<- c(0,1,0,0,-means[i],rep(0,n.var[k]-5))
      }

      else {
        mu[[i]]<- c(2,0,2,0,means[i],rep(0,n.var[k]-5))
      }
}
      x1[[i]] <- rmvnorm(g[i], mean =mu[[i]],sigma=liscovar[[i]])

    }
    else{
      k=k+1
      covar=matrix(NA, nrow=n.var[k], ncol=n.var[k], byrow=FALSE)
      covar[lower.tri(covar,diag=F)]<-corr[i]*sqrt(sds2[i]*sds2[i])
      covar[upper.tri(covar)]<-corr[i]*sqrt(sds2[i]*sds2[i])
      diag(covar)<-sqrt(sds2[i]*sds2[i])
      x1[[i]] <- rmvnorm(sum(g), mean = rep(means[i],n.var[k]),sigma=covar)
    }
  }

  x2<-list.remove(x1,c(length(g)+1))
  grp <- factor(rep(LETTERS[1:length(g)],g))
  output <- do.call(rbind,lapply(x2,matrix,ncol=n.var[1],byrow=F))
  x3<-x1[[length(g)+1]]
  x4<-as.data.frame(cbind(output,x3,grp))

  for(i in 1:length(g)){
    lisx[[i]]<- subset(x4,x4$grp==i, select=-grp)
  }

  return(list(lisx=lisx,grp=grp,x=data.frame(x4),mu=mu,liscovar=liscovar))

}
