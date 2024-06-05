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
  means<-vector();
  M<-length(g)
  x1<-NULL
  liscovar<-vector("list", M)
  mu<-NULL
  lisx<-NULL
  mu2<-NULL

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
      covar=matrix(NA, nrow=n.var[1], ncol=n.var[1], byrow=FALSE)
      covar[lower.tri(covar,diag=F)]<-corr[i]*sqrt(sds2[i]*sds2[i])
      covar[upper.tri(covar)]<-corr[i]*sqrt(sds2[i]*sds2[i])
      diag(covar)<-sqrt(sds2[i]*sds2[i])
      liscovar[[i]]<-covar


      if(i==1){
        mu[[i]]<- seq(from =means[i],to=means[i]+2,length=n.var[1])
        x1[[i]] <- rmvnorm(g[i], mean =mu[[i]],sigma=liscovar[[i]])
      }
      else if(i==2){
        mu[[i]]<- seq(from =means[i],to=means[i]+2,length=n.var[1])
        x1[[i]] <- rmvnorm(g[i], mean =mu[[i]],sigma=liscovar[[i]])
      }
      else if(i==3){
        mu[[i]]<- seq(from =means[i],to=means[i]+2,length=n.var[1])
        x1[[i]] <- rmvnorm(g[i], mean =mu[[i]],sigma=liscovar[[i]])
      }
      else if(i==4){
        mu[[i]]<- seq(from =means[i],to=means[i]+2,length=n.var[1])
        x1[[i]] <- rmvnorm(g[i], mean =mu[[i]],sigma=liscovar[[i]])
      }
      else if(i==5){
        mu[[i]]<- seq(from =means[i],to=means[i]+2,length=n.var[1])
        x1[[i]] <- rmvnorm(g[i], mean =mu[[i]],sigma=liscovar[[i]])
      }

      #mu[[i]]<- rep(means[i],n.var[1])
      #x1[[i]] <- rmvnorm(g[i], mean =mu[[i]],sigma=liscovar[[i]])

    }
    else{

      covar=matrix(NA, nrow=n.var[2], ncol=n.var[2], byrow=FALSE)
      covar[lower.tri(covar,diag=F)]<-corr[i]*sqrt(sds2[i]*sds2[i])
      covar[upper.tri(covar)]<-corr[i]*sqrt(sds2[i]*sds2[i])
      diag(covar)<-sqrt(sds2[i]*sds2[i])
      mu2<-rep(means[i],n.var[2])
      x1[[i]] <- rmvnorm(sum(g), mean = mu2,sigma=covar)
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

  return(list(lisx=lisx,grp=as.factor(as.numeric(grp)),x=data.frame(x4),mu=mu,liscovar=liscovar,mu2=mu2))
}

