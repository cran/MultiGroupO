#' Performs a Single Dimension Exploration (SDE) analysis in the
#' given multigroup data matrix.
#' Show SDE graphical output.
#'
#' @param mydata data set
#' @param groups is a vector of classes
#' @param plt grafics
#' @import gplots
#' @import ggrepel
#' @import ggplot2
#' @import qgraph
#' @import mgm
#' @return \strong{return an grafics} .
#'
#'
#' @examples
#'sim.list2<-fun.sim(g=c(20,50,10),mean1=0.5,d=0,sds2=c(1,1,1,1),corr=c(0.1,0.5,0.5,0),
#'n.var=c(20,20))
#'datos2 <- as.data.frame(sim.list2$x)
#'datos2<-subset(datos2,select=-grp)
#'grupos <- sim.list2$grp
#'grupos<-factor(grupos,labels=c(1,2,3))
#'sde.method(mydata=datos2,groups=grupos,plt=FALSE)
#' @export
sde.method<-function(mydata,groups,plt=FALSE){
  X<-list();Xn<-list()
  n<-list();nn<-list();
  u<-list();un<-list();
  A<-list();z<-list();
  B<-list();BB<-list();
  gcor<-list();
    for(i in 1:length(levels(groups))){
      X[[i]]<-as.matrix(mydata[groups==i,])
      Xn[[i]]<-as.matrix(mydata[groups!=i,])

      n[i]<-nrow(as.matrix(X[[i]]))
      nn[i]<-nrow(as.matrix(Xn[[i]]))

      u[[i]]<-as.matrix(rep(1,n[i]),n[i],1)
      un[[i]]<-as.matrix(rep(1,nn[i]),nn[i],1)

      A[[i]]<- X[[i]]-u[[i]]%*%t(un[[i]])%*%Xn[[i]]/nn[[i]]
      z[[i]]<-t(A[[i]])%*%u[[i]]

      B[[i]]<-z[[i]]%*%t(z[[i]])


      if(plt==TRUE){
      #      BB<-Reduce("+",B)
      BB<-as.matrix(as.data.frame(B[[i]]))

      res<-eigen(BB)
      eigenvals<-res$values
      eigenvecs<-res$vectors
      Y<-as.matrix(mydata)%*%eigenvecs
      as.matrix(cbind(Y[,i],as.matrix(mydata)))

      datos2<-as.matrix(cbind(Y[,i],as.matrix(mydata)))
      nc<-dim(datos2)[2]
      type<-rep("g",nc)
      level<-rep(1,nc)

      group_list <- list("grupo1v"=c(1),
                         "grupo2v"=c(2,4,5,6),
                         "grupo3v" = c(7),
                         "grupo3v"=c(3))
      # define nice colors
      group_cols <- c("#E35959","#8FC45A","#4B71B3","#E8ED61")


      fit <- mgm(datos2, type,level, lamda.sel="EBIC", k=2)
      gcor[[i]]<-fit$pairwise$wadj
      qgraph(fit$pairwise$wadj,
      edge.color = fit$pairwise$edgecolor,
      layout = "spring"
      )
      }


    }

  r<-list(B)
  r
  class(r) = c("sde")
  return(r)
 #G <- do.call(rbind, lapply(gcor, matrix, ncol = ncol(as.matrix(mydata)),byrow = F))
 #   G<-as.matrix(G)
 }
