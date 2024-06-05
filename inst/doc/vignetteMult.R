## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  fig.width=6,
  fig.height=5,
  collapse = TRUE,
  comment = "#"
)


## ----setup--------------------------------------------------------------------
library(MultiGroupO)
sim.list<-fun.sim(g=c(30,30,30),mean1=2,d=0,sds2=c(1,1,1,1),corr=c(0.5,0.5,0.5,0),
                      n.var=c(20,50))
    datos1 <- as.data.frame(sim.list$x)
    datos1<-subset(datos1,select=-grp)
    grupos <- sim.list$grp
     pca(datos1,grupos,Plot=TRUE,center=TRUE,scale=FALSE)
     mat.to.diag1<-new.cov(datos1,cls=grupos,A=diag(ncol(datos1)))
      mgpca(mat.to.diag=mat.to.diag1,mat.x=as.matrix(datos1),cls=grupos,Plot=T,ncomp=2,center = TRUE,scale = TRUE)

## ----setup2-------------------------------------------------------------------
library(MultiGroupO)
sim.list<-fun.sim(g=c(30,30,30),mean1=2,d=0,sds2=c(1,1,1,1),corr=c(0.5,0.5,0.5,0),
                  n.var=c(20,50))
datos1 <- as.data.frame(sim.list$x)
datos1<-subset(datos1,select=-grp)
grupos <- sim.list$grp
mdr(group=grupos,data.x=sim.list$`lisx`,c=2)

