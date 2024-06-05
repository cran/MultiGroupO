## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
  )


## ----setup--------------------------------------------------------------------
library(plsgenomics)
library(MultiGroupO)
data(SRBCT)
 mydata<-SRBCT$X
 mydata<-mydata[,1:5]
 groups<-as.factor(SRBCT$Y)
 pca(datos=mydata,grupos=groups,Plot=TRUE,center=TRUE,scale=TRUE)
      

## ----setup1-------------------------------------------------------------------
library(plsgenomics)
data(SRBCT)
mydata<-SRBCT$X
mydata<-mydata[,1:5]
groups<-as.factor(SRBCT$Y)
mat.to.diag1<-new.cov(x=mydata,cls=groups,A=diag(ncol(mydata)))
mgpca(mat.to.diag=mat.to.diag1,mat.x=as.matrix(mydata),cls=groups,Plot=TRUE,ncomp=2,center = TRUE,scale = TRUE)

## ----setup2-------------------------------------------------------------------
library(plsgenomics)
data(SRBCT)
mydata<-SRBCT$X
mydata<-mydata[,1:5]
groups<-as.factor(SRBCT$Y)
mydata<-split(as.data.frame(mydata),groups)
mdr(group=groups,data.x=mydata,c=2)

