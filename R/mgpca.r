#' Performs a Multigroup PCA analysis in the
#' given multigroup data matrix.
#' Show mgpca graphical output.
#'
#'
#' @importFrom ggplot2 ggplot geom_point scale_color_brewer theme geom_density scale_fill_brewer aes xlab ylab labs
#' @importFrom gridExtra grid.arrange
#' @importFrom stats dist
#' @import plsgenomics
#' @param mat.to.diag is a matrix with the data
#' @param mat.x is a vector of classes
#' @param cls group
#' @param Plot grafics output of mgpca
#' @param ncomp number of component
#' @param center is a logical vector where TRUE is center (whether the variables should be
#' shifted to be zero centered), if is FALSE, is original data.
#' @param scale is a logical vector where TRUE is scale (indicating whether
#' the variables should be scaled), if is FALSE, is original data.
#'
#'
#'
#' @return \strong{If simplify == TRUE} class values.
#'         \strong{If simplify == FALSE}, the result is a list of length
#'         \code{nsim} data.tables.
#'
#' @examples
#' library(plsgenomics)
#' data(SRBCT)
#' mydata<-SRBCT$X
#' mydata<-mydata[1:50,1:5]
#' groups<-as.factor(SRBCT$Y)[1:50]
#' mat.to.diag1<-new.cov(x=mydata,cls=groups,A=diag(ncol(mydata)))
#' mgpca(mat.to.diag=mat.to.diag1,mat.x=as.matrix(mydata),
#' cls=groups,Plot=TRUE,ncomp=2,center = TRUE,scale = TRUE)
#' @export
mgpca<-function(mat.to.diag,mat.x,cls,Plot=TRUE,ncomp=2,center = TRUE,scale = TRUE)
{
  V1<-NULL
  V2<-NULL
  V3<-NULL
  mat.x<-scale(mat.x,center = center,scale = scale)
  s <- svd(mat.to.diag)
  u<-s$u
  PCA.points <- mat.x%*%u
  PCA.1.points <- mat.x%*%s$u[,1]
  PCA.2.points <- mat.x%*%s$u[,2]


  loadings = u
  values<-(s$d)
  group<-cls
  pcmgdatos<-as.data.frame(PCA.points)
  pcmgdatos<-cbind(pcmgdatos,group)


  # Calculo de sumas de interdistancias al cuadrado global y entre grupos.
  interd_global_x <- sum(as.matrix(dist(PCA.points)^2))

  interd_total <- as.matrix(dist(PCA.points)^2)
  sum_interd_total  <-sum(interd_total)

  #forma 1
  sum.within <- 0
  for (i in levels(cls)){
    sum.within <- sum.within + sum(interd_total[cls==i,cls==i])
  }

  sum_interd_between <- sum_interd_total  - sum.within


  pctg.1 <- sum(as.matrix(dist(PCA.1.points)^2))/sum_interd_total
  pctg.2 <- sum(as.matrix(dist(PCA.2.points)^2))/sum_interd_total



  percentage1 <- round(c(pctg.1,pctg.2)*100,2)
  percentage <- paste( "(", paste( as.character(percentage1), "%", ")", sep="") )

  if(Plot){
    scatterPlot1 <- ggplot(pcmgdatos,aes(V1,V2, color=group)) +
      geom_point(size = 3) + labs(x = paste("PC1",percentage[1]),
                                  y=paste("PC2",percentage[2]))+
      scale_color_brewer(palette="Set2")


      scatterPlot2 <- ggplot(pcmgdatos,aes(V1, V3, color=group)) +
      geom_point(size = 3) + labs(x = paste("PC1",percentage[1]),
                                  y=paste("PC3",percentage[3]))+
      scale_color_brewer(palette="Set2")


    scatterPlot3 <- ggplot(pcmgdatos,aes(V2, V3, color=group)) +
      geom_point(size = 3) + labs(x = paste("PC2",percentage[2]),
                                  y=paste("PC3",percentage[3]))+
      scale_color_brewer(palette="Set2")


    xdensity <- ggplot(pcmgdatos, aes(V1, fill=group)) +
      geom_density(alpha=.5) +   labs(x = "PC1")+
      scale_fill_brewer(palette="Set2")



    ydensity <- ggplot(pcmgdatos, aes(V2, fill=group)) +
      geom_density(alpha=.5) +  labs(x = "PC2")+
      scale_fill_brewer(palette="Set2")



    zdensity <- ggplot(pcmgdatos, aes(V3, fill=group)) +
      geom_density(alpha=.5) +  labs(x = "PC3")+
      scale_fill_brewer(palette="Set2")

    nt <- theme(legend.position='hidden')
    grid_arrange_shared_legend(xdensity+nt,  zdensity+nt,ydensity+nt,scatterPlot1, scatterPlot2+nt,scatterPlot3+nt,  nrow=2,ncol=3)

  }

  r<-list(variates=PCA.points,loadings=loadings,prop_expl_var=percentage1)
  r
  class(r) = c("mgpca")
  return(r)

}
