#' Performs a principal components analysis in the
#' given data matrix.
#' Show PCA graphical output.
#'
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot geom_point scale_color_brewer theme geom_density scale_fill_brewer aes xlab ylab
#' @importFrom gridExtra grid.arrange
#' @param datos is a matrix with the data
#' @param grupos is a vector of classes
#' @param Plot vector logic for grafic
#' @param center data set center by columns
#' @param scale  data set scaled by columns
#'
#' @return \strong{return an grafics}.
#'
#' @examples
#' library(plsgenomics)
#' data(SRBCT)
#' mydata<-SRBCT$X
#' mydata<-mydata[1:30,1:20]
#' groups<-as.factor(SRBCT$Y)[1:30]
#' pca(datos=mydata,grupos=groups,Plot=TRUE,center=TRUE,scale=FALSE)
#' @export
pca<-function(datos,grupos,Plot=TRUE,center=TRUE,scale=FALSE)
{
  PC1<-NULL
  PC2<-NULL
  PC3<-NULL
  pcs<-prcomp(datos,center =center, scale =scale)
  pcdatos<-as.data.frame(pcs$x)
  pcdatos<-cbind(pcdatos,grupos)
  group<-pcdatos$grupos

  percentage <- round((pcs$sdev^{2} / sum(pcs$sdev^{2}))*100, 2)
  percentage <- paste(colnames(pcdatos), "(", paste( as.character(percentage),
                                                     "%", ")", sep="") )

  if(Plot)
  {
    scatterPlot1 <- ggplot(pcdatos,aes(PC1, PC2, color=group)) +
      geom_point(size = 3) + xlab(percentage[1]) + ylab(percentage[2])+
      scale_color_brewer(palette="Accent")+
      theme(legend.position="bottom")

    scatterPlot2 <- ggplot(pcdatos,aes(PC1, PC3, color=group)) +
      geom_point(size = 3) + xlab(percentage[1]) + ylab(percentage[3])+
      scale_color_brewer(palette="Accent")+
      theme(legend.position="bottom")

    scatterPlot3 <- ggplot(pcdatos,aes(PC2, PC3, color=group)) +
      geom_point(size = 3) + xlab(percentage[2]) + ylab(percentage[3])+
      scale_color_brewer(palette="Accent")+
      theme(legend.position="bottom")

    xdensity <- ggplot(pcdatos, aes(PC1, fill=group)) +
      geom_density(alpha=.5) +
      scale_fill_brewer(palette="Accent")+
      theme(legend.position = "none")

    ydensity <- ggplot(pcdatos, aes(PC2, fill=group)) +
      geom_density(alpha=.5) +
      scale_fill_brewer(palette="Accent")+
      theme(legend.position = "none")

    zdensity <- ggplot(pcdatos, aes(PC3, fill=group)) +
      geom_density(alpha=.5) +
      scale_fill_brewer(palette="Accent")+
      theme(legend.position = "none")

    grid.arrange(xdensity,ydensity,zdensity, scatterPlot1,scatterPlot2,scatterPlot3,
                 ncol=3, nrow=2, widths=c(2,2,2), heights=c(2, 2.5))
  }
  # calculate explained variance
  prop_expl_var <- round(pcs$sdev^{2} / sum(pcs$sdev^{2}), 2)

  r<-list(loadings=pcs$rotation,variates=pcs$x,prop_expl_var=prop_expl_var)
  r
  class(r) = c("pca")
  return(r)
}
