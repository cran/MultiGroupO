% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BIplot.r
\name{BIplot}
\alias{BIplot}
\title{biplot methods}
\usage{
BIplot(
  variates,
  loadings,
  prop_expl_var,
  comp = c(1, 2),
  group = NULL,
  rownamevar = T,
  rownameload = T
)
}
\arguments{
\item{variates}{is the size of groups}

\item{loadings}{is a vector of classes}

\item{prop_expl_var}{data set}

\item{comp}{component numeric}

\item{group}{is a vector of groups}

\item{rownamevar}{is a logical vector where TRUE is the
label of the observations, if is FALSE, is index.}

\item{rownameload}{is a logical vector where TRUE is the
label of the vectors of loadings, if is FALSE, is index.}
}
\value{
\strong{return an grafics} .
}
\description{
biplot methods
}
\examples{
library(datasets)
obj<-pca(datos=iris[,-5],grupos=iris[,5],Plot=FALSE,center=TRUE,scale=TRUE)
BIplot(variates=obj$variates,loadings=obj$loadings,
      prop_expl_var=obj$prop_expl_var,comp=c(1,2),
      group=factor(as.numeric(iris[,5])),rownamevar=FALSE,rownameload=TRUE)

}
