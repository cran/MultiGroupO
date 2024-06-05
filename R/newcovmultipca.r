#' Function for the new covariance matrix in the multigroup PCA method
#'
#' Generates covariance matrix...
#'
#' @importFrom expm sqrtm
#' @importFrom stats cov
#' @import plsgenomics
#' @param x is a matrix with the data
#' @param cls is a vector of classes
#' @param A is a symmetric and positive definite matrix associated to inner product respect to the base of its vectorial space.
#'
#' @return \strong{return an grafics}.

#'
#' @examples
#' library(plsgenomics)
#' data(SRBCT)
#' mydata<-SRBCT$X
#' mydata<-mydata[1:50,1:20]
#' groups<-as.factor(SRBCT$Y)[1:50]
#' new.cov(x=mydata,cls=groups,A=diag(ncol(mydata)))
#' @export
new.cov<- function(x,cls,A)
{
  k <- nrow(x)
  n <-table(cls)
  n.cls <-length(n)

  Ss.p <- by(x,cls,function (y) cov(y)*((nrow(y)-1)*(k-nrow(y))))

   mean.x <- apply(x,2,mean)

  Bs <- by(x,cls,function (y)
  {mean.y <- apply(y,2,mean)
  (mean.y-mean.x)%*%t(mean.y-mean.x)*nrow(y)
  })

  B <-Reduce("+",Bs)
  B <- (1/k)*B

  Ss.all <-Reduce("+",Ss.p)
  sqA <- sqrtm(A)
  res<- sqA %*% (Ss.all+ ((k^2)*B))%*% sqA

  return(res)

}
