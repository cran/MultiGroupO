#' Performs a Multigroup Dimensionality Reduction (MDR) analysis in the
#' given multigroup data matrix.
#' Show MDR graphical output.
#'
#'
#'
#' @importFrom ggplot2 ggplot geom_point scale_color_brewer theme geom_density scale_fill_brewer aes xlab ylab labs geom_blank
#' @importFrom gridExtra grid.arrange
#' @param group is a vector of classes
#' @param data.x quantitative data set
#' @param c component numeric
#' @param Plot grafics output of MDR
#'
#' @return \strong{return an grafics} .
#'
#'
#' @examples
#' sim.list<-fun.sim(g=c(50,50,50),mean1=2,d=0,sds2=c(1,1,1,1),
#' corr=c(0.5,0.5,0.5,0),n.var=c(30,30))
#'
#' mdr(group=as.factor(sim.list$grp),
#' data.x=sim.list$`lisx`,c=2)
#'
#' @export
mdr<-function(group,data.x,c,Plot=T)
{
  g<-  as.vector(summary(group))
  if(is.list(data.x))
  {
    mg<-NULL
    ln1<-NULL
    A<-NULL
    z<-NULL
    C<-NULL
    CC<-NULL
    ys<-NULL
    Dim1<-NULL
    Dim2<-NULL
    groups<-NULL
    if(length(g)==2){
      for(l in 1:length(g)){
        mg[[l]]<-as.matrix(data.x[[l]])
        ln1[[l]]<-matrix(rep(1,nrow(mg[[l]])),ncol=1)
      }

      A<-mg[[2]]-1/nrow(mg[[1]])*ln1[[2]]%*%t(ln1[[1]])%*%mg[[1]]
      z<-t(ln1[[2]])%*%A
      C<-t(z)%*%z

      v<-eigen(C)$vectors[,1:c]
      loadings<-v

      ys1<-mg[[1]]%*%v

      y<-ys1
      datapcs<-cbind(y,group)
      colnames(datapcs)<-c("Dim1","Dim2","group")
      group<-as.factor(datapcs$group)
      percentage2 <- round((eigen(C)$values / sum(eigen(C)$values))*100, 2)
      percentage1 <- round((eigen(C)$values / sum(eigen(C)$values))*100, 2)
      percentage <- paste(colnames(datapcs), "(", paste( as.character(percentage1),
                                                         "%", ")", sep="") )
      group<-factor(datapcs$group,labels=c(1,2))
      groups=group
      datapcs<-cbind(datapcs,groups)
      if(Plot){
        scatterPlot1 <- ggplot(datapcs,aes(Dim1, Dim2, color=groups)) +
          geom_point(size = 3) + xlab(percentage[1]) + ylab(percentage[2])+
          scale_color_brewer(palette="Set2",labels = as.factor(1:length(g)))+
          theme(legend.position="bottom")



        xdensity <- ggplot(datapcs, aes(Dim1, fill=groups)) +
          geom_density(alpha=.5) +
          scale_fill_brewer(palette="Set2")+
          theme(legend.position = "none")

        ydensity <- ggplot(datapcs, aes(Dim2, fill=groups)) +
          geom_density(alpha=.5) +
          scale_fill_brewer(palette="Set2")+
          theme(legend.position = "none")

        blankPlot <- ggplot()+geom_blank(aes(1,1)) +
          cowplot::theme_nothing()

        grid<-grid.arrange(xdensity, blankPlot, scatterPlot1, ydensity,
                           ncol=2, nrow=2, widths=c(4, 3), heights=c(2.5, 4))
      }

      r<-list(variates=y,loadings=loadings,prop_expl_dim=percentage2)
      r
      class(r) = c("mdr")
      return(r)
    }


    else if(length(g)>2){

      s=1
      for(l in 1:length(g)){
        mg[[l]]<-as.matrix(data.x[[l]])
        ln1[[l]]<-matrix(rep(1,nrow(mg[[l]])),ncol=1)
      }

      for(w in 1:length(g)){
        for(r in 1:length(g)){
          if(w>r){
            A[[s]]<-mg[[w]]-1/nrow(mg[[r]])*ln1[[w]]%*%t(ln1[[r]])%*%mg[[r]]
            s=s+1
          }
        }
      }

      W=1
      for(w in 1:length(g)){
        for(r in 1:length(g)){
          if(w>r){
            z[[W]]<-t(ln1[[w]])%*%A[[W]]
            W=W+1
          }
        }
      }
      for(q in 1:length(g)){
        C[[q]]<-t(z[[q]])%*%z[[q]]
      }

      CC<-Reduce("+",C)
      v<-eigen(CC)$vectors[,1:c]
      loadings<-v

      for(p in 1:length(g)){
        ys[[p]]<-mg[[p]]%*%v
      }

      y <- as.data.frame(do.call(rbind,lapply(ys,matrix,ncol=c,byrow=F)))

      datapcs<-cbind(y,group)
      colnames(datapcs)<-c("Dim1","Dim2","group")
      group<-as.factor(datapcs$group)
      percentage2 <- round(eigen(CC)$values / sum(eigen(CC)$values)*100, 2)
      percentage1 <- round(eigen(CC)$values / sum(eigen(CC)$values)*100, 2)
      percentage <- paste(colnames(datapcs), "(", paste( as.character(percentage1),
                                                         "%", ")", sep="") )
      if(Plot){
        scatterPlot1 <- ggplot(datapcs,aes(Dim1, Dim2, color=group)) +
          geom_point(size = 3) + xlab(percentage[1]) + ylab(percentage[2])+
          scale_color_brewer(palette="Accent",labels = as.factor(1:length(g)))+
          theme(legend.position="bottom")



        xdensity <- ggplot(datapcs, aes(Dim1, fill=group)) +
          geom_density(alpha=.5) +
          scale_fill_brewer(palette="Set2")+
          theme(legend.position = "none")

        ydensity <- ggplot(datapcs, aes(Dim2, fill=group)) +
          geom_density(alpha=.5) +
          scale_fill_brewer(palette="Set2")+
          theme(legend.position = "none")

        blankPlot <- ggplot()+geom_blank(aes(1,1)) +
          cowplot::theme_nothing()

        grid<-grid.arrange(xdensity, blankPlot, scatterPlot1, ydensity,
                           ncol=2, nrow=2, widths=c(4, 3), heights=c(2.5, 4))
      }

      r<-list(variates=y,loadings=loadings,prop_expl_dim=percentage2)
      r
      class(r) = c("mdr")
      return(r)
    }
  }

  else{

    g<-  as.vector(summary(group))
    pcmgdatos<-as.data.frame(cbind(data.x,group))
    lisx<-NULL
    for(i in 1:length(g)){
      lisx[[i]]<- subset(pcmgdatos,pcmgdatos$group==i, select=-group)
    }
    data.x1<-lisx
    {
      mg<-NULL
      ln1<-NULL
      A<-NULL
      z<-NULL
      C<-NULL
      CC<-NULL
      ys<-NULL
      Dim1<-NULL
      Dim2<-NULL
      # A1<-NULL
      # z1<-NULL
      # C1<-NULL
      # C2<-NULL
      grupos<-NULL
      if(length(g)==2){
        for(l in 1:length(g)){
          mg[[l]]<-as.matrix(data.x1[[l]])
          ln1[[l]]<-matrix(rep(1,nrow(mg[[l]])),ncol=1)
        }

        A<-mg[[2]]-1/nrow(mg[[1]])*ln1[[2]]%*%t(ln1[[1]])%*%mg[[1]]
        z<-t(ln1[[2]])%*%A
        C<-t(z)%*%z

        #A1<-mg[[1]]-1/nrow(mg[[2]])*ln1[[1]]%*%t(ln1[[2]])%*%mg[[2]]
        #z1<-t(ln1[[1]])%*%A1
        #C1<-t(z1)%*%z1

        #C<-C2+C1


        v<-eigen(C)$vectors[,1:c]
        loadings<-v
        #ys2<-mg[[2]]%*%v
        ys1<-mg[[1]]%*%v

        # y <- rbind(ys1,ys2)
        y<-ys1
        datapcs<-cbind(y,group)
        colnames(datapcs)<-c("Dim1","Dim2","group")
        group<-as.factor(datapcs$group)
        percentage2 <- round((eigen(C)$values / sum(eigen(C)$values))*100, 2)
        percentage1 <- round((eigen(C)$values / sum(eigen(C)$values))*100, 2)
        percentage <- paste(colnames(datapcs), "(", paste( as.character(percentage1),
                                                           "%", ")", sep="") )
        group<-factor(datapcs$group,labels=c(1,2))
        groups=group
        datapcs<-cbind(datapcs,groups)
        if(Plot){
          scatterPlot1 <- ggplot(datapcs,aes(Dim1, Dim2, color=groups)) +
            geom_point(size = 3) + xlab(percentage[1]) + ylab(percentage[2])+
            scale_color_brewer(palette="Set2",labels = as.factor(1:length(g)))+
            theme(legend.position="bottom")



          xdensity <- ggplot(datapcs, aes(Dim1, fill=groups)) +
            geom_density(alpha=.5) +
            scale_fill_brewer(palette="Set2")+
            theme(legend.position = "none")

          ydensity <- ggplot(datapcs, aes(Dim2, fill=groups)) +
            geom_density(alpha=.5) +
            scale_fill_brewer(palette="Set2")+
            theme(legend.position = "none")

          blankPlot <- ggplot()+geom_blank(aes(1,1)) +
            cowplot::theme_nothing()

          grid<-grid.arrange(xdensity, blankPlot, scatterPlot1, ydensity,
                             ncol=2, nrow=2, widths=c(4, 3), heights=c(2.5, 4))
        }

        r<-list(variates=y,loadings=loadings,prop_expl_dim=percentage2)
        r
        class(r) = c("mdr")
        return(r)
      }


      else if(length(g)>2){

        s=1
        for(l in 1:length(g)){
          mg[[l]]<-as.matrix(data.x1[[l]])
          ln1[[l]]<-matrix(rep(1,nrow(mg[[l]])),ncol=1)
        }

        for(w in 1:length(g)){
          for(r in 1:length(g)){
            if(w>r){
              A[[s]]<-mg[[w]]-1/nrow(mg[[r]])*ln1[[w]]%*%t(ln1[[r]])%*%mg[[r]]
              s=s+1
            }
          }
        }

        W=1
        for(w in 1:length(g)){
          for(r in 1:length(g)){
            if(w>r){
              z[[W]]<-t(ln1[[w]])%*%A[[W]]
              W=W+1
            }
          }
        }
        for(q in 1:length(g)){
          C[[q]]<-t(z[[q]])%*%z[[q]]
        }

        CC<-Reduce("+",C)
        v<-eigen(CC)$vectors[,1:c]
        loadings<-v

        for(p in 1:length(g)){
          ys[[p]]<-mg[[p]]%*%v
        }

        y <- as.data.frame(do.call(rbind,lapply(ys,matrix,ncol=c,byrow=F)))

        datapcs<-cbind(y,group)
        colnames(datapcs)<-c("Dim1","Dim2","group")
        group<-as.factor(datapcs$group)
        percentage2 <- round(eigen(CC)$values / sum(eigen(CC)$values)*100, 2)
        percentage1 <- round(eigen(CC)$values / sum(eigen(CC)$values)*100, 2)
        percentage <- paste(colnames(datapcs), "(", paste( as.character(percentage1),
                                                           "%", ")", sep="") )
        if(Plot){
          scatterPlot1 <- ggplot(datapcs,aes(Dim1, Dim2, color=group)) +
            geom_point(size = 3) + xlab(percentage[1]) + ylab(percentage[2])+
            scale_color_brewer(palette="Set2",labels = factor(1:length(g)))+
            theme(legend.position="bottom")



          xdensity <- ggplot(datapcs, aes(Dim1, fill=group)) +
            geom_density(alpha=.5) +
            scale_fill_brewer(palette="Set2")+
            theme(legend.position = "none")

          ydensity <- ggplot(datapcs, aes(Dim2, fill=group)) +
            geom_density(alpha=.5) +
            scale_fill_brewer(palette="Set2")+
            theme(legend.position = "none")

          blankPlot <- ggplot()+geom_blank(aes(1,1)) +
            cowplot::theme_nothing()

          grid<-grid.arrange(xdensity, blankPlot, scatterPlot1, ydensity,
                             ncol=2, nrow=2, widths=c(4, 3), heights=c(2.5, 4))
        }

        r<-list(variates=y,loadings=loadings,prop_expl_dim=percentage2)
        r
        class(r) = c("mdr")
        return(r)
      }
    }
  }
}
