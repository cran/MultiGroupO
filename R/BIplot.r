#' biplot methods
#'

#' @param variates is the size of groups
#' @param loadings is a vector of classes
#' @param prop_expl_var data set
#' @param comp component numeric
#' @param group is a vector of groups
#' @param rownamevar is a logical vector where TRUE is the
#' label of the observations, if is FALSE, is index.
#' @import lemon
#' @param rownameload is a logical vector where TRUE is the
#' label of the vectors of loadings, if is FALSE, is index.
#' @return \strong{return an grafics} .
#'
#'
#' @examples
#' library(datasets)
#'obj<-pca(datos=iris[,-5],grupos=iris[,5],Plot=FALSE,center=TRUE,scale=TRUE)
#'BIplot(variates=obj$variates,loadings=obj$loadings,
#'       prop_expl_var=obj$prop_expl_var,comp=c(1,2),
#'       group=factor(as.numeric(iris[,5])),rownamevar=FALSE,rownameload=TRUE)
#'
#' @export
BIplot<-function(variates,loadings,prop_expl_var,comp = c(1,2),group=NULL,rownamevar=T,rownameload=T){

  scaler <- max(variates, na.rm = TRUE)/max(loadings, na.rm = TRUE)
  loadings <- loadings*scaler
  if(rownamevar==T){
    ind.names <- rownames(variates)
  }
  else{rownames(variates)<-as.character(1:nrow(variates))
  ind.names <- rownames(variates)}



  if(rownameload==T){
    var.labels <- rownames(loadings)
  }
  else{rownames(loadings)<-1:nrow(loadings)
  var.labels <- rownames(loadings)}

  PCs <- paste0('component_', comp)
  expl_vars <- round(prop_expl_var[comp]*100,2)
  axes.titles <- sprintf("%s   (%s%%)", PCs, expl_vars)

  gg_biplot <-
    ggplot() +
    theme_classic() +
    labs(x = axes.titles[1],
         y = axes.titles[2])


  gg_biplot <-
    gg_biplot + geom_segment(
      aes(
        x = 0,
        y = 0,
        xend = loadings[,comp[1]],
        yend = loadings[,comp[2]],
      ),
      col = 'red',
      arrow = arrow(length = unit(0.2, "cm")),
      size = 1,
      show.legend = FALSE) +
    geom_text_repel(aes(
      x = loadings[,comp[1]],
      y = loadings[,comp[2]],
      label = var.labels,
      angle = 0,
    ),
    col = 'blue',
    size = 3,
    box.padding = 0.3,
    point.padding = 0.1)

  gg_biplot <- gg_biplot + geom_vline(xintercept = 0, size = 0.3, col = 'grey75')
  gg_biplot <- gg_biplot +  geom_hline(yintercept = 0, size = 0.3, col = 'grey75')

  gg_biplot <- gg_biplot +
    geom_point(aes(x = variates[, comp[1]],
                   y = variates[, comp[2]],col=group, shape=group
    ),
    size = 2
    )+
    scale_color_brewer(palette="Set2",labels =as.factor(group)) +
    stat_ellipse(geom = "polygon",
                 aes(x = variates[, comp[1]],
                     y = variates[, comp[2]], fill=group),
                 alpha = 0.25)+
    scale_fill_brewer(palette="Set2")

  gg_biplot <- gg_biplot +
    geom_text_repel(mapping = aes(x = variates[, comp[1]],
                                  y = variates[, comp[2]],
                                  label = ind.names

    ),
    size = 2,
    show.legend = T,
    box.padding = 0.1,
    point.padding = 0.1
    )

  gg_biplot <- gg_biplot + scale_y_continuous(sec.axis = sec_axis(~.*1/scaler)) +
    scale_x_continuous(sec.axis = sec_axis(~.*1/scaler))

  # custom box around legend
  gg_biplot<- gg_biplot+guides(color=FALSE)+theme(legend.background = element_rect(fill = "grey90"),
                                                  legend.key.size = unit(0.5, "cm"),
                                                  legend.key.width = unit(0.5,"cm") )

  gg_biplot
}
