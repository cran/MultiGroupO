0#' biplot methods
#'
#' @param variates is the size of groups
#' @param loadings is a vector of classes
#' @param prop_expl_var data set
#' @param comp component numeric
#' @param group is a vector of groups
#' @param rownamevar is a logical vector where TRUE is the
#' label of the observations, if is FALSE, is index.
#' @param rownameload is a logical vector where TRUE is the
#' label of the vectors of loadings, if is FALSE, is index.
#' @return \strong{return an grafics} .
#'
#'
#' @examples
#' library(datasets)
#'obj<-pca(datos=iris[,-5],grupos=iris[,5],Plot=FALSE,center=TRUE,scale=TRUE)
#'BIplot(variates=obj$variates,loadings=obj$loadings,
#'prop_expl_var=obj$prop_expl_var,comp=c(1,2),
#'group=factor(as.numeric(iris[,5])),rownamevar=FALSE,rownameload=FALSE)
#'
#' @export
BIplot<-function(variates,loadings,prop_expl_var,comp = c(1,2),group=NULL,rownamevar=T,rownameload=T){


  scaler <- max(variates, na.rm = TRUE)/max(loadings, na.rm = TRUE)
  loadings <- loadings*scaler
  if(rownamevar==T){
    ind.names <- rownames(variates)
  }
  else{rownames(variates)<-1:nrow(variates)
  ind.names <- rownames(variates)}



  if(rownameload==T){
    var.labels <- rownames(loadings)
  }
  else{rownames(loadings)<-1:nrow(loadings)
  var.labels <- rownames(loadings)}


  PCs <- paste0('component_', comp)
  expl_vars <- prop_expl_var[comp]*100
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
      size = 0.5,
      show.legend = FALSE
    )

  angle <- rep(0, nrow(loadings))

  angle <- atan(loadings[, comp[2]]/loadings[,comp[1]]) * 360/(2 * pi)

  gg_biplot <-
    gg_biplot + geom_text_repel(
      aes(
        x = loadings[,comp[1]],
        y = loadings[,comp[2]],
        label = var.labels,
        angle = 0,
        # hjust = ifelse(loadings[, "PC1"] > 0, 1, 0),
        #  vjust = ifelse(loadings[, "PC2"] > 0, 1, 0)
      ),
      col = 'blue',
      size = 3,
      box.padding = 0.1,
      point.padding = 0.1
    )

  gg_biplot <- gg_biplot + geom_vline(xintercept = 0, size = 0.3, col = 'grey75')
  gg_biplot <- gg_biplot +  geom_hline(yintercept = 0, size = 0.3, col = 'grey75')

  gg_biplot <- gg_biplot +
    geom_point(aes(x = variates[, comp[1]],
                   y = variates[, comp[2]]
    ),
    col=group
    #col="blue"
    ,size = 2
    )



  gg_biplot <- gg_biplot +
    geom_text_repel(mapping = aes(x = variates[, comp[1]],
                                  y = variates[, comp[2]],
                                  label = ind.names

    ),
    size = 2,
    show.legend = T)

  gg_biplot <- gg_biplot + scale_y_continuous(sec.axis = sec_axis(~.*1/scaler)) +
    scale_x_continuous(sec.axis = sec_axis(~.*1/scaler))
  gg_biplot

}


