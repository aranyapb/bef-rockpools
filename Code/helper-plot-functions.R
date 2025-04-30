#' @title cond_eff_plot
#' 
#' @description Code to plot the conditional effect plots that we generate for
#' the supplementary material and reuse many times.
#'
#' @param data_raw - raw data.frame
#' @param data_pred - data.frame with the expected value and confidence intervals
#' @param x_var - name of the x-variable in both datasets
#' @param y_var - name of the y-variable in both datasets
#' @param xlab - x-axis label
#' @param ylab - y-axis label
#' @param taxa - all, active, passive
#' @param labels - the labels to add to each of the three plots

# load relevant libraries
library(ggplot2)
library(ggbeeswarm)

cond_eff_plot <- function(data_raw, data_pred, 
                          x_var, y_var, size_var = NA, 
                          xlab, ylab, 
                          taxa, labels = c("a", "b", "c")) {
  
  # get correct labels
  if(taxa == "all") {
    title <- "All taxa"
    fig <- labels[1]
    ylab <- ylab
    xlab <- " "
  } else if(taxa == "active") {
    title <- "Active dispersers"
    fig <- labels[2]
    ylab <- NA
    xlab <- xlab
  } else if(taxa == "passive") {
    title <- "Passive dispersers"
    fig <- labels[3]
    ylab <- NA
    xlab <- " "
  }
  
  px <- 
    ggplot() +
    geom_quasirandom(data = data_raw, 
                     mapping = aes_string(x = x_var, 
                                          y = y_var, 
                                          colour = "Inselberg",
                                          size = if(!is.na(size_var)){size_var}else{NULL}), 
                     width = 0.025, shape = 1, stroke = 0.5, alpha = 0.5) +
    geom_line(data = data_pred, mapping = aes_string(x = x_var, y = "predicted_mean")) +
    geom_ribbon(data = data_pred, 
                mapping = aes_string(x = x_var, ymax = "upper", ymin = "lower"), alpha = 0.1) +
    ggtitle(title) +
    scale_size_continuous(range = c(1, 2.5)) +
    ylab(if(!is.na(as.character(ylab))){ylab}else{NULL}) +
    xlab(xlab) +
    scale_colour_manual(values = wes_palette(name = "Darjeeling1", n = 10, type = "continuous")) +
    guides(size = "none", colour = guide_legend(override.aes = list(size = 2.5, alpha = 1, shape = 1, stroke = 1),
                                                nrow = 1))+
    theme_meta() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.title.x = element_text(colour = "black"))
  
  return(list(px, fig))
  
}

### END
