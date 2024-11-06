#'
#' @title Combine plots from the slope variation analysis
#' 
#' @description Combines the plots generated when examining whether the
#' slope varies by inselberg and the correlates of the different inselberg
#' slopes with other environmental variables
#'

# load relevant libraries
library(ggpubr)
library(ggplot2)

# get a list of files
files <- list.files("Figures/")

# get the file names for fig. 6
fig6 <- files[ grepl(pattern = "_6[a-z]", files) ]

# load the files for fig. 6
fig6 <- lapply(fig6, function(x) {
  
  y <- readRDS(paste0("Figures/", x)) 
  z <- 
    y +
    scale_y_continuous(limits = c(-8.2, 7.2))
  
  return(z)
  
  } )

# remove background from legend
fig6 <- lapply(fig6, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
fig6 <- 
  ggarrange(plotlist = fig6, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom")
plot(fig6)

ggsave("Figures/fig_6.pdf", fig6, units = "cm", width = 20,
       height = 9.5)


# get the file names for fig. 7
fig7 <- files[ grepl(pattern = "_7[a-z]", files) ]

# load the files for fig. 7
fig7 <- lapply(fig7, function(x) {
  
  y <- readRDS(paste0("Figures/", x)) 
  z <- 
    y +
    scale_y_continuous(limits = c(-10, 8))
  
  return(z)
  
} )

fig7 <- lapply(fig7, function(x) {
x + theme(legend.key = element_blank())  # remove background from legend icons 
})


# combine using ggarrange
fig7 <- 
  ggpubr::ggarrange(plotlist = fig7,
                    ncol = 3, nrow = 3,
                    labels = c("(a)", " ", " ", "(b)", " ", " ", "(c)", " ", " "),
                    heights = c(1.05, 1, 1),
                    widths = c(1, 1, 1),
                    font.label = list(size = 11, face = "plain"),
                    align = c("none"),
                    common.legend = TRUE,
                    legend = "bottom")

plot(fig7)

ggsave("Figures/fig_7.pdf", fig7, 
       units = "cm", width = 21, height = 22)

### END
