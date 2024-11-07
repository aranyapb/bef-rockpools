#'
#' @title Combine plots from the SEM analysis
#' 
#' @description Combines the plots generated using the SEM analysis.
#' 

# load relevant libraries
library(ggpubr)
library(ggplot2)

# get a list of .rds plots
files <- list.files("Figures/")

# main text: fig. 5

# get relevant file names
fig5 <- files[ grepl(pattern = "_5[a-z]", files) ]
print(fig5)

# load the files for fig. 5
fig5 <- lapply(fig5, function(x) readRDS(paste0("Figures/", x) ) )

# remove background from legend
fig5 <- lapply(fig5, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
fig5 <- 
  ggarrange(plotlist = fig5, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom") 
plot(fig5)

ggsave("Figures/fig_5.pdf", fig5, units = "cm", width = 20,
       height = 9.5)

# appendix: fig. s1.2

# get relevant file names
figs1.2 <- files[ grepl(pattern = "_s1.2[a-z]", files) ]
print(figs1.2)

# load the files for fig. s1.2
figs1.2 <- lapply(figs1.2, function(x) readRDS(paste0("Figures/", x) ) )

# remove background from legend
figs1.2 <- lapply(figs1.2, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
figs1.2 <- 
  ggarrange(plotlist = figs1.2, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom")
plot(figs1.2)

ggsave("Figures/fig_s1.2.pdf", figs1.2, units = "cm", width = 20,
       height = 9.5)

# appendix: fig. S1.3

# get relevant file names
figs1.3 <- files[ grepl(pattern = "_s1.3[a-z]", files) ]
print(figs1.3)

# load the files for fig. s1.3
figs1.3 <- lapply(figs1.3, function(x) readRDS(paste0("Figures/", x) ) )

# remove titles for plots 4 to 6
for(i in 4:6) {
  figs1.3[[i]] <- figs1.3[[i]] + ggtitle(" ")
}

# remove background from legend
figs1.3 <- lapply(figs1.3, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
figs1.3 <- 
  ggarrange(plotlist = figs1.3, ncol = 3, nrow = 2, labels = c("(a)", " ", " ",
                                                             "(b)", " ", " "),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom")
plot(figs1.3)

ggsave("Figures/fig_s1.3.pdf", figs1.3, units = "cm", width = 20,
       height = 19)

# appendix: fig. S1.4

# get relevant file names
figs1.4 <- files[ grepl(pattern = "_s1.4[a-z]", files) ]
print(figs1.4)

# load the files for fig. s6
figs1.4 <- lapply(figs1.4, function(x) readRDS(paste0("Figures/", x) ) )

# remove titles for plots 4 to 6
for(i in 4:6) {
  figs1.4[[i]] <- 
    figs1.4[[i]] + ggtitle(" ")
}

# remove background from legend
figs1.4 <- lapply(figs1.4, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
figs1.4 <- 
  ggarrange(plotlist = figs1.4, ncol = 3, nrow = 2, labels = c("(a)", " ", " ",
                                                             "(b)", " ", " "),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom")
plot(figs1.4)

ggsave("Figures/fig_s1.4.pdf", figs1.4, units = "cm", width = 20,
       height = 19)

# appendix: fig. s1.5

# get relevant file names
figs1.5 <- files[ grepl(pattern = "_s1.5[a-z]", files) ]
print(figs1.5)

# load the files for fig. s1.5
figs1.5 <- lapply(figs1.5, function(x) readRDS(paste0("Figures/", x) ) )

# remove titles for plots 4 to 6
for(i in 4:6) {
  figs1.5[[i]] <- figs1.5[[i]] + ggtitle(" ")
}

# remove background from legend
figs1.5 <- lapply(figs1.5, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
figs1.5 <- 
  ggarrange(plotlist = figs1.5, ncol = 3, nrow = 2, labels = c("(a)", " ", " ",
                                                             "(b)", " ", " "),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom")
print(figs1.5)

ggsave("Figures/fig_s1.5.pdf", figs1.5, units = "cm", width = 20,
       height = 19)

# appendix: fig. s1.6

# get relevant file names
figs1.6 <- files[ grepl(pattern = "_s1.6[a-z]", files) ]
print(figs1.6)

# load the files for fig. 8
figs1.6 <- lapply(figs1.6, function(x) readRDS(paste0("Figures/", x) ) )

# remove background from legend
figs1.6 <- lapply(figs1.6, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
figs1.6 <- 
  ggarrange(plotlist = figs1.6, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom")
plot(figs1.6)

ggsave("Figures/fig_s1.6.pdf", figs1.6, units = "cm", width = 20,
       height = 9.5)

### END
