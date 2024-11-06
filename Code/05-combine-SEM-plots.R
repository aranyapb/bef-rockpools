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

# appendix: fig. s4 (S1.2)

# get relevant file names
figs4 <- files[ grepl(pattern = "_s4[a-z]", files) ]
print(figs4)

# load the files for fig. s4
figs4 <- lapply(figs4, function(x) readRDS(paste0("Figures/", x) ) )

# remove background from legend
figs4 <- lapply(figs4, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
figs4 <- 
  ggarrange(plotlist = figs4, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom")
plot(figs4)

ggsave("Figures/fig_s4.pdf", figs4, units = "cm", width = 20,
       height = 9.5)

# appendix: fig. s5 (S1.3)

# get relevant file names
figs5 <- files[ grepl(pattern = "_s5[a-z]", files) ]
print(figs5)

# load the files for fig. s5
figs5 <- lapply(figs5, function(x) readRDS(paste0("Figures/", x) ) )

# remove titles for plots 4 to 6
for(i in 4:6) {
  figs5[[i]] <- figs5[[i]] + ggtitle(" ")
}

# remove background from legend
figs5 <- lapply(figs5, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
figs5 <- 
  ggarrange(plotlist = figs5, ncol = 3, nrow = 2, labels = c("(a)", " ", " ",
                                                             "(b)", " ", " "),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom")
plot(figs5)

ggsave("Figures/fig_s5.pdf", figs5, units = "cm", width = 20,
       height = 19)

# appendix: fig. s6 (S1.4)

# get relevant file names
figs6 <- files[ grepl(pattern = "_s6[a-z]", files) ]
print(figs6)

# load the files for fig. s6
figs6 <- lapply(figs6, function(x) readRDS(paste0("Figures/", x) ) )

# remove titles for plots 4 to 6
for(i in 4:6) {
  figs6[[i]] <- 
    figs6[[i]] + ggtitle(" ")
}

# remove background from legend
figs6 <- lapply(figs6, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
figs6 <- 
  ggarrange(plotlist = figs6, ncol = 3, nrow = 2, labels = c("(a)", " ", " ",
                                                             "(b)", " ", " "),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom")
plot(figs6)

ggsave("Figures/fig_s6.pdf", figs6, units = "cm", width = 20,
       height = 19)

# appendix: fig. s7 (S1.5)

# get relevant file names
figs7 <- files[ grepl(pattern = "_s7[a-z]", files) ]
print(figs7)

# load the files for fig. s7
figs7 <- lapply(figs7, function(x) readRDS(paste0("Figures/", x) ) )

# remove titles for plots 4 to 6
for(i in 4:6) {
  figs7[[i]] <- figs7[[i]] + ggtitle(" ")
}

# remove background from legend
figs7 <- lapply(figs7, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
figs7 <- 
  ggarrange(plotlist = figs7, ncol = 3, nrow = 2, labels = c("(a)", " ", " ",
                                                             "(b)", " ", " "),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom")
print(figs7)

ggsave("Figures/fig_s7.png", figs7, units = "cm", width = 20,
       height = 19)

# appendix: fig. s8 (S1.6)

# get relevant file names
figs8 <- files[ grepl(pattern = "_s8[a-z]", files) ]
print(figs8)

# load the files for fig. 8
figs8 <- lapply(figs8, function(x) readRDS(paste0("Figures/", x) ) )

# remove background from legend
figs8 <- lapply(figs8, function(x) {
  x + theme(legend.key = element_blank())
})

# combine using ggarrange
figs8 <- 
  ggarrange(plotlist = figs8, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 11, face = "plain"),
            widths = c(1.1, 1, 1),
            common.legend = TRUE, legend = "bottom")
plot(figs8)

ggsave("Figures/fig_s8.pdf", figs8, units = "cm", width = 20,
       height = 9.5)

### END
