#'
#' @title select_taxa
#' 
#' @description Function to choose a group of taxa to analyse i.e. all, active
#' or passive dispersers interactively
#'

select_taxa <- function() {
  
  # set the options and question
  options <- c("all", "active", "passive")
  question <- "Which group of organisms are you analysing?"
  
  # ask the user which taxon group
  switch( menu(options, title = question), "all", "active", "passive" )
  
}