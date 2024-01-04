################
############## Day 2
########################

setwd("C:/Users/prile/OneDrive - Washington State University (email.wsu.edu)/PhD_Documents/Courses/BIOL_592_TimeSeries")

## Packages and libraries:
#More efficient way of checking and installing packages:
packages <- c("devtools", "learnr", "stats", "TMB", "MARSS", "marssTMB", "datasets", "magrittr", "forecast")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
#load packages:
invisible(lapply(packages, library, character.only = TRUE))

install.packages('marssTMB',
                 repos = c('https://atsa-es.r-universe.dev','https://cloud.r-project.org'))
library(marssTMB)

#####
## Lab 4 # Univariate State Space Models ###
#####
# MARSS package requires matrix values

#










