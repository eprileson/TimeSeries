###################################################################################
######################## Time Course Series - WSU Pullman ####################
###############################################################################

setwd("C:/Users/prile/OneDrive - Washington State University (email.wsu.edu)/PhD_Documents/Courses/BIOL_592_TimeSeries")

## Packages and libraries:
#More efficient way of checking and installing packages:
packages <- c("devtools", "learnr", "stats", "MARSS", "datasets")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
#load packages:
invisible(lapply(packages, library, character.only = TRUE))


#####
##### Actual lectures at Pullman course

#make an annual time series dataset:
xt <- round(rnorm(20), 2)  #annual values

#specialized ts data function to work with ts as special data type
ts(xt, start = 2001, end = 2020, frequency = 1)

#specialized plotting function; 
plot.ts(xt, col = "blue", las = 1)

#monthly time series data:
mt <- round(rnorm(36), 2)

ts(mt, start = c(2001, 1), end = c(2003, 12), frequency = 12) #or last arg could be detlat = 1/12 so change in time = 1 / 12 months











