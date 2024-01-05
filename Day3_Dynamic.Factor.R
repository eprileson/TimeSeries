################
############## Day 3
########################

setwd("C:/Users/prile/OneDrive - Washington State University (email.wsu.edu)/PhD_Documents/Courses/BIOL_592_TimeSeries")

## Packages and libraries:
#More efficient way of checking and installing packages:
packages <- c("devtools", "learnr", "stats", "TMB", "MARSS", "marssTMB", "datasets", "magrittr", "forecast", "MASS")
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
## Lab 6 Dynamic factor analysis
##########################
#DFA: try to explain temporal variation in set of n observed time series
# using linear combos of a set of m hidden random walks

#DFA constraints:
# a -> first m values are set to 0 or the whole vector is set to 0
# top right corner of Z matrix set to 0, unless Z is a nx1 matrix
# Q (process error variance) is set to the identity matrix Im ( or diag(1))

#different error structures for the observation error variance (R)
# could be unconstrained with all unique values, or could change assumptions to have all independently different error variances
# or with same variance but not independent, so they have same covariance

## Lake Washington Phytoplankton:
#Load data
data(lakeWAplankton, package = "MARSS")

#want transformed data so 0s are NAs and data is scaled / z-scored:
all_dat <- lakeWAplanktonTrans

#use first 10 years
yr_first <- 1980
yr_last <- 1989
plank_dat <- all_dat[all_dat[, "Year"] >= yr_first &
                       all_dat[, "Year"] <=yr_last,]

#create vector of phytoplankton groups:
phytoplankton <- c("Cryptomonas", "Diatoms", "Greens", "Unicells", "Other.algae")

#get only phytoplankton:
dat_1980 <- plank_dat[, phytoplankton]

#transpose, so time is columns and rows are plankton obs:
dat_1980 <- t(dat_1980)

#get number of time series:
N_ts <- dim(dat_1980)[1] #5 time series

#length of ts:
TT <- dim(dat_1980)[2] #120

##demean data (so that the a value = 0)
y_bar <- apply(dat_1980, 1, mean, na.rm = T) #gives mean of each phytoplanton

#subtract the mean:
dat <- dat_1980 - y_bar

#assign new col names from old to new dataset:
spp <- rownames(dat_1980)
rownames(dat) <- spp
head(dat)

##create time series plots of all 5 phyto groups:
clr <- c("brown", "blue", "darkgreen", "darkred", "purple")

## initialize a counter
cnt <- 1 #for the for loop to know where it is and move ahead by 1

## set up plotting space & make plots
par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 0, 0, 0))

#Plot:
for (i in spp){
  plot(dat[i,], bty = "L", xaxt = "n", pch = 16,
       xlab = "",
       ylab = "Abundance index", col = clr[cnt], type = "b")
  axis(1, 12 * (0:dim(dat_1980)[2]) + 1, yr_first + 0:dim(dat_1980)[2])
  title(i)
  cnt <- cnt + 1
}



















