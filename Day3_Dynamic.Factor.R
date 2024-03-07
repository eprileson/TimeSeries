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

## now to fit DFA models:
#as before, need to create matrix with model list with forms for 
#each vector and matrix:
#OBSERVATION MODEL:
#5 observed time series and 3 hidden states, so Z will be 5x3 matrix:
#Z matrix:
Z_vals <- list("z11",  0  ,  0  ,
               "z21","z22",  0  ,
               "z31","z32","z33",
               "z41","z42","z43",
               "z51","z52","z53")
ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 3, byrow = TRUE)

#a is the offset
aa <- "zero"

##DD and d are strength of covariate and covariate
DD <- "zero" #ie. matrix(0, mm, 1)
dd <- "zero" #matrix(0, 1, wk_last)

#variance covariance matrix for obs errors
RR <- "diagonal and unequal"

## the process model - for DFA, this equals the random walk / effects structure
#num of processes:
mm <- 3

## B is an identity matrix (1)
BB <- "identity" #diag(mm)

#uu is col vector of zeros
uu <- "zero"

# cc is covariates
CC <- "zero"
cc <- "zero"

#QQ is identity matrix for process error variance
QQ <- "identity" #diag(mm)

## now fit the model:

mod_list <- list(
  #Observation model  
  Z = ZZ,
  A = aa,
  D = DD,
  d = dd,
  R = RR,
  #process model / random walks
  B = BB,
  U = uu,
  C = cc,
  Q = QQ
)

#list with model initials, everything set to 0
init_list <- list(x0 = matrix(rep(0, mm), mm, 1))

#model control params: add extra iterations in algorithm
con_list <- list(maxit = 3000, allow.degen = T)

##fit MARSS:
dfa_1 <- MARSS(y = dat,
               model = mod_list,
               inits = init_list,
               control = con_list)

##Interpreting MARSS output:
#Z outputs are the loadings of each observed time series on the 3 hidden states (5 ts, each with 3 states)
#remember that three of the Zs are 0s

#estimates of the processes are in the dfa_1$states output:
dfa_1$states #shows random values for each time point

## ROtating trends + loadings:
#use rotation matrix H

## get the estimated ZZ as a matrix
Z_est <- coef(dfa_1, type = "matrix")$Z

## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat

##rotate both Z and x:
#rotate Z factor loadings:
Z_rot <- Z_est %*% H_inv

#rotate processes:
proc_rot <- solve(H_inv) %*% dfa_1$states

## plot labels
ylbl <- phytoplankton
w_ts <- seq(dim(dat)[2])

## set up plot area
layout(matrix(c(1,2,3,4,5,6), mm, 2), widths = c(2,1))
par(mai = c(0.5, 0.5, 0.5, 0.1), omi = c(0, 0, 0, 0))

## plot the processes
for(i in 1:mm) {
  ylm <- c(-1, 1) * max(abs(proc_rot[i,]))
  ## set up plot area
  plot(w_ts,proc_rot[i,], type = "n", bty = "L",
       ylim = ylm, xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts, proc_rot[i,], lwd = 2)
  lines(w_ts, proc_rot[i,], lwd = 2)
  ## add panel labels
  mtext(paste("State",i), side = 3, line = 0.5)
  axis(1, 12 * (0:dim(dat_1980)[2]) + 1, yr_first + 0:dim(dat_1980)[2])
}

## plot the loadings
## set minimum loading to actually plot
minZ <- 0
## set y-limits for plots
ylm <- c(-1, 1) * max(abs(Z_rot))
## cycle through the states
for(i in 1:mm) {
  plot(x = c(1:N_ts)[abs(Z_rot[,i]) > minZ],
       y = as.vector(Z_rot[abs(Z_rot[,i]) > minZ, i]),
       type = "h",
       lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylm,
       xlim = c(0.5, N_ts + 0.5), col = clr)
  for(j in 1:N_ts) {
    ## plot names above or below loadings
    if(Z_rot[j,i] > minZ) {
      text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, col = clr[j])
    }
    if(Z_rot[j,i] < -minZ) {
      text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, col = clr[j])
    }
    ## add horizontal line at 0
    abline(h = 0, lwd = 1.5, col = "gray")
  } 
  ## add labels
  mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
}

#each factor loading distribution affects the respective xt state in diff ways, hence
#the different plots of each state and dis of factor loadings

##possible cross correlation, maybe due to seasonal cycles:
## set up plot area
par(mfrow = c(1,1),mai = c(0.9,0.9,0.1,0.1))

## plot CCF's
ccf(proc_rot[1,],proc_rot[2,], lag.max = 12, main="")


## Plotting the data with the model fits:

#create a function that returns the model fitted values with +/- 95% CIs

get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
  ## empty list for results
  fits <- list()
  ## extra stuff for var() calcs
  Ey <- MARSS:::MARSShatyt(MLEobj)
  ## model params
  ZZ <- coef(MLEobj, type="matrix")$Z
  ## number of obs ts
  nn <- dim(Ey$ytT)[1]
  ## number of time steps
  TT <- dim(Ey$ytT)[2]
  ## get the inverse of the rotation matrix
  H_inv <- varimax(ZZ)$rotmat
  ## check for covars
  if(!is.null(dd)) {
    DD <- coef(MLEobj, type = "matrix")$D
    ## model expectation
    fits$ex <- ZZ %*% H_inv %*% MLEobj$states + DD %*% dd
  } else {
    ## model expectation
    fits$ex <- ZZ %*% H_inv %*% MLEobj$states
  }
  ## variance in model fits
  VtT <- MARSSkfss(MLEobj)$VtT
  VV <- NULL
  for(tt in 1:TT) {
    RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[,,tt] %*% t(ZZ)
    SS <- Ey$yxtT[,,tt] - Ey$ytT[,tt,drop = FALSE] %*% t(MLEobj$states[,tt,drop = FALSE])
    VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
  }
  SE <- sqrt(VV)
  ## upper & lower (1-alpha)% CI
  fits$up <- qnorm(1-alpha/2)*SE + fits$ex
  fits$lo <- qnorm(alpha/2)*SE + fits$ex
  return(fits)
}

#now plot the ts of the 5 plankton groups with mean of the DFA fits and CIs
## get model fits & CI's
mod_fit <- get_DFA_fits(dfa_1)

## set up plotting area
par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 0, 0, 0))

## plot the fits
for(i in 1:N_ts) {
  up <- mod_fit$up[i,]
  mn <- mod_fit$ex[i,]
  lo <- mod_fit$lo[i,]
  plot(w_ts, mn, type = "n", ylim = c(min(lo), max(up)),
       cex.lab = 1.2,
       xlab = "", ylab = ylbl[i], xaxt = "n")
  axis(1, 12 * (0:dim(dat_1980)[2]) + 1, yr_first + 0:dim(dat_1980)[2])
  points(w_ts,dat[i,], pch = 16, col = clr[i])
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}


## covariates in the DFA model
#adding in the D and dt factors, eg. temp or phosphorous

#covariates cannot have missing vals; should be inputs, not data
## get temperature and phosphorus covariates
temp <- t(plank_dat[,"Temp", drop = FALSE]) #transposes temp column into ts matrix
TP <- t(plank_dat[,"TP", drop = FALSE]) #transposes phosp column into ts matrix

#fit three models, one with each separate covariate, and then one with both
## set up the model definition
mod_list = list(m = 3, R = "diagonal and unequal")

## fit model with temperature
dfa_temp <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                  control = con_list, covariates = temp)

## fit model with phosphorus
dfa_TP <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                control = con_list, covariates = TP)

## fit model with both temperature & phosphorus
dfa_both <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                  control = con_list, covariates = rbind(temp, TP))

#the form = "dfa" allows you to just add in values for covariates and 
#MARSS builds the models automatically with 0s in the corners when you 
#add in the "dfa" form

#print AICc values from the MARSS models:
print(cbind(model = c("nocovars", "Temp", "TP", "Temp and TP"),
            AICc = round(c(dfa_1$AICc, dfa_temp$AICc, dfa_TP$AICc, dfa_both$AICc))),
      quote = FALSE)

#the temp model or both covars seem to best match data: lowest AICs

#can compare with dummy covariates sine and cosine seasonal effects:

## create dummy sine and cosine waves
cos_t <- cos(2 * pi * seq(TT) / 12)
sin_t <- sin(2 * pi * seq(TT) / 12)

## combine sine & cosine into one matrix
dd <- rbind(cos_t, sin_t)

## fit model with dummy vars
dfa_seas <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE,
                  control = con_list, covariates = dd)
#actually, these sine and cosine functions as covariates fit better!

#get model fits and CIs and plot:
mod_fit <- get_DFA_fits(dfa_seas, dd = dd)

## plot the fits
par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 0, 0, 0))
for(i in 1:N_ts) {
  up <- mod_fit$up[i,]
  mn <- mod_fit$ex[i,]
  lo <- mod_fit$lo[i,]
  plot(w_ts, mn, type = "n",
       xlab = "", ylab = ylbl[i],
       xaxt = "n", cex.lab = 1.2,
       ylim = c(min(lo), max(up)))
  axis(1, 12 * (0:dim(dat_1980)[2]) + 1, yr_first + 0:dim(dat_1980)[2])
  points(w_ts, dat[i,], pch = 16, col = clr[i])
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}



















