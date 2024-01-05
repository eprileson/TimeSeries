################
############## Day 2
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
## Lab 4 # Univariate State Space Models ###
#####
# MARSS package requires matrix values

#
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

#Simulate some data for the AR(1) model - use this to represent the density dependent pop growht model or Gompertz stochastic model where u = 0
set.seed(592)
#num of time steps
TT <- 40
#strength of dens dependence (i.e, the b term)
bb <- 0.5
#ts of process errors with sd = 1
ww <- rnorm(TT, 0, 1) #white noise with data, mean = 0, sd = 1

#set initial state, set x0 to w0
xx <- ww

#for loop over time steps:
for (t in 2:TT){
  xx[t] <- bb * xx[t-1] + ww[t]
}

#now need to add observer error to the true state (data = truth + error)
vv <- rnorm(TT, 0, 0.5) # v ~ N(0, 0.5)

#Obs. data:
yy <- xx + vv #this is the simulated data but now with vv errors

#using the MARSS model: needs any data to be in matrix form; need to define parameters u (setting = 0), b, Z, a (also = 0)
#each param needs to be in matrix form using a list to define them:

mod_list <- list(
  #first the state model:
  B = matrix("b"), #AR(1)
  U = matrix(0), #defining the u term as 0 matrix
  Q = matrix("q"), #equals variance of process errors
  #observation model:
  Z = matrix(1), #Maps states (from xt in state model) to this model
  A = matrix(0), #similar to u term, we set it to 0
  R = matrix("r") #variance of obs errors
)

#define the data as a N (row, the observed time series data) x T (cols, the time steps) matrix:
YY <- matrix(yy, nrow = 1, ncol = TT)

#fit model with MARSS(); it uses the y because we can only create an observer model with the hopes that it matches onto the theoretical state model
mod_fit <- MARSS(y = YY, model = mod_list)
#output:
#MARSS fit is Estimation method: kem 
#Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
#Estimation converged in 39 iterations. 
#Log-likelihood: -70.74115 
#AIC: 149.4823   AICc: 150.6252   

#Estimate of the parameters
#R.r      1.074 #Process error variance est
#B.b      0.528 #the AR(1) coefficient
#Q.q      0.798 #The observation error variance est
#x0.x0   -0.795 #Initial state estimate
#Initial states (x0) defined at t=0

##Extracting info from marss
str(mod_fit)
mod_fits <- t(mod_fit$states) #transpose the lists of the estimates and ses to make into single col vector
mod_fits_SE<- t(mod_fit$states.se)

## model parameters:
#use coef() to extract the parameters from model:
params_est <- coef(mod_fit, type = "vector") #gets B, Q, R, and x0
#get just the b factor:
params_est["B.b"] %>%
  round(2)
mod_fit$coef

##Approximate 95% confidence intervals using the t distribution (fat tailed Normal dist), using qt() function
#upper CI:
mod_fits_CI.hi <- mod_fits + qt(p = 0.975, df = TT -1)* mod_fits_SE #take upper CI (i.e. upper 2.5% with n-1 df and multiply by SE)
#lower CI:
mod_fits_CI.low <- mod_fits - qt(p = 0.975, df = TT - 1)* mod_fits_SE

## Plotting States + CIs
#this plots the true states (using the simulated data) and model fit 
par(mfrow = c(1,1))
par(mai = c(1.2, 1, 0.3, 0), omi = c(0, 0, 0.5, 1))
plot.ts(xx, lwd = 2, type = "o", pch = 16, col = "#488fdf",
        las = 1, ylim = c(min(xx,yy), max(xx,yy)),
        ylab = expression(italic(x[t])~~or~~italic(y[t])))
# lines(yy, lwd = 2, type = "o", pch = 16, cex = 1.5, col = "#844870")
lines(mod_fits, lwd = 2, type = "l", cex = 1.5, col = "black")
lines(mod_fits_CI.hi, lwd = 2, type = "l", cex = 1.5, col = "gray")
lines(mod_fits_CI.low, lwd = 2, type = "l", cex = 1.5, col = "gray")

## model fit with observations only, but not true values(so uses yy = info + error from observation model)
par(mai = c(1.2, 1, 0.3, 0), omi = c(0, 0, 0.5, 1))
plot.ts(yy, lwd = 2, type = "o", pch = 16, cex = 1.5, col = "#844870", las = 1,
        ylim = c(min(xx, yy), max(xx, yy)), ylab = expression(italic(x[t])~~or~~italic(y[t])))
lines(mod_fits, lwd = 2, type = "l", cex = 1.5, col = "black")
lines(mod_fits_CI.hi, lwd = 2, type = "l", cex = 1.5, col = "grey")
lines(mod_fits_CI.low, lwd = 2, type = "l", cex = 1.5, col = "grey")

## model with both observations and states;
par(mai = c(1.2, 1, 0.3, 0), omi = c(0, 0, 0.5, 1))
plot.ts(xx, lwd = 2, type = "o", pch = 16, col = "#488fdf",
        las = 1, ylim = c(min(xx,yy), max(xx,yy)),
        ylab = expression(italic(x[t])~~or~~italic(y[t])))
lines(yy, lwd = 2, type = "o", pch = 16, cex = 1.5, col = "#844870")
lines(mod_fits, lwd = 2, type = "l", cex = 1.5, col = "black")
lines(mod_fits_CI_hi, lwd = 2, type = "l", cex = 1.5, col = "gray")
lines(mod_fits_CI_low, lwd = 2, type = "l", cex = 1.5, col = "gray")

## Model summaries:
#tidy(mod_fit) to get estimated params with CIs
#fitted(mod_fit) to get model estimates of mean y
#plot(mod_fit) and ggplot2::autoplot(mod_fit) to get series of plots of states, data, diagnostics

plot(mod_fit) 
ggplot2::autoplot(mod_fit)  #diagnostic plots: residuals variance, qq plot, and acf plot



###################
#################################
## Lab 5: Multivariate State Space models #########
###############
### Animal movement
library(MASS)
#use mvnorm() function from {MASS} pkg to simulate multivariate normally dist data
TT <- 40 #time steps = 40
##var-cov matrix for the process error in the state model:
QQ <- diag(c(0.009, 0.009)) # no covariance

#transpose process errors to match equational form?
ww <- mvrnorm(n = TT, mu = c(0,0), Sigma = QQ) %>% t() #mu = mean, specified to be 0 in all models where we subtract u

#initialize state vector:
xx <- ww

#set starting geo points, lat lons:
xx[,1] <- c(45, 120)

#random walk simulation: 1st row is lat, 2nd is long: this represents true data, not observed
for (t in 2:TT){
  xx[,t] <- xx[, t-1] + ww[, t]
}

#plot the true location from random walk
par(mai = c(1.2, 1, 0.3, 0),
    omi = c(0, 0, 0.5, 1))
plot(xx[2,], xx[1, ],
     pch = 16, type = "o", col = "blue",
     xlab = "Longitude (W)", ylab = "Latitude (N)")
## add start & end points
points(xx[2, 1], xx[1, 1],
       pch = 5, cex = 2, col = "blue")
points(xx[2, TT], xx[1, TT],
       pch = 0, cex = 2, col = "blue")


## now add in observational error for the movement data:
#create var-cov matrix for observed error (R values); no covariance
RR <- diag(c(0.004, 0.004))

#obs errors transposed to match the state equation:
vv <- mvrnorm(n = TT, mu = c(0,0), Sigma = RR) %>% t()

#add errors to true locations:
yy <- xx + vv

## plot the true locations
par(mai = c(1.2, 1, 0.3, 0),
    omi = c(0, 0, 0.5, 1))
plot(xx[2,], xx[1, ],
     xlim = range(xx[2, ], yy[2, ]), ylim = range(xx[1, ], yy[1, ]),
     pch = 16, type = "o", col = "blue",
     xlab = "Longitude (W)", ylab = "Latitude (N)")
## add start & end points
points(xx[2, 1], xx[1, 1],
       pch = 5, cex = 2, col = "blue")
points(xx[2, TT], xx[1, TT],
       pch = 0, cex = 2, col = "blue")
## add the obs
lines(yy[2,], yy[1,],
      type = "o", pch = 16, col = "darkgray")
points(yy[2, 1], yy[1, 1],
       pch = 5, cex = 2, col = "darkgray")
points(yy[2, TT], yy[1, TT],
       pch = 0, cex = 2, col = "darkgray")


## now fitting the model, given the observed and state data above:
#need to explicitly define the forms for all vectors and matrices in the full MARSS model
#state: matrices(xt = Bxt-1 + u + Cct-k + wt)
#obs: matrices(yt = Zxt + a + Ddt-h + vt)

#for wandering animal, need to set vectors + matrices to either 0 or 1
#create lists to generate matrices that could have both char and numeric classes:
#2x2 matrix with empty lists:
M <- matrix(list(), 2, 2)
#chars on diag, numbers on off-diag:
diag(M) <- c("a", "b")
M[1, 2] <- 1
M[2, 1] <- 0
M

#now we know how to create that type of mtrx, make the model for MARSS()
mod_moveL <- list(
  #state model first:
  B = diag(2), #2x2 identity matrix, i.e. (1)
  U = matrix(0, 2, 1), #0s since we are removing u in a 2x1 matrix
  C = matrix(0, 2, 1), #no covariate effects, so also just 0s
  c = matrix(0), #1x1 covariate 
  Q = matrix(c(list("q"),list(0),list(0),list("q")),2,2),
#Obs model
  Z = diag(2), #2x2 identity matrix, i.e. (1)
  A = matrix(0, 2, 1), #same as U essentially , offsets = 0
  D = matrix(0, 2, 1), #Mtrx of covariate effects = 0
  d = matrix(0),
  R = matrix(c(list("r"), list(0), list(0), list("r")), 2, 2)
)

#Now use MARSS to fit the RW model on the simulated movement data:
mod_rw <- MARSS(y = yy, model = mod_moveL) #yy = the observed data

#compare estimated values from MARSS model with the simulated values:
#extract estimated states and SE:
xx_hat <- mod_rw$states

#plot estimated states with the data from before:
par(mai = c(1.2, 1, 0.3, 0),
    omi = c(0, 0, 0.5, 1))
plot(xx[2,], xx[1, ],
     xlim = range(xx[2, ], yy[2, ]), ylim = range(xx[1, ], yy[1, ]),
     pch = 16, type = "o", col = "blue",
     xlab = "Longitude (W)", ylab = "Latitude (N)")
## add start & end points
points(xx[2, 1], xx[1, 1],
       pch = 5, cex = 2, col = "blue")
points(xx[2, TT], xx[1, TT],
       pch = 0, cex = 2, col = "blue")
## add the obs
lines(yy[2,], yy[1,],
      type = "o", pch = 16, col = "darkgray")
points(yy[2, 1], yy[1, 1],
       pch = 5, cex = 2, col = "darkgray")
points(yy[2, TT], yy[1, TT],
       pch = 0, cex = 2, col = "darkgray")
#Now add estimated states:
lines(xx_hat[2, ], xx_hat[1,],
      type = "o", pch = 16, col = "darkred")
#add the same start _ end points
points(xx_hat[2,1], xx_hat[1,1],
       pch = 5, cex = 2, col = "darkred") #start point
points(xx_hat[2, TT], xx_hat[1, TT],
       pch = 0, cex = 2, col = "darkred")

#the estimated track is closer than observed, but still not spot on suggesting that
#the estimates of process errors and observational errors are contributing to these diffs:

##Migration: want to estimate rates of movemnet (eg. km/day) in adition to the actual movement
#biased random walk, includes bias in movement direction

#state model: matrices(xt = xt-1 + u + wt) process errors dist in MVN(0, Q)
#Observation model: matrices(yt = xt + vt) process errors dist in MVN(0, R)

# first create the animal track from 2d random walk
TT <- 40
QQ <- diag(c(0.009, 0.009))
#bias term for migration; let's do north west
uu <- matrix(c(0.1, 0.02), 2, 1) #2x1 matrix with bias term of north +0.1 and west 0.02

#process errors transposed for w term:
ww <- mvrnorm(n = TT, mu = c(0,0), Sigma = QQ) %>% t() #transpose 100x2 matrix into 2x100

#initialize state vector:
xx <- ww
#set starting position:
xx[,1] <- c(45, 120)

#biased random walk sim:
for (t in 2:TT){
  xx[,t] <- xx[,t-1] + uu + ww[,t]
}

## plot the true locations
par(mai = c(1.2, 1, 0.3, 0),
    omi = c(0, 0, 0.5, 1))
plot(xx[2,], xx[1, ],
     pch = 16, type = "o", col = "blue",
     xlab = "Longitude (W)", ylab = "Latitude (N)")
## add start & end points
points(xx[2, 1], xx[1, 1],
       pch = 5, cex = 2, col = "blue")
points(xx[2, TT], xx[1, TT],
       pch = 0, cex = 2, col = "blue")

#Now add the observation error to true locations:
RR <- diag(c(0.004, 0.004))
#Obs errors transposed:
vv <- mvrnorm(n = TT, mu = c(0,0), Sigma = RR) %>% t()

# add errors to true location to get obs locations:
yy <- xx + vv

#plot observed on top of true:
par(mai = c(1.2, 1, 0.3, 0),
    omi = c(0, 0, 0.5, 1))
plot(xx[2,], xx[1, ],
     pch = 16, type = "o", col = "blue",
     xlab = "Longitude (W)", ylab = "Latitude (N)")
## add start & end points
points(xx[2, 1], xx[1, 1],
       pch = 5, cex = 2, col = "blue")
points(xx[2, TT], xx[1, TT],
       pch = 0, cex = 2, col = "blue")
#add the observed:
lines(yy[2,], yy[1,],
      type = "o", pch = 16, col = "darkgrey")
points(yy[2,1], yy[1,1],
       pch = 5, cex = 2, col = "darkgrey")
points(yy[2,TT], yy[1, TT],
       pch = 0, cex = 2, col = "darkgrey")


#now fit the model with MARSS:
#need to include all estimates as before

mod_migr <- list(
  #state model first:
  B = diag(2), #2x2 identity matrix, i.e. (1)
  U = matrix(c("NS", "EW"),2, 1), #unconstrained matrix for the bias effect (NW migration)
  C = matrix(0, 2, 1), #no covariate effects, so also just 0s
  c = matrix(0), #1x1 covariate 
  Q = matrix(c(list("q"),list(0),list(0),list("q")),2,2),
  #Obs model
  Z = diag(2), #2x2 identity matrix, i.e. (1)
  A = matrix(0, 2, 1), #same as U essentially , offsets = 0
  D = matrix(0, 2, 1), #Mtrx of covariate effects = 0
  d = matrix(0),
  R = matrix(c(list("r"), list(0), list(0), list("r")), 2, 2)
)

#now fit with the MARSS() function:
mod_mfit <- MARSS(y = yy, model = mod_migr) #Note the bias estimates, larger for N, smaller for W

#extract states from estimates and add to the plot:
xx_hat <- mod_mfit$states

#plot with true, observed
par(mai = c(1.2, 1, 0.3, 0),
    omi = c(0, 0, 0.5, 1))
plot(xx[2,], xx[1, ],
     pch = 16, type = "o", col = "blue",
     xlab = "Longitude (W)", ylab = "Latitude (N)")
## add start & end points
points(xx[2, 1], xx[1, 1],
       pch = 5, cex = 2, col = "blue")
points(xx[2, TT], xx[1, TT],
       pch = 0, cex = 2, col = "blue")
#add the observed:
lines(yy[2,], yy[1,],
      type = "o", pch = 16, col = "darkgrey")
points(yy[2,1], yy[1,1],
       pch = 5, cex = 2, col = "darkgrey")
points(yy[2,TT], yy[1, TT],
       pch = 0, cex = 2, col = "darkgrey")
#add the estimates:
lines(xx_hat[2,], xx_hat[1,],
      type = "o", pch = 16, col = "darkred")
points(xx_hat[2, 1], xx_hat[1,1],
       pch = 5, cex = 2, col = "darkred")
points(xx_hat[2,TT], xx_hat[1,TT],
       pch = 0, cex = 2, col = "darkred")
#pretty cool; much stronger estimates, likely because of the small bias estimator?


##home range, now use AR(1)model or mean-reverting model
#doesn't allow for unrestricted movement or directional movement (ie stationary model)
#add 2 terms: the B matrix (the coefficient that allows xt-1 to vary in lat and lon direction)
#so this B term is a matrix with lat and lon on diagonal; but don't have to be diff for other models

#create animal track again:
TT <- 40

#generate process error:
QQ <- diag(c(0.009, 0.009))

ww <- mvrnorm(n = TT, mu = c(0,0), Sigma = QQ) %>% t()

BB <- matrix(c(0.7, 0, 0, 0.4), 2, 2)

xx <- ww

#set initial point:
xx[,1] <- c(0, 0)

## simulate the true state data: need to multiply the B matrix by xt-1 term
for (t in 2:TT){
  xx[,t] <- BB %*% xx[,t-1] + ww[,t]
}

#plot true locations:
par(mai = c(1.2, 1, 0.3, 0),
    omi = c(0, 0, 0.5, 1))
plot(xx[2,], xx[1, ],
     pch = 16, type = "o", col = "blue",
     xlab = "Longitude (W)", ylab = "Latitude (N)")
## add start & end points
points(xx[2, 1], xx[1, 1],
       pch = 5, cex = 2, col = "blue")
points(xx[2, TT], xx[1, TT],
       pch = 0, cex = 2, col = "blue")


## now add observation data;
RR <- diag(c(0.004, 0.004))

vv <- mvrnorm(n = TT, mu = c(0,0), Sigma = RR) %>% t()

yy <- xx + vv

#plot:
lines(yy[2,], yy[1,],
      type = "o", pch = 16, col = "darkgrey")
points(yy[2,1], yy[1,1],
       pch = 5, cex = 2, col = "darkgrey")
points(yy[2, TT], yy[1, TT],
       pch = 0, cex = 2, col = "darkgrey")


#now fit MARSS model:
mod_Home <- list(
  #state:
  B = "diagonal and equal", # 2x2 matrix with mean reversion on diag i.e, the coefficient in matrix form times the xt-1 matrix term
  U = "zero",
  C = "zero",                  # 2x1 vector of covariate effects = 0
  c = "zero",                  # 1x1 matrix of covariates (none)
  Q = "diagonal and equal",    # 2x2 covariance matrix
  ## obs model
  Z = "identity",              # 2x2 identity matrix
  A = "zero",                  # 2x1 vector of offsets = 0
  D = "zero",                  # 2x1 vector of covariate effects = 0
  d = "zero",                  # 1x1 matrix of covariates (none)
  R = "diagonal and equal"     # 2x2 covariance matrix
)

#now fit the model:
#first define control list, so that you define maximum iterations of MARSS
con_list <- list(maxit = 2000)
home_fit <- MARSS(y = yy, model = mod_Home)

#extract estimates:
xx_hat <- home_fit$states

## add the estimated states to plot
lines(xx_hat[2,], xx_hat[1,],
      type = "o", pch = 16, col = "darkred")
points(xx_hat[2, 1], xx_hat[1, 1],
       pch = 5, cex = 2, col = "darkred")
points(xx_hat[2, TT], xx_hat[1, TT],
       pch = 0, cex = 2, col = "darkred")

#follows the observed track much more closely; 
#the MARSS estimate decided to go for an R value much closer to 0 than for Q, so it is close to on top of the observed track



#############################
##### Lab 6: Estimating species interactions ###
##########################################3
## set the years
yr_set <- seq(1960, 2011)

## select the wolf & moose data & log-transform them
royale_data <- isleRoyal[isleRoyal[,"Year"] %in% yr_set,
                         c("Wolf", "Moose")] %>%
  log()

## set of years
tt <- isleRoyal[, "Year"]

## plot the data
par(mai = c(1.2, 1, 0.3, 0),
    omi = c(0, 0, 0.5, 1))
matplot(yr_set, royale_data,
        type = "o", pch = 16, lwd = 2, lty = "solid",
        col = c("darkgray", "brown"),
        ylab = "Log of counts", xlab = "Year"
)
## add a legend
legend(x = "center", legend = c("Moose", "Wolves"),
       col = c("brown", "darkgray"),lwd = 2, bty = "n")


## fit the model using list() again to use in MARSS function:

wolf_moose.list <- list(
  #state:
  B = "unconstrained", #diff items in each element
  U = "unequal",
  C = "zero",
  c = "zero",
  Q = "diagonal and unequal",
#observation model
  Z = "identity", #identity matrix
  A = "zero",
  D = "zero",
  d = "zero",
  R = "diagonal and unequal"
)  

royale_data <- t(royale_data) #transpose data to use for MARSS

mod_first <- MARSS(royale_data, model = wolf_moose.list)
#model did not converge after 500 iterations because both U and B are being estimated

##remove the mean to standardize variance, using z score:
royale_data_z <- zscore(royale_data)

#update the model list for MARSS by setting u = 0, then refit
wolf_moose.list$U <- "zero"

mod_second <- MARSS(royale_data_z, model = wolf_moose.list)
#still not converged, reinspect log counts of moose and wolves; look for patterns like RW or white noise:
#looks more like white noise

#set R (the observation error term) to 0 assuming that wolves and moose are easy to count fully
wolf_moose.list$R <- "zero"

#update model
mod_third <- MARSS(royale_data_z, model = wolf_moose.list)

#this works and you can gather various effects of moose - wolf and wolf - moose interactions from the B terms

##Now examine interaction terms / effects:
#use coef() to extract elements of B:

B_hat <- coef(mod_third, type = "matrix")$B

#reconfigure row names:
rownames(B_hat) <- colnames(B_hat) <- rownames(royale_data)
#inspect values
print(B_hat, digits = 2)
#these values are not actually the density dependence factors as this is really 1 - bi,j. We therefore need to subtract the B values from 1, or the identity matrix minus the B values

#get strength of density dependence by subtracting diagonal elements of B from 1
#calculate
round(diag(2) - B_hat, 2) #identity matrix - B matrix; round to two decimals

## adding covariates to model:
#add in temperature and precip 3 yr averages to the model as C terms:
#only looking at env effects on moose, not wolves, so the C matrix has 0s in the top row representing no effect on wolves

#firsr prepare covariate data:

## names of climatic covariates
covar_names <- c("jan.feb.ave.temp", "jan.feb.ave.precip", "july.sept.ave.temp")

## set the years for the covariates
yr_covars <- seq(1959, 2010) #includes previous year from start of dataset and one below end point

## select the appropriate covariates
covars <- isleRoyal[isleRoyal[,"Year"] %in% yr_covars, covar_names]

## transpose the covariates & z-score for MARSS()
covars_z <- t(covars) %>% zscore()

## rename covariates to match model
rownames(covars_z) <- c("WT", "WP", "ST")

##multicollinearity:
#since temp and precip are highly correlated, some variance may be incorrectly attributed to some of covariates and give incorrect estimates
#create function to estimate correlation:
cor_fun <- function(x, y){
  text(0.5, 0.5, format(cor(x, y), digits = 2), cex = 2)
}

pairs(t(covars_z), lower.panel = cor_fun)
#not much correlation between winter temp and precip

#update the MARSS List with the covariate term added for C:
wolf_moose.list <- list(
  #state:
  B = "unconstrained", #diff items in each element
  U = "unequal",
  C = matrix(list(  0,    0,    0,
                "WT", "WP", "ST"),
         nrow = 2, ncol = 3, byrow = TRUE),
  c = covars_z,
  Q = "diagonal and unequal",
  #observation model
  Z = "identity", #identity matrix
  A = "zero",
  D = "zero",
  d = "zero",
  R = "zero"
)  

#Now refit model:
mod_fourth <- MARSS(royale_data_z, model = wolf_moose.list)

#get the new interaction estimates:
B_hat <- coef(mod_fourth, type = "matrix")$B

#reconfigure row names:
rownames(B_hat) <- colnames(B_hat) <- rownames(royale_data)
#inspect values
print(B_hat, digits = 2)
#these values are not actually the density dependence factors as this is really 1 - bi,j. We therefore need to subtract the B values from 1, or the identity matrix minus the B values

#get strength of density dependence by subtracting diagonal elements of B from 1
#calculate
round(diag(2) - B_hat, 2) #identity matrix - B matrix; round to two decimals



