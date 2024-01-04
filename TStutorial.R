################################################################3
#####################3 Time Series Tutorial ############################################]
########################################################

setwd("C:/Users/prile/OneDrive - Washington State University (email.wsu.edu)/PhD_Documents/Courses/BIOL_592_TimeSeries")

## Packages and libraries:
#More efficient way of checking and installing packages:
packages <- c("devtools", "learnr", "stats", "MARSS", "datasets", "forecast")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
#load packages:
invisible(lapply(packages, library, character.only = TRUE))

#tutorial to go with Applied time series analysis practice;
# https://atsa-es.github.io/atsa-labs/chap-mlr.html

#matrix tutorial shiny App
devtools::install_github("nwfsc-timeseries/atsalibrary")
learnr::run_tutorial("matrix", package = "atsalibrary")


####Practice with Matrix math:

#matrix multiplication:
#number of cols in matrix on left must = number of rows in matrix on right
# use the %*% symbol to multiply matrices
# 
A = matrix(1:6, 2, 3) #2 rows, 3 cols
B = matrix(1:6, 3, 2) #3 rows, 2 cols
A %*% B #gives 2 x 2 matrix since it starts with rows from left
#or
B %*% A #going the other way, it starts with 3 rows, so gives a 3x3 matrix
#B%*%A does not work because num rows != num cols

#matrix addition: matrices must have same dimensions
A+A
#B+A does not work, but using the t() function (transpose) to switch num rows + cols, so a 2x3 matrix becomes 3x2

A + t(B)

A%*%A #doesn't work, throws error;

A%*%t(A) #works since we use transpose which conforms to matrix mult rules

#subsetting or replacing elements in a matrix:
#use [] brackets
A = matrix(1, 3, 3) #3x3 matrix with 1s in each position
A[1:2, 2] = 2
A
A[c(1, 3), c(1, 3)] = 2 #replace 1s with 2s at position 1.1, 1.3
A

#inverse of matrices- for an inverse to exist, matrix has to be square
#variance - covariance matrices exist

#uses the solve() function:
A <- diag(3, 3) + matrix(1, 3, 3)
A
invA <- solve(A)
invA %*% A

#can use chol2inv() to invert matrix:
A <- diag(3, 3) + matrix(1, 3, 3)
invA <- chol2inv(chol(A))
invA %*% A


mat1 <- matrix(1, 4, 3)
mat2 <- matrix(2, 3, 4)
mat1 %*% mat2


#practice probs: matrix math:

#1. 4x3 matrix w nums 1:3 in each col. Repeat with nums 1:4 in each row

A <- matrix(1:3, 3, 4)
A

B <- matrix(1:4, 3, 4, byrow = TRUE)
B

#2. extract elements in 1st and 2nd rows and 2st and 2nd cols:

C <- A[1:2, 1:2]
C

#3: 4x3 matrix w nums 1:12 by row
byrow <- matrix(1:12, 4, 3, byrow = T)
byrow

#4: Extract 3rd row from above; show how you end up w/ vector and 
#how you ened up with 1x3 matrix:

#vector answer:
v <- byrow[3, ]
v
class(v)

m <- byrow[3,, drop = FALSE]
m
class(m)

#5. 4x3 matrix that is all 1s except for 2 in the 2,3 position:
ones <- matrix(1, 4, 3)
ones[2,3] = 2
ones

#6 take transpose of #5:
t(ones) #i.e., 2 is now in 3, 2 position (3rd row, 2nd col)

#7: 4x4 matrix with 1:4 on diag:
four <- diag(1:4)
four

#8 5x5 identity matrix:
five <- diag(1, 5)

#9 replace diag in five with 2
two <- diag(2, 5)
two

#10 2s on diag, 1s on off diag
two_off <- matrix(1, 5,5)
diag(two_off) <- 2

#11. take inverse of ten:
inv <- solve(two_off)
inv

#12. 3x3 matrix, first 9 letters of alphabet
alpha <- matrix(letters[1:9], 3, 3)
alpha

#13: replace diag with word cat:
alpha[1, 1] = "c"
alpha[2, 2] = "a"
alpha[3,3] = "t"

#14 4x3 matrix, all 1s, multiply by 3x4 matrix with all 2s
ones <- matrix(1, 4, 3)
twos <- matrix(2, 3, 4)
ones %*% twos

#15 AA is 4x3 matrix, AA Possible? = no
# AAt = 3x4 = yes (col of left = col of right)

A <- matrix(1, 4, 3)
At <- t(A)
A %*% At 

#16  ??
A <- matrix(1:9, 3, 3)
B <- matrix(c(1,0,0,0,1,1,1,1,0), 3, 3)
B
A %*% B

#17 same as above, but C = 2A
A <- matrix(1:9, 3, 3)
B <- diag(2,3)
B
A %*% B

#18 build b matrix so that AB = C and row sums for each row
A <- matrix(1:9, 3, 3)
B <- matrix(1, 3, 1)
B
A %*% B

#19 build B matrix to compute column sums of A with BA = C
A <- matrix(1:9, 3, 3)
B <- matrix(1, 3, 1)

t(B) %*% A

### 
# chapter 2: linear regression in matrix form

#other practice:
data(AirPassengers)
AP <- AirPassengers
layout(1:2) #can use instead of par(mfrow)
plot(aggregate(AP)) #aggregate to get annual trend, remove seasonal
boxplot(AP ~ cycle(AP)) #summary of each season's values:

#use window function to view parts of a ts"

win <- window(AP, start = c(1949, 1), end = c(1955, 5))
Newtime <- time(win)
layout(1:1)
plot(win); abline(reg = lm(win ~ Newtime))




#import data
data(stackloss, package = "datasets")
dat <- stackloss[1:4, ] #subset first 4 rows
dat

#simple linear model:
fit <- lm(stack.loss ~ Air.Flow, data = dat)
#equivalent to this form: stack.loss(i) = intercept + slope*air.flow + error(i) ; where e is a function of N(0, variance) 

#now matrix form:
#putting data into matrix form:
stack_m <- matrix(dat$stack.loss, 4, 1)
air.m <- matrix(dat$Air.Flow, 4, 1)

#matrix form should be: y = Zx + e
Z <- model.matrix(fit) #: makes matrix of air flow but also provides intercept of 1
Z[1:4, ]

#what is R doing behind the scenes with lm() on the data?
y <- matrix(dat$stack.loss, ncol = 1)
Z = cbind(1, dat$Air.Flow)
solve(t(Z) %*% Z) %*% t(Z) %*% y #-11.615, 0.641

#match with lm() output:
coef(lm(stack.loss ~ Air.Flow, data = dat)) #same output

## now try with multiple variables, now Z = matrix with 4 cols
fit1.mult <- lm(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., data = dat)

#Of the same form as reg linear model but uses a matrix with a ones column, then a col for each data input and type as 
#element value (row, col). Since 3 variables, matrix is 3 cols * the matrix of the slope for each var and one value for the intercept

Z <- model.matrix(fit1.mult)

#repeat to see what lm does under the hood:

y <- matrix(dat$stack.loss, ncol = 1)
Z <- cbind(1, dat$Air.Flow, dat$Water.Temp, dat$Acid.Conc.)
solve(t(Z) %*% Z) %*% t(Z) %*% y

coef(lm(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., data = dat))
#same output

##can write matrix model in transpose form: y = Dd
y <- matrix(dat$stack.loss, nrow = 1)
d <- rbind(1, dat$Air.Flow, dat$Water.Temp, dat$Acid.Conc.)
y %*% t(d) %*% solve(d %*% t(d))

#works because the transpose of Z (the variables in the matrix model) can
#still multiply through

#form2: y = Z(x) + a + e; where (x) = a one column vector of each of the variables
#and where Z is a diagonal matrix with Beta on the diagonal; each new variable adds a new 
#diagonal matrix with beta on diag with rows and cols = n*k with  n = number of instances of each var and k = # of variables

#eg, single matrix for all vars
x <- matrix(c(dat$Air.Flow, dat$Water.Temp, dat$Acid.Conc.), ,1)
Z <- matrix(0, 4, 4)
diag(Z) = 'B'
#form 2 more typical for being read by humans

#add new col for 'region' as an intercept that might be affecting data
dat <- cbind(dat, reg = rep(c('n', 's'), 4)[1:4])
dat

#now model:
fit2 <- lm(stack.loss ~ -1 + Air.Flow + reg, data = dat) #-1 gets rid of the original y intercept so we have just the regions as intercepts
coef(fit2)

#write in form 1 matrix model:
Z <- model.matrix(fit2) #fits 0s and 1s for the intercepts as separate cols; could add them together to show how R is doing this under the hood
Z[1:4,]

#solve for params:
y <- matrix(dat$stack.loss, ncol =1)
solve(t(Z) %*% Z) %*% t(Z) %*% y #Matches params from model below
coef(fit2) 

#could also write in form 2 (where Z = the diag with slope on diagonal)

## groups of betas:
#now add in cols to dataset with slopes from effects of diff owners, S, and A:
#add in operator variable as factor to dataframe
dat <- cbind(dat, owner = c('s', 'a'))
dat

#fit new model
fit3 <- (lm(stack.loss ~ -1 + Air.Flow:owner + reg, data = dat)) # : adds the owner operator to the model

#now we have a more complex model with mult intercepts and groups of Betas
#using model matrix will save time

Z <- model.matrix(fit3)
Z[1:4,] #again, R codes n and s or a and s params as 1s or 0s depending on if present or not

#solve for params to see if matches coefficients of model:
y <- matrix(dat$stack.loss, ncol = 1)
solve(t(Z) %*% Z) %*% t(Z) %*% y


## seasonal factor
#add in seasonal effect imagining data were taken consecutively by quarter:

dat <- cbind(dat, qtr = paste(rep("qtr", 4), 1:4, sep = ""))
dat

#model with intercept included
fit4 <- lm(stack.loss ~ -1 + qtr, data = dat)

Z <- model.matrix(fit4)
Z[1:4, ]

#model w/o intercepts:
fit5  <- lm(stack.loss ~ qtr, data = dat)
Z <- model.matrix(fit5)
Z[1:4, ]

#Pretty basic since there's no explanatory variables, just diff intercepts

## now include seasonal intercept + variables w/ full dataset
data(stackloss, package = "datasets")
fulldat <- stackloss
n <- nrow(fulldat)
fulldat <- cbind(fulldat, owner = rep(c("sue", "aneesh", "joe"),
                      n)[1:n], qtr = paste("qtr", rep(1:4, n)[1:n], sep = ""),
                      reg = rep(c("n","s"), n)[1:n])
#fit basic model with airflow:
# stackloss (i) = alpha(j) + Beta(j,k)Air(i) + e(i) ; here, k is the owner and j is the quarter
#there should be 4x3 betas (4 quarters and 3 diff owners)

fit7 <- lm(stack.loss ~ -1 + qtr + Air.Flow:qtr:owner, data = fulldat)

#lets look at Z for form1 of model matrix:
Z <- model.matrix(fit7)

##Models with confounding params:
fit8 <- lm(stack.loss ~ -1 + Air.Flow + reg + qtr, data = fulldat)
Z <- model.matrix(fit8)
Z #in this model, reg and qtr are left as NAs since the model can't estimate both intercepts at the same time
# in this case, one of the factors affecting intercept(owner or region) was non-identifiable


##solving for params of model form 2:
#write form2 of model (with Z = diagonal matrix with beta on diag) into vector:
#vec takes matrix and stacks columns ontop of each other in one column
#eg
A <- matrix(1:6, nrow = 2, byrow = T)
vecA = matrix(A, ncol = 1)

#transpose the x (the explanatory var) and multiply by identity matrix with 1 on diag, then multiply by vector of Z
# t(x) %*% In %*% vecZ  #then make permutation and multiply

#code:
#make y and x matrices
y <- matrix(dat$stack.loss, ncol = 1)
x <- matrix(c(1, dat$Air.Flow), ncol = 1)

#Make Z matrix
n <- nrow(dat)
k <- 1
#Z has n rows and 1 col for intercept (alpha) and n cols for datapts
#list allows combining characters and numbers
Z <- matrix(list(0), n, k*n+1)
Z[, 1] <- "alpha"
diag(Z[1:n, 1 + 1:n]) = "beta"

#creates permutation function
P <- MARSS:::convert.model.mat(Z)$free[,,1]
M <- kronecker(t(x), diag(n)) %*% P
solve(t(M) %*% M) %*% t(M) %*% y

#chapter 2 problems:
#load dataset:
data(airquality, packages = "datasets")
#remove NAs:
airquality <- na.omit(airquality)
#make month factor:
airquality$Month <- as.factor(airquality$Month)
#add region factor:
airquality$region <- rep(c("north", "south"), 60)[1:111]

#use only first 5 rows for hw dataset:
homeworkdat <- airquality[1:5, ]

#Q1. 
fit <- lm(Ozone ~ Wind + Temp, data = homeworkdat)
Z <- model.matrix(fit)
Z
#Q2. write out R code to create y and Z matrices, then solve for x (params)

#a. 
y <- matrix(homeworkdat$Ozone, ncol = 1)
Z <- cbind(1, homeworkdat$Wind, homeworkdat$Temp)

x <- solve(t(Z) %*% Z) %*% t(Z) %*% y
x
#confirm with mod:
coef(fit)

#Q3. Add -1 to model (ie. an intercept model)
#a what changes:
fit1 <- lm(Ozone ~ -1 + Wind + Temp, data = homeworkdat)
coef(fit1) #removes intercepts, so now only slopes / betas in model that are slightly diff

#b. write out form 1 as equation: y matrix = zmatrix * alpha intercept matrix*x matrix + e matrix
Z <- model.matrix(fit1) 
y <- matrix(homeworkdat$Ozone, ncol = 1)
x <- solve(t(Z) %*% Z) %*% t(Z) %*% y
#c
coef(fit1)

#Q4. for mod in Q1, write in form 2 eq
y <- matrix(homeworkdat$Ozone, ncol = 1)
x <- matrix(c(1, homeworkdat$Wind, homeworkdat$Temp), ncol = 1)

n <- nrow(homeworkdat)
diag(n)
k <- 1
#z matrix to combine character + numeric vectors
Z <- matrix(list(0), n, k*1+1)
Z[, 1] = "alpha" #first col = alpha
diag(Z) = "beta" #makes diagonals of sq matrix = beta

#create permutation matrix
P <- MARSS:::convert.model.mat(Z)$free[,,1]
M <- kronecker(t(x),diag(n))%*%P
solve(t(M) %*% M) %*% t(M) %*% y


#Q5: form 1 and solve for params
#a. ozone matrix =  region (alpha) (0 or 1 for n and s) + e
fit <- lm(Ozone ~ -1 + region, data = homeworkdat)
y <- matrix(homeworkdat$Ozone, ncol = 1)
Z <- model.matrix(fit)
Z <- cbind(c(1,0,1,0,1), c(0,1,0,1,0))

solve(t(Z) %*% Z) %*% t(Z)%*%y
coef(fit)

#Q6. same model as Q5, write form 2, write out x and Z, solve for params
y <- matrix(homeworkdat$Ozone, ncol = 1)
x <- matrix(c(1, homeworkdat$region), ncol = 1)
n <- nrow(homeworkdat)
k <- 1
Z <- matrix(list(0), n, k*1+1)
Z[,1] = "alpha"
diag(Z) = "beta"

P <- MARSS:::convert.model.mat(Z)$free[,,1]
M <- kronecker(t(x),diag(n))%*% P
solve(t(M) %*% M) %*% t(M) %*% y

#Q7 - form 2: lm(Ozone ~ Temp:Region, data = homeworkdat)

y <- matrix(homeworkdat$Ozone, ncol = 1)
x <- matrix(homeworkdat$Temp, ncol = 1)
Z <- cbind(c(1,0,1,0,1), c(0,1,0,1,0))


### Chapter 3: Time Series introduction
#example time series
data(WWWusage, package = "datasets")
par(mai = c(0.9, 0.9, 0.1, 0.1), omi = c(0,0,0,0))
#Plot
plot.ts(WWWusage, ylab = "", las = 1, col = "blue", lwd = 2)

#second example:
data(lynx, package = "datasets")
par(mai = c(0.0, 0.9, 0.1, 0.1), omi = c(0,0,0,0))
plot.ts(lynx, ylab = "", las = 1, col = "blue", lwd  = 2)

#simple ts model of random noise: white noise:
# Xt ~ N(0, 1)
par(mai = c(0.9, 0.9, 0.1, 0.1), omi = c(0, 0, 0, 0))
matplot(ww, type = "l", lty = "solid", las = 1, ylab = expression(italic(x[t])), 
xlab = "Time", col = gray(0.5, 0.4)) 

#random walk model:
#Xt = Xt-1 + wt ; with wt ~ N(0,1)
par(mai = c(0.9, 0.9, 0.1, 0.1), omi = c(0, 0, 0, 0))
matplot(apply(ww, 2, cumsum), type = "l", lty = "solid", las = 1,
        ylab = expression(italic(x[t])), xlab = "Time", col = gray(0.5, 0.4))



## classical decomposition:
#model time sereis as combo of trend (mt), seasonal component: (st), and remainder (et)
# Xt = mt + st + et
#first, need trend (mt)
#to get trend (mt), could use moving average
#Mt = 1/(2a +1)* Xt-1 + Xt + Xt + 1 + Xt + n  ;  the 1/(2a+1) = lambda
#can change a value to add more values to average on each side for the model

#st = seasonal effect, use subtraction
#st = Xt - mt;  this st includes the et (remainder) with it
#instead, can also calc mean seasonal effect; plot
seas_2 <- decompose(xx)$seasonal
par(mai = c(0.9, 0.9, 0.1, 0.1), omi = c(0, 0, 0, 0))
plot.ts(seas_2, las = 1, ylab = "")

#getting the remainder (et)
#estimate by subtraction:
#et  = xt - mt - st
#plot error:
ee <- decompose(xx)$random
par(mai = c(0.9, 0.9, 0.1, 0.1), omi = c(0, 0, 0, 0))
plot.ts(ee, las = 1, ylab = "")

## decomposition of log transformed data:
#log transformed airline data:
lx <- log(AirPassengers)
par(mai = c(0.9, 0.9, 0.1, 0.1), omi = c(0, 0, 0, 0))
plot.ts(lx, las = 1, ylab = "")


## Chapter 4: basic time series functions:
#load data:
data(NHTemp, package = "atsalibrary")
Temp <- NHTemp
data(MLCO2, package = "atsalibrary") #CO2 data from Mauna Loa obs
CO2 <- MLCO2
data(hourlyphyto, package = "atsalibrary") #phytoplankton countdata
phyto_dat <- hourlyphyto

#R converts dataframes into time series objects, ts using ts() function
#ts() takes two args: frequency (samples/cycle), and start (or first sample timepoint)
co2 <- ts(data = CO2$ppm, frequency = 12, #monthly observations
          start = c(CO2[1, "year"], CO2[1, "month"])) #firs obs = first month / year from df

#plot data:
plot.ts(co2, ylab = expression(paste("CO"[2], " (ppm)")))

#combine and plot multiple ts together:
#first make temp ts object:
temp_ts <- ts(data = Temp$Value, frequency = 12, start = c(1880, 1))
#need to line up time series since they begin at diff times
#use ts.intersect function:
dat_int <- ts.intersect(co2, temp_ts)
dim(dat_int)
#use reg plot function to plot both stacked with alt y axes
plot(dat_int, main = "", yax.flip = TRUE)

## decomposition of time series:
#recall that modeling ts with 3 parts: Xt = mt + st + et

#1 start estimating trends, mt, using moving averages or filters (windows of the data with various length averages: 1/1+2a)
#can use the filter() function within stats::filter to estimating moving avgs and linear filters:

#first need to show filter weights:
fltr <- c(1/2, rep(1, times = 11), 1/2)/12  #this adds 1/2 to both extreme ends so that filtered value at time t matches up with original observation

#estimtae of tremd (mt):
co2_trend <- stats::filter(co2, filter = fltr, method = "convo", sides = 2) #sides = both sides of the mean

#Plot trend:
plot.ts(co2_trend, ylab = "Trend", cex = 1)

#2 estimate seasonal effects:
#use subtraction: st = xt - mt ; this would include remainder / error
co2_seasonal <- co2 - co2_trend

#plot seasonal effects:
plot.ts(co2_seasonal, ylab = "Seasonal effect", xlab = "Month", cex = 1)

#can get overall seasonal effect by averaging across 12 months:
ll <- length(co2_seasonal) #total length of seasonal data
ff <- frequency(co2_seasonal) #12 months
periods <- ll %/% ff #num periods or years

index <- seq(1, 11, by = ff)-1 #index of cumulative month
mm <- numeric(ff) #Mean of month:
for (i in 1:ff) {
  mm[i] <- mean(co2_seasonal[index + i], na.rm = TRUE)
}
#subtract mean to make overall mean = 0
mm <- mm - mean(mm)

#plot monthly seasonal effects:
plot.ts(mm, ylab = "Seasonal effect", xlab = "Month", cex = 1)

#create ts object for season
co2_seas_ts <- ts(rep(mm, periods + 1)[seq(ll)], start = start(co2_seasonal),
frequency = ff)

#complete model with error: e = Xt - mt - st
co2_err <- co2 - co2_trend - co2_seas_ts

#full plot:
plot(cbind(co2, co2_trend, co2_seas_ts, co2_err), main = "",
     yax.flip = TRUE)

## decompose function: does all the above (trend, seasonal, error) breakdown all in one
co2_decomp <- decompose(co2) #if the time series is multiplicative (i.e., not linear or changes each time unit, then can use the arg: type = multiplicative)
str(co2_decomp)   

#plot all:
plot(co2_decomp, yax.flip = TRUE)


##differencing to remove a trend or seasonal effects
#using diff() function; takes the data as first arg, the lag (lag at which to difference), and differences; e.g., first order difference removes a linear trend, second order diff removes quadratic trend, etc..
#removing seasonal trend, use lag set = to the period, eg. 1 month

#twice-diff the co2 data:
co2_d2 <- diff(co2, differences = 2) #removes the trend
plot(co2_d2, ylab = expression(paste(nabla^2, "CO"[2]))) #but seasonal effect remains

#difference the differnece data to remove seasonal by setting the lag to 12, or per period measure
co2_d2_d12 <- diff(co2_d2, lag = 12)
plot(co2_d2_d12, ylab = expression(paste(nabla, "(", nabla^2, "CO"[2])))
#now, remaining variation likely is due just to random errors


## autocorrelation Rk:
#use this function for plotting acf()
plot.acf <- function(ACFobj) {
  rr <- ACFobj$acf[-1]
  kk <- length(rr)
  nn <- ACFobj$n.used
  plot(seq(kk), rr, type = "h", lwd = 2, yaxs = "i", xaxs = "i", 
       ylim = c(floor(min(rr)), 1), xlim = c(0, kk + 1), xlab = "Lag", 
       ylab = "Correlation", las = 1)
  abline(h = -1/nn + c(-2, 2)/sqrt(nn), lty = "dashed", col = "blue")
  abline(h = 0)
}


##Simulating white noise: WN is assumed if error terms following time series model adjustment have a mean of zero and are independent, identitically distributed residuals
#first generate 100 random samples:
set.seed(123)
#random normal variates
GWN <- rnorm(n = 1e6, mean = 5, sd = 0.2)
#random poisson variates:
PWN <- rpois(n = 5e5, lambda = 20)

#plot these dists:
par(mfrow = c(1,2)) #1 row, 2 cols
plot.ts(GWN)
abline(h = 5, col = "blue", lty = "dashed")
plot.ts(PWN)
abline(h = 20, col = "blue", lty = "dashed")

#test autocorrelation (for WN series, should be zero for lags >=1)
par(mfrow = c(1, 2))
#Plot autocorrelation of normal variats with mean 
acf(GWN, main = "", lag.max = 20)
#Plot poisson autocorr with normal variates
acf(PWN, main = "", lag.max = 20) #plots show > 0, but not signif different, which with 20 samples is something we'd see by chance
#when we use much larger sample, the autocorr is essentially 0

##random walks: are the most simple non-stationary timeseries model:
#RN = a time series {Xt} where Xt = Xt-1 + wt;  where wt = discrete white noise series with independent identically distributed residuals w/ mean = 0
#usually assuming that wt ~ N(0, q), or Gaussian distributed

#simulate random walk:
set.seed(123)
TT <- 100 #Length of time series
#initialize {x_t} and {w_t}
xx <- ww <- rnorm(n = TT, mean = 0, sd = 1)

#compute values 2 : TT
for (t in 2:TT) {
  xx[t] <- xx[t-1] + ww[t]
}

#now plot the RW time series + the ACF
par(mfrow = c(1, 2))
#plot line:
plot.ts(xx, ylab = expression(italic(x[t])))
#Plot acf
plot.acf(acf(xx, plot = FALSE)) #autocorr is high even after several rounds of sampling

#alternate simulation of RW:
#can substitute Xt-2 + wt + wt-1 for the eq. Xt-1
#this results in a series as t goes from 1 to T, sum of all wt, or the sum of all random errors up to time point T

x2 <- cumsum(ww) #the sum / series equation of all random errors for each time point
#plot both time series to compare the random walk types:
par(mfrow = c(1, 2))
plot.ts(xx, ylab = expression(italic(x[t])))
plot.ts(x2, ylab = expression(italic(x[t]))) #both plots work to show RW


## Autoregressive models:
#ARIMA = Autogressive integrated, moving average




