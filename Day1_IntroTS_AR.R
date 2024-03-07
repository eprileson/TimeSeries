###################################################################################
######################## Time Course Series - WSU Pullman ####################
###############################################################################

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


######### 
## Day 1: 
### First lab section ##

#datasets:
## Atmospheric CO2 measured on Mauna Loa, Hawai'i
CO2 <- read.csv("https://raw.githubusercontent.com/SOE592/website/main/lectures/day_01/data/ML_CO2.csv")
head(CO2)
## Northern hemisphere temperature anomolies
NH_temp <- read.csv("https://raw.githubusercontent.com/SOE592/website/main/lectures/day_01/data/NH_temp.csv")

#make ts objects of datasets:
co2 <- ts(CO2$ppm, start = c(CO2[1, "year"], CO2[1, "month"]), frequency = 12)
temp <- ts(NH_temp$Value, start = c(1880, 1), frequency = 12)


#line up time indices; use intersect and union:
dat_int <- ts.intersect(co2, temp)
dim(dat_int)
str(dat_int)

#try union:
dat_union <- ts.union(co2, temp)
dim(dat_union) #Larger, as it includes all years,, but first 70 years)
str(dat_union)
head(dat_union)

#plot both ts
plot.ts(temp)
plot.ts(co2, ylab = expression(paste("CO"[2], " (ppm)")), col ="blue")

#plot both (the intersection):
plot.ts(dat_int)

#basic plot: with no title and y axes flipped for readability
plot(dat_int, main = "", yax.flip = T)


## decomposition of ts:
#1. estimate trend mt
#use filter function to get moving avg:
#create filter weight (i.e. lambda) for monthly avgs acros year
filtr <- c(1/2, #for first half for month 1
           rep(1, times = 11), #create replicates 11 times
           1/2) / 12 #for 2nd half for month 13 to get avg, integer

#use filter function to get trend using moving avg calculated above:

mt <- filter(co2, filter = filtr, method = "convo", sides = 2)
#plot trend:
plot.ts(mt)

#2 estimate seasonal effects; st:
# basic subtraction; still contains seasonal effect
st <- co2 - mt

#plot:
plot(st, ylab = "seasonal effect plus error", cex = 1)

#now calculate the mean seasonal effect:
## length of ts
ll <- length(st)

## frequency (ie, 12)
ff <- frequency(st)

## number of periods (years); %/% is integer division
periods <- ll %/% ff

## index of cumulative month
index <- seq(1, ll, by = ff) - 1

## get mean by month
mm <- numeric(ff)
for (i in 1:ff) {
  mm[i] <- mean(st[index + i], na.rm = TRUE)
}

## subtract mean to make overall mean = 0
st_mean <- mm - mean(mm)

#plot the avg seasonal:
plot.ts (st_mean, ylab = "Seasonal effect", xlab = "Month", cex = 1)

#now repeat the mean average seasonal effect for all years

st_avg <- rep(st_mean, periods + 1)[seq(ll)]

#now turn the replicated avg into ts:
st_avgTS <- ts(st_avg, start = start(st), frequency = ff)
              
#plot:
plot(st_avgTS)

#3: now get remainder /  random errors
et <- co2 - mt - st_avgTS

#plot all four together in one:
plot(cbind(co2, mt, st_avgTS, et), main = "", yax.flip = TRUE)

## now do same thing but in one function:
co2_decomp <- decompose(co2)
plot(co2_decomp, yax.flip = T)


## differencing: to remove individual parts from the time series:
#using differences arg set to 2, takes away both seas and trend
co2_diff <- diff(co2, differences = 2)
plot(co2_diff, ylab = expression(paste(nabla^2, "CO"[2])))

#this actually still has the seas component, need to set lags:
#difference the differenced time series:
co2_diff2 <- diff(co2_diff, lag = 12)
plot(co2_diff2, ylab = expression(paste(nabla^2, "CO"[2])))


## Autocorrelation: relationship of variable to itself lagged at interal; autocovariance standardized
#estimate autocorr for co2 data out to 36 months:
acf(co2, lag.max = 36) #atomatically creates plot

#estimate acf for straight line:
nn <- 100 #create number
tt <- seq(nn) #create list of numbers, i.e. a time series that increase by 1 (i.e. slope of 1)
#plot both line and acf
par(mfrow = 1:2)
plot.ts(tt, ylab = expression(italic(x[t])))
acf(tt)

##sinoidal pattern and acf
# create sine wave
tt <- sin(2 * pi * seq(nn) / 12)

## set up plot area
par(mfrow = c(1, 2))

## plot line
plot.ts(tt, ylab = expression(italic(x[t])))

## get ACF
acf(tt)

## sine or cosine waves can be used to estimate seasonal effects
## create sine wave with trend:
sst <- sin(2*pi*seq(nn)/12) - seq(nn)/50
par(mfrow = c(1:2))
    
plot.ts(sst, ylab = expression(italic(x[t])))
acf(sst)    

## partial autocorrelation:
#takes auto corr of xt and xt-1 with time points in between removed
par(mfrow = c(1, 1))
pacf(co2, lag.max = 36)


##cross correlation: linear relationship between datasets and their relationships 
#first see how data intersect:
both <- ts.intersect(sunspot.year, lynx)
dimnames(both)[[2]][1] <- "sunspots" #changes name to sunspots

plot.ts(both, yax.flip = T, main = "")


#compute ccf:
sunspots <- both[, "sunspots"]
lynx <- both[, "lynx"]

ccf(sunspots, log(lynx), ylab = "cross-correlation", main = "") #lynx on y axis, log scale;
#shows that lynx number are much lower 3 - 5 yrs after high sunspot activity


### #############
### Lab 2 ####
######

##simulating data: start with white noise
#White noise:
## set the seed for the random number generator so we all get the same results
set.seed(592)

## length of the time series
nn <- 100

## set the variance = 2 --> SD = sqrt(2)
sigma <- sqrt(2)

## draw random values
ww <- rnorm(n = nn, mean = 0, sd = sigma)

#Plot white noise:
plot.ts(ww, ylab = expression(italic(w)[t]))

#use pipe operator to do same thing:
ww %>%
  plot.ts(ylab = expression(italic(w)[t]))

#estimate correlation of ts with time shifted version of itself, ie. autocorrelation:
acf(ww) #lag 5 shows a signif deviation, but this is spurious because by chance we'd expect 1/20 points to be a false neg given alpha threshold of 0.05

#repeat with time series up to 1000
set.seed(592)

## length of the time series
nn <- 100

## set the variance = 2 --> SD = sqrt(2)
sigma <- sqrt(2)

## draw random values
ww <- rnorm(n = 100, mean = 0, sd = sigma)
ww %>%
  plot.ts(ylab = expression(italic(w)[t]))
acf(ww) 

##Random walks:
#sim random walk:

#initialize vector for storing time series to be equal to white noise term:
xx <- ww
#specify initial value:
xx[1] <- ww[1]

#Loop over data in time steps 2 - 100:
for (t in 2:nn){
  xx[t] <- xx[t-1] + ww[t]
}  #this creates the random walk with formula: xt = xt-1 + wt

#now plot the random walk:
plot.ts(x = xx, ylab = expression(italic(x)[t]))

#acf of random walk:
acf(xx) #slow decrease over lag times

## Alternate Random walk method: RW = cumulative sum of white noise sequence
wn <- cumsum(ww)

#check correlation
cor(xx, wn) #1:1, exactly the same


## Biased random walks:
#simulate biased random walk:
#first create set to hold data that is the same
xb <- ww
#set bias = 1
uu <- 0.9
# set initial value:
xb[1] <- uu + ww[1]

#loop over all values to simulate rw:
for (t in 2:100){
  xb[t] <- xb[t-1] + uu + ww[t]
}

#Plot
plot.ts(xb, ylab = expression(italic(x)[t]))

#acf for bias RW:
acf(xb) #slow descent of ac, similar to linear trend; lasts longer than acf for simple RW

## Autoregressive models

#simulating stationary AR(p) model where phi is < abs(1)
#start with AR(1) model and small phi coefficient
arsm <- list(order = c(1, 0, 0), ar = 0.1)

#now create ar model with large phi:
arlg <- list(order = c(1, 0, 0), ar = 0.9)

#now simulate the AR model outputs:
armaSM <- arima.sim(model = arsm, n = 50, sd = 0.1)
armaLG <- arima.sim(model = arlg, n = 50, sd = 0.1)

#Plot both model simulations with small and large coef
par(mfrow = c(1, 2))
#create y limits to match the two plots:
ylm <- c(min(armaSM, armaLG), max(armaSM, armaLG))

plot.ts(armaSM, ylim = ylm, ylab = expression(italic(x)[t]), main = "Small Phi", las = 1)
plot.ts(armaLG, ylim = ylm, ylab = expression(italic(x)[t]), main = "Large Phi", las = 1)
#as model approaches 0 (small phi), the model resembles white noise = stationary in mean and variabce; as model approaches 1 (large phi), model resembles Random walk

##Try a different AR1 models with same magnitude coefficient but diff signs:
arneg <- list(order = c(1, 0, 0), ar = -0.5)
arpos <- list(order = c(1, 0, 0), ar = 0.5)
#simulate
armaneg <- arima.sim(model = arneg, n = 50, sd = 0.1)
armapos <- arima.sim(model = arpos, n = 50, sd = 0.1)

#plot:
plot.ts(armaneg, ylab = expression(italic(x)[t]), main = "Negative Phi", las = 1)
plot.ts(armapos, ylab = expression(italic(x)[t]), main = "Positive Phi", las = 1)
#seems negative phi is choppier, oscillating more quickly

## Correlation of AR(p) processes
#create 4 diff AR(p) models of differing orders:
ar_p_val <- c(0.2, -0.3, -0.1, 0.7)
ar_mods <- list() #create empty list to store vals:

#now create for loop to create each model, each with its own number of orders (ie. first model is first order, 2nd 2 orders, etc.)

for (p in 1:4){
  ar_mods[[p]] <- arima.sim(n = 10000, list(ar = ar_p_val[1:p])) #use each number in values to simulate diff AR models
}

#now plot:
par(mfrow = c(1, 1))

for (p in 1:4){
  plot.ts(ar_mods[[p]][1:50], ylab = paste("AR(", p, ")", sep = " "))
  acf(ar_mods[[p]], lag.max = 12)
  pacf(ar_mods[[p]], lag.max = 12, ylab = "PACF")
}
  #PACF identifies the order of the AR model (lag = order)

## Moving average models:

#MA(q) = weighted sum of white noise random error + q most recent errors

##simulating moving average processes:
## list description for MA(1) model with small coef
MA_sm <- list(order = c(0, 0, 1), ma = 0.2)

## list description for MA(1) model with large coef
MA_lg <- list(order = c(0, 0, 1), ma = 0.8)

## list description for MA(1) model with negative coef
MA_neg <- list(order = c(0, 0, 1), ma = -0.5)

## simulate MA(1)
MA1_sm <- arima.sim(n = 50, model = MA_sm, sd = 0.1)
MA1_lg <- arima.sim(n = 50, model = MA_lg, sd = 0.1)
MA1_neg <- arima.sim(n = 50, model = MA_neg, sd = 0.1)

#plot:

par(mfrow = c(1, 3))

plot.ts(MA1_sm, ylab = expression(italic(x)[italic(t)]), main = "Small coefficient")
plot.ts(MA1_lg, ylab = expression(italic(x)[italic(t)]), main = "Large coefficient")
plot.ts(MA1_neg, ylab = expression(italic(x)[italic(t)]), main = "Negative coefficient")
#plots aren't too different as MA models are variations on white noise

#correlation ofMA processes:

#set up as the AR:
ma_q_val <- c(0.7, 0.2, -0.1, -0.3)
MA_mods <- list()

for (q in 1:4){
  MA_mods[[q]] <- arima.sim(n = 10000, list(ma = ma_q_val[1:q]))
}

#now plot with acfs and pacfs to look for correlational structure:
par(mfrow = c(1,1))

for (q in 1:4){
  plot.ts(MA_mods[[q]][1:50], ylab = paste("MA (", q, ")", sep = ""))
  acf(MA_mods[[q]], lag.max = 12)
  pacf(MA_mods[[q]], lag.max = 12)
}
#acf goes to 0 for lags greater than q, but not for PACF

## Fitting ARMA(p,q) models:


#simulate an ARMA (2, 2) model:
arm2 <- list(order = c(2, 0, 2), ar = c(0.2, -0.7), ma = c(0.2, 0.7))

#set mean:
mu <- 5

#simulate process + the mean to get a univariate time series
armasim2 <- arima.sim(n = 100, model = arm2) + mu

#estimate parameters with arima:
arima(x = armasim2, order = c(2, 0, 2))



#searching over model orders:
#write a script to loop over various orders of ARMA(p,q) models:

## empty list to store model fits
ARMA_res <- list()

## set counter for model index
cc <- 1

## loop over AR
for (p in 0:3) {
  ## loop over MA
  for (q in 0:3) {
    ARMA_res[[cc]] <- arima(x = armasim2, order = c(p, 0, q))
    cc <- cc + 1
  }
}
#now there are multiple models with various orders that you can compare 
#use AIC to compare:
## get AIC values for model evaluation
ARMA_AIC <- sapply(ARMA_res, function(x) x$aic)

## model with lowest AIC is the best
ARMA_res[[which(ARMA_AIC == min(ARMA_AIC))]]

#can do this using the auto.arima() function:

## find best ARMA(p,q) model
forecast::auto.arima(armasim2, start.p = 0, max.p = 3, start.q = 0, max.q = 3)




  






