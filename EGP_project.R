##############################################################################33
########## Project - 592 - Dros Population Time Series Modeling #########
#########################################################################
###########


setwd("C:/Users/prile/OneDrive - Washington State University (email.wsu.edu)/PhD_Documents/Courses/BIOL_592_TimeSeries")

## Packages and libraries:
#More efficient way of checking and installing packages:
packages <- c("lubridate", "reshape2", "devtools", "learnr", "stats", "TMB", 
              "MARSS", "marssTMB", "datasets", "magrittr", "tidyr", "forecast", 
              "ggplot2", "viridis", "esquisse", "MASS", "AICcmodavg")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
#load packages:
invisible(lapply(packages, library, character.only = TRUE))

install.packages('marssTMB',
                 repos = c('https://atsa-es.r-universe.dev','https://cloud.r-project.org'))
library(marssTMB)

####
#### Goals: model time series of Drosophila populations outdoors

####
#### First, create the time series data
###
##data wrangling:
# read in the population data
dros <- read.csv("dros_pop.data.csv", header = T)
head(dros)
dim(dros)
names(dros)
#change col names:
colnames(dros) <- c("cage.number", "treatment", "07242023", "08012023", "08092023",    
                    "08162023", "08232023", "08302023", "09062023" , "09132023",   
                    "09202023", "09302023", "10032023", "11122023", "10192023",
                    "10282023", "11022023", "11102023", "11162023", "11232023",
                    "12042023")
  

#remove last two cols and the remaining NA rows:
dros <- dros[-c(40:999), -c(22:23)]

#make treatment a factor:
dros$treatment <- as.factor(dros$treatment)
class(dros$treatment)

#need to go from wide to long format for dates:
dros <- dros %>%
  pivot_longer(cols = c("07242023", "08012023", "08092023",    
                        "08162023", "08232023", "08302023", "09062023" , "09132023",   
                        "09202023", "09302023", "10032023", "11122023", "10192023",
                        "10282023", "11022023", "11102023", "11162023", "11232023",
                        "12042023"),
               values_to = "population", names_to = "Date") #%>%
  #pivot_wider(names_from = "treatment", values_from = "population") #Now df is grouped by each cage treatment

#change date to date class:
dros$Date <- mdy(dros$Date)
#split into month, day, and year

#basic plot:
#esquisser set up
esquisser(data = dros)

ggplot(dros) +
  aes(x = Date, y = population, colour = treatment) +
  geom_point(shape = "circle", size = 2, stat = "summary") +
  geom_errorbar(width = 0.05, stat = "summary")+
  scale_color_viridis_d(option = "plasma", direction = -1) +
  scale_y_continuous(trans = "log") +
  labs(y = "Population size (log)", color = "Treatment") +
  theme_classic()
  

#data exploration:

hist(log(dros$population))




## now create time series object from the df:
dros.ts <- ts(dros, start = dros[1, "Date"], frequency = 1)

str(dros.ts)

#plot time series:
plot.ts(dros.ts, ylab = "Population", xlab = "Date")


#####
######
## Data Simulation for Model 1: Biased Random Walk ####
set.seed(592)
library(MASS)
#use mvnorm() function from {MASS} pkg to simulate multivariate normally dist data
TT <- 19 #time steps

#set strength of density dependence (highly variable in Dros)

bb <- 0.5

#bias term for RW:
uu <- 0 

##var-cov matrix for the process error in the state model:
QQ <- diag(c(0.01, 0.01, 0.01, 0.01)) # no covariance, simulated process error vals

#transpose process errors to match equational form
ww <- mvrnorm(n = TT, mu = c(0,0,0,0), Sigma = QQ) %>% t() #mu = mean, specified to be 0 in all models where we subtract u

#initialize state vector:
xx <- ww

#sinmulate data as random walk with bias
for (t in 2:TT){
  xx[,t] <- xx[, t-1] + uu + ww[, t]
}

#now add observation error with var-covar matrix
RR <- diag(c(0.005, 0.005, 0.005, 0.005))

#observation errors from mvn with 0 as mean and var as matrix above; transposed to match vals in state
vv <- mvrnorm(n = TT, mu = c(0,0,0,0), Sigma = RR) %>% t()

yy <- xx + vv


#make the Z model matrix:
#start with 0s
Z.modelA <- matrix(0, 19, 4)

#then add elements that are 1 in col 1
Z.modelA[c(1, 5, 9, 13, 17),1] <- 1

#elements that are 1 in col 2:
Z.modelA[c(2, 6, 10, 14, 18), 2] <- 1

#elements that are 1 in col 3:
Z.modelA[c(3, 7, 11, 15, 19), 3] <- 1

#elements that are 1 in col 4:
Z.modelA[c(4, 8, 12, 16), 4] <- 1
str(Z.modelA)

### setting up model parameters for MARSS: need to make everything a matrix using a list for each item
# from the state and obs models:

model_listA <- list(
  #state model:
  B = diag(4), #identity matrix
  U = matrix(0, 4, 1), # 4x1 0 matrix for bias factor
  C = matrix(0, 4, 1), #no predictor for the first model
  c = matrix(0),
  Q = matrix(c(list("q"),list(0),list(0),list(0), #Q covariance matrix for process errors
               list(0),list("q"),list(0),list(0),
               list(0),list(0),list("q"),list(0),
               list(0),list(0),list(0),list("q")), 4, 4),
  #obs model
  Z = diag(4), #Z matrix to map obs onto state
  A = matrix(0, 4, 1), #scaling coeff A, set to 0
  D = matrix(0, 4, 1), #explanatory coeff effect none here set to 0
  d = matrix(0), #explanatory variable, none for this model
  R = matrix(c(list("r"),list(0),list(0),list(0), #Q covariance matrix for process errors
               list(0),list("r"),list(0),list(0),
               list(0),list(0),list("r"),list(0),
               list(0),list(0),list(0),list("r")), 4, 4)
  
)

##now fit the random walk with bias to the simulated data with MARSS

#need to set the observations (y) to a defined N (rows) x T (cols) matrix, where
# T = time steps, so 14:
YY <- matrix(yy, nrow = 4, ncol = TT)


#fit model:
mod_rwA <- MARSS(y = YY, model = model_listA)

##analyses:
xx_hatA <- mod_rwA$states

##plot estimates:
ggplot(data = mod_rw$states, aes(x = xx, y = xx_hat))+
  geom_point()




plot(xx[2,], xx[1],
     xlim = range(xx[1:4,], yy[1:4, ], xx_hat[1:4, ]),
     ylim = range(xx[1:4, ], yy[1:4, ], xx_hat[1:4, ]),
     pch = 16, type = "o", col = "blue",
     xlab = "Time Points", ylab = "Fly population")
## add start & end points
points(xx[4,], xx[3,], xx[2, 1], xx[1, 1],
       pch = 5, cex = 2, col = "blue")
points(xx[4, TT],xx[3, TT], xx[2, TT], xx[1, TT],
       pch = 0, cex = 2, col = "blue")
## add the obs
lines(yy[4,], yy[3,],yy[2,], yy[1,],
      type = "o", pch = 16, col = "darkgray")
points(yy[4,], yy[3,], yy[2, 1], yy[1, 1],
       pch = 5, cex = 2, col = "darkgray")
points(yy[2, TT], yy[1, TT],
       pch = 0, cex = 2, col = "darkgray")
## add the estimated states
lines(xx_hat[4,], xx_hat[3,],xx_hat[2,], xx_hat[1,],
      type = "o", pch = 16, col = "darkred")
points(xx_hat[4,1], xx_hat[3,1], xx_hat[2, 1], xx_hat[1, 1],
       pch = 5, cex = 2, col = "darkred")
points(xx_hat[4, TT], xx_hat[3, TT], xx_hat[2, TT], xx_hat[1, TT],
       pch = 0, cex = 2, col = "darkred")



#####
######
## Data Simulation for Model 2: AR(1) model with 1 covariate ####
set.seed(592)
library(MASS)
#use mvnorm() function from {MASS} pkg to simulate multivariate normally dist data
TT_b <- 19 #time steps

#set strength of density dependence (highly variable in Dros)

bb <- 0.5

#bias term for RW:
uu <- 0 

##var-cov matrix for the process error in the state model:
QQ <- diag(c(0.01, 0.01, 0.01, 0.01)) # no covariance, simulated process error vals

#transpose process errors to match equational form
ww <- mvrnorm(n = TT_b, mu = c(0,0,0,0), Sigma = QQ) %>% t() #mu = mean, specified to be 0 in all models where we subtract u

#initialize state vector:
xx <- ww

#sinmulate data as random walk with bias
for (t in 2:TT_b){
  xx[,t] <- xx[, t-1] + uu + ww[, t]
}

#now add observation error with var-covar matrix
RR <- diag(c(0.005, 0.005, 0.005, 0.005))

#observation errors from mvn with 0 as mean and var as matrix above; transposed to match vals in state
vv <- mvrnorm(n = TT_b, mu = c(0,0,0,0), Sigma = RR) %>% t()

yy <- xx + vv


### setting up model parameters for MARSS: need to make everything a matrix using a list for each item
# from the state and obs models:

#make the Z model matrix:
#start with 0s
Z.modelB <- matrix(0, 19, 4)

#then add elements that are 1 in col 1
Z.modelB[c(1, 5, 9, 13, 17),1] <- 1

#elements that are 1 in col 2:
Z.modelB[c(2, 6, 10, 14, 18), 2] <- 1

#elements that are 1 in col 3:
Z.modelB[c(3, 7, 11, 15, 19), 3] <- 1

#elements that are 1 in col 4:
Z.modelB[c(4, 8, 12, 16), 4] <- 1

model_listB <- list(
  #state model:
  B = diag(4), #identity matrix
  U = matrix(0, 4, 1), # 4x1 0 matrix for bias factor
  C = matrix(0.5, 4, 1), #population predictor strength
  c = matrix(1, 1, 1),
  Q = matrix(c(list("q"),list(0),list(0),list(0), #Q covariance matrix for process errors
               list(0),list("q"),list(0),list(0),
               list(0),list(0),list("q"),list(0),
               list(0),list(0),list(0),list("q")), 4, 4),
  #obs model
  Z = diag(4),
  A = matrix(0, 4, 1), #scaling coeff A, set to 0
  D = matrix(0, 4, 1), #explanatory coeff effect none here set to 0
  d = matrix(0), #explanatory variable, none for this model
  R = matrix(c(list("r"),list(0),list(0),list(0), #Q covariance matrix for process errors
               list(0),list("r"),list(0),list(0),
               list(0),list(0),list("r"),list(0),
               list(0),list(0),list(0),list("r")), 4, 4)
  
)


#need to set the observations (y) to a defined N (rows) x T (cols) matrix, where
# T = time steps, so 14:
YY_b <- matrix(yy, nrow = 4, ncol = TT_b)

##now fit the random walk with bias to the simulated data with MARSS

#fit model with the response as the observation matrix and the model as the list of matrices 
#defined above:
mod_rwB <- MARSS(y = YY_b, model = model_listB)

##analyses:
xx_hatB <- mod_rwB$states
  

#model selection:
lapply(mods, AICc)




plot(xx[4,], xx[3,], xx[2,], xx[1],
     xlim = range(xx[1:4,], yy[1:4, ], xx_hat[1:4, ]),
     ylim = range(xx[1:4, ], yy[1:4, ], xx_hat[1:4, ]),
     pch = 16, type = "o", col = "blue",
     xlab = "Time Points", ylab = "Fly population")
## add start & end points
points(xx[4,], xx[3,], xx[2, 1], xx[1, 1],
       pch = 5, cex = 2, col = "blue")
points(xx[4, TT],xx[3, TT], xx[2, TT], xx[1, TT],
       pch = 0, cex = 2, col = "blue")
## add the obs
lines(yy[4,], yy[3,],yy[2,], yy[1,],
      type = "o", pch = 16, col = "darkgray")
points(yy[4,], yy[3,], yy[2, 1], yy[1, 1],
       pch = 5, cex = 2, col = "darkgray")
points(yy[2, TT], yy[1, TT],
       pch = 0, cex = 2, col = "darkgray")
## add the estimated states
lines(xx_hat[4,], xx_hat[3,],xx_hat[2,], xx_hat[1,],
      type = "o", pch = 16, col = "darkred")
points(xx_hat[4,1], xx_hat[3,1], xx_hat[2, 1], xx_hat[1, 1],
       pch = 5, cex = 2, col = "darkred")
points(xx_hat[4, TT], xx_hat[3, TT], xx_hat[2, TT], xx_hat[1, TT],
       pch = 0, cex = 2, col = "darkred")


