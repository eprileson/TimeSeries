##############################################################################33
########## Project - 592 - Dros Population Time Series Modeling #########
#########################################################################
###########

setwd("C:/Users/prile/OneDrive - Washington State University (email.wsu.edu)/PhD_Documents/Courses/BIOL_592_TimeSeries/TimeSeries")

#latest version of R:
# Install and load the necessary package
install.packages("installr")
library(installr)

# Check for updates and update R if a new version is available
if (installr::updateR() == FALSE) {
  cat("R is already up to date.\n")
} else {
  cat("R has been successfully updated to the latest version.\n")
}

## Packages and libraries:
#More efficient way of checking and installing packages:
packages <- c("lubridate", "reshape2", "devtools", "learnr", "stats", "TMB", 
              "MARSS", "marssTMB", "datasets", "magrittr", "tidyr", "dplyr", "forecast", 
              "ggplot2", "viridis", "esquisse", "cowplot", "corrplot", "PerformanceAnalytics", "MCMCglmm", "MASS", "AICcmodavg")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
#load packages:
invisible(lapply(packages, library, character.only = TRUE))


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

#log transform pop data, ignore 0s
dros$log_pop <- ifelse(dros$population != 0, log(dros$population), 0)

#change date to date class:
dros$Date <- mdy(dros$Date)

#basic plot for raw data:
baseplot <- 
ggplot(data = dros, aes(x = Date, y = log_pop, colour = treatment, group = treatment)) +
  geom_line(stat = "summary", fun = mean, size = 1)+
  geom_point(size = 3, stat = "summary", fun = "mean") +
  scale_color_viridis_d(option = "plasma") +
  facet_wrap(~treatment, nrow = 2, ncol = 2)+
  labs(y = bquote("Population size (log "[e]*")"), color = "Population type") +
  guides(color = "none")+
  theme_classic()+
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 12)
  )
#geom_line()+
  
#geom_errorbar(width = 0.05, fun.data = "mean_se", stat = "summary")+
  

#now instead of grouping by cage, let's group by avg pop / population type:
#this will limit the total rows to the number of observations for each pop, so 76
drosAvg <- dros %>%
  group_by(Date, treatment) %>%
  summarise(pop_avg = mean(log_pop))

#now pivot back wider to get each Timepoint as a col
#could also use t() transpose once it's a matrix
drosAvg <- drosAvg %>%
  pivot_wider(names_from = "Date", values_from = "pop_avg")

#remove name of pop column for matrix transform
drosAvg <- drosAvg[,-c(1)]

tmp1 <- rownames(drosAvg)
drosAvg <- as.data.frame(drosAvg)
rownames(drosAvg) <- c("E", "LRC", "PA", "S")


#####
######
## Data Simulation for Model 1: Biased Random Walk ####
set.seed(592)
#use mvnorm() function from {MASS} pkg to simulate multivariate normally dist data
TT <- 19 #time steps

#set strength of density dependence (highly variable in Dros)

bb <- diag(c(0.5, 0.5, 0.5, 0.5))

#bias term for RW:
uu <- 0 

##var-cov matrix for the process error in the state model:
QQ <- diag(c(0.1, 0.1, 0.1, 0.1)) # no covariance, simulated process error vals

#transpose process errors to match equational form
ww <- mvrnorm(n = TT, mu = c(0,0,0,0), Sigma = QQ) %>% t() #mu = mean, specified to be 0 in all models where we subtract u

#initialize state vector:
xx <- ww

#sinmulate data as random walk with bias
for (t in 2:TT){
  xx[,t] <- bb %*% xx[, t-1] + uu + ww[, t]
}

#now add observation error with var-covar matrix
RR <- diag(c(0.05, 0.05, 0.05, 0.05))

#observation errors from mvn with 0 as mean and var as matrix above; transposed to match vals in state
vv <- mvrnorm(n = TT, mu = c(0,0,0,0), Sigma = RR) %>% t()

yy <- xx + vv


## this simulated data worked fine, so now on to adding in my own data from the 
#dros dataset

### setting up model parameters for MARSS: need to make everything a matrix using a list for each item
# from the state and obs models:

#simulated model list:
model_listSim <- list(
  #state model:
  B = "diagonal and equal", #identity matrix
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

YYsim <- matrix(yy, nrow = 4, ncol = 19)
simMod <- MARSS(YYsim, model_listSim)
#plot simulated data
autoplot(simModA, plot.type = "fitted.ytT")+
  theme_classic()#this seems to plot the simulated data
#well


## Now implement the real data:"
#the real data model list:
model_listA1 <- list(
  #state model:
  B = "diagonal and equal", #identity matrix
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



#make new model with params with U as 'unequal'
model_listA2 <- list(
  #state model:
  B = "diagonal and equal", #identity matrix
  U = "unequal", # 4x1 0 matrix for bias factor
  C = matrix(0, 4, 1), #no predictor for the first model
  c = matrix(0),
  Q = "unconstrained",
  #obs model
  Z = factor(c(1,2,3,3)), #Z matrix to map obs onto state
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
# T = time steps, so 19. Note this was done above in data wrangling
#remove date col first:
#YY <- matrix(dros, nrow = 4, ncol = 19)
YY <- as.matrix(drosAvg) #the matrix() was changing the df to a list: "matrix" "array" and not a matrix that R recognized
class(YY)

#fit model:
mod_rwA1 <- MARSS(y = YY, model = model_listA1, method = "BFGS") #was able to run 
summary(mod_rwA1)

#fit second model:
mod_rwA2 <- MARSS(y = YY, model = model_listA2, method = "BFGS")

#diff Q iterations
#model with Q = "diagonal and unequal" --> AICc = 290.96
#Model with Q = "unconstrained" --> AICc = 171.4151 
#Model with Q = "identity" --> AICc = 280.5394


##analyses:
xx_hatA <- mod_rwA2$ytT

#extract coef:
params <- MARSSparamCIs(mod_rwA2, alpha = 0.05)
upper <- matrix(params$ytT + ((params$ytT.se)*1.96))
lower <- matrix(params$ytT - ((params$ytT.se)*1.96))

##plot estimates:
plot(mod_rwA2)

#Final plot, model 1
autoplot(mod_rwA2, plot.type = "fitted.ytT")+ # plot for each time series 
  geom_line(color =cols, size = 1)+
  geom_point(size = 2, color = cols)+
  scale_color_viridis_d()+
  labs(title = "", caption = "", x = "Time step", y = bquote("Population (Log "[e]*")"))+
  ylim(-0.5, 11)+
  theme_classic()+
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 12)
  )  
#not able to color the CIs and line
#geom_ribbon(aes(ymin = lower, ymax = upper, fill = params$states),color = cols, alpha = 0.25)+
  
cols <- c("#0E0887","#0E0887","#0E0887","#0E0887","#0E0887","#0E0887","#0E0887","#0E0887","#0E0887","#0E0887","#0E0887","#0E0887","#0E0887","#0E0887","#0E0887","#0E0887", "#0E0887","#0E0887","#0E0887",
          "#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D","#9D179D",
          "#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953","#EB7953",
          "#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A","#F2FA1A")
          

  

##Forecasting:
  
#try the forecast package:
mod_rwA2 %>%
  forecast(h = 8) %>% #set number of time steps beyond data to 4
  autoplot()
  
  #doesn't tell us much
  
  
  



#####
######
## Data Model 2: AR(1) model with 1 covariate ####


#make a matrix for your covariate of temperature:
#exported out the original csv file to now add the temperature data:
#write.csv(dros, "dros_ts.csv")

#now with temp data added, re-import dros data:
#use the drosT df as a covariate matrix for temperature
drosT <- read.csv("dros_ts.csv", header = T)
head(drosT)

#need to group and summarize temp:
drosT <- drosT %>%
  group_by(Date) %>%
  summarise(temp = mean(temperature))

#remove nonTemp cols
drosT <- drosT[,-c(1)]
head(drosT)

#change to matrix by wide with each col representing a time stamp
drosT <- as.matrix(drosT)
drosT <- t(drosT)

#create the parameter list for the MARSS model with each parameter its own vector or matrix
model_listB <- list(
  #state model:
  B = diag(4), #identity matrix
  U = "unequal", # 4x1 matrix for bias factor to be unequal for each state
  Q = "unconstrained",
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

#create the C matrix for the covariate of temperature; add it into the 
#list for the model; C matrix is just the vector of temps, while c is the matrix of 
#temps
temp <- matrix(drosT,nrow=1)
C1 <- matrix(c("temp1","temp2", "temp3", "temp4"),4,1)
model_listB$C <- C1
model_listB$c <- temp



#fit model with the response as the observation matrix (no with the covariates included) and the model as the list of matrices 
#defined above:
mod_rwB <- MARSS(y = YY, model = model_listB, method = "BFGS")

##analyses:
summary(mod_rwB)
xx_hatB <- mod_rwB$states
xx_hatBSE <- mod_rwB$states.se

paramsB <- MARSSparamCIs(mod_rwB, alpha = 0.05)

paramsB$ytT <- as.data.frame(paramsB$ytT)

low <- matrix(paramsB$par.lowCI$x0, nrow = 4, ncol = 19)
up <- matrix(paramsB$par.upCI$x0, nrow = 4, ncol = 19)
lowCI <- paramsB$ytT - low
upCI <- paramsB$ytT + up

#visualize
plot(mod_rwB)+
  geom_ribbon()
geom_ribbon(aes(ymin= lowCI, ymax=upCI), fill="grey") +

autoplot(mod_rwB, plot.type = "fitted.ytT")+
  labs(title = "", caption = "", x = "Time", y = bquote("Population (Log"[e]*" )"))+
  theme_classic()+
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 12)
  )

# plot for each time series 
mod_rwB$states

#model selection:
AICs <- c(mod_rwA2$AICc,mod_rwB$AICc, mod_rwC$AICc)
AICs

## the first covariate model looks like it is not capturing the up and down nature
# of the data. I'll try a seasonal effect first:

#data needs to be a time series:
drosTS <- drosT #make new object for this analysis

drosTS <- ts(t(drosTS), frequency = 19)
#add seasonal effect with sine/cosine pairs:
drosSin <- forecast::fourier(drosTS, 1) %>%
  t()

#only include sin:
drosSin <- drosSin[-c(2),]

#now add this in as the new c covariate
#first new model list:

model_listC <- list(
  #state model:
  B = diag(4), #identity matrix
  U = "unequal", # 4x1 0 matrix for bias factor
  Q = "unconstrained",
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

#create the C matrix for the covariate of temperature; add it into the 
#list for the model; C matrix is just the vector of temps, while c is the matrix of 
#temps
tempS <- matrix(drosSin,nrow=1)
C2 <- matrix(c("temp1","temp2", "temp3", "temp4"),4,1)
model_listC$C <- C2
model_listC$c <- tempS

#fit model with the response as the observation matrix (no with the covariates included) and the model as the list of matrices 
#defined above:
mod_rwC <- MARSS(y = YY, model = model_listC, method = "BFGS")

##analyses:
summary(mod_rwC)
xx_hatB <- mod_rwC$states
xx_hatBSE <- mod_rwC$states.se

paramsC <- MARSSparamCIs(mod_rwC, alpha = 0.05)

#visualize
plot(mod_rwA1)

autoplot(mod_rwC, plot.type = "fitted.ytT")+ # plot for each time series 
  labs(title = "", caption = "", x = "Time", y = bquote("Population (Log"[e]*" )"))+
  theme_classic()+
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 12)
  )



### fourth model: now try with cosine as seasonal temp effect:

#add seasonal effect with sine/cosine pairs:
drosCos <- forecast::fourier(drosTS, 1) %>%
  t()

#only include sin:
drosCos <- drosCos[-c(1),]

#now add this in as the new c covariate
#first new model list:

model_listD <- list(
  #state model:
  B = diag(4), #identity matrix
  U = "unequal", # 4x1 0 matrix for bias factor
  Q = "unconstrained",
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

#create the C matrix for the covariate of temperature; add it into the 
#list for the model; C matrix is just the vector of temps, while c is the matrix of 
#temps
tempD <- matrix(drosCos,nrow=1)
C4 <- matrix(c("temp1","temp2", "temp3", "temp4"),4,1)
model_listD$C <- C4
model_listD$c <- tempD

#fit model with the response as the observation matrix (no with the covariates included) and the model as the list of matrices 
#defined above:
mod_rwD <- MARSS(y = YY, model = model_listD, method = "BFGS")


summary(mod_rwD)
paramsD <- MARSSparamCIs(mod_rwD, alpha = 0.05)
paramsD$states.se


autoplot(mod_rwD, plot.type = "fitted.ytt")+
  labs(title = "", caption = "", x = "Time", y = bquote("Population (Log "[e]*" )"))+
  theme_classic()+
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 12)
  )



###model selection:
AICs <- c(mod_rwA2$AICc,mod_rwB$AICc, mod_rwC$AICc, mod_rwD$AICc)
mod_rwA2$AICc - mod_rwB$AICc











plot(xx[4,], xx[3,]mod_rwBplot(xx[4,], xx[3,], xx[2,], xx[1],
     xlim = range(xx[1:4,], yy[1:4, ], xx_hat[1:4, ]),
     ylim = range(xx[1:4, ], yy[1:4, ], xx_hat[1:4, ]),
     pch = 16, type = "o", col = "blue",
     xlab = "Time Points", ylab = "Fly population"))
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


####
### End

####







### Unused stuff:
### simulated data:
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



#make the Z model matrix:
#start with 0s
Z.modelA <- matrix(0, 4, 4)

#then add elements that are 1 in col 1
Z.modelA[1, c(1)] <- 1
Z.modelA[2,c(3)] <- 1
Z.modelA[3, c(3)] <- 1
Z.modelA[4, c(4)] <- 1


1000010000100001000
0100001000010000100
0010000100001000010
0001000010000100001

#elements that are 1 in col 2:
Z.modelA[c(2, 6, 10, 14, 18), 2] <- 1

#elements that are 1 in col 3:
Z.modelA[c(3, 7, 11, 15, 19), 3] <- 1

#elements that are 1 in col 4:
Z.modelA[c(4, 8, 12, 16), 4] <- 1
str(Z.modelA)

