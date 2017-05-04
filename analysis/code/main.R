############################################################
### ECON 293 Homework 1

# Luis Armona and Jack Blundell
# Stanford University
# Spring 2017


# Section numbers correspond to assignment page

############################################################

# set your working directory

setwd("C:/Users/Jack/Documents/Git/Athey ML homework 1/AtheyMLhw1") # Jack
#setwd('/home/luis/AtheyMLhw1') #Luis
# clear things in RStudio

rm(list = ls())


# Call packages

library(ggplot2)
library(dplyr)
library(reshape2)
library(glmnet)
library(plotmo)
library(balanceHD)

# set seed
set.seed(12345)

############################################################

# Load data
fname <- 'analysis/input/charitable_withdummyvariables.csv'
char <- read.csv(fname)
attach(char) # attach so don't have to call each time


### Exploratory analysis

dim(char) # 50,083 obs, 63 vars
names(char)
head(char) # Look at first few entries of each var

# Treatment
summary(treatment) # Anyone who got any of the 27 treatments (3 match x 3 match size x 3 reccomended amount)
mean(treatment) # 67% treated


# Gives at all
summary(out_gavedum)

# Giving
summary(out_amountgive) # amount given. Highly skewed
hist(out_amountgive)

############################################################
### 1. Regression for average treatment effect

reg.ols <- lm(out_amountgive ~ treatment) 
summary(reg.ols) # show results, significant at 90% but not 95% level
# Consistent with Table 4 of paper
confint(reg.ols, level=0.95) # CI

#probit regression
gave.probit <- glm(out_gavedum ~ treatment,family=binomial(link='probit'))
#convert coef to derivative
marginal.effect <- mean(dnorm(predict(gave_probit, type = "link")))*coef(gave_probit)
print(marginal.effect)

############################################################
### 2. Dropping some observations

# Try to find some drop of observations based on observables that changes the treatment effect

# look at some potential variables

summary(page18_39)
hist(page18_39[page18_39_missing!=1])

summary(perbush)
hist(perbush[perbush_missing!=1])

# Generate restricted dataset, dropping all those missing key covariates
char.res <- char[ which(page18_39!=-999
                         & perbush!=-999
                         & median_hhincome!=-999), ] # drop all those with missings of key variables
detach(char)
attach(char.res) # attach so don't have to call each time

#char_res$drop <- 0 # variable telling us to drop or not

# Make threshold rule for dropping (alternatively do with random variable)
#char_res$thres <- perbush #+ 0.1*(1+perbush)^2 - 0.1*page18_39*perbush #- page18_39 - page18_39^2 + perbush^2
#summary(char_res$thres)
#char_res$drop[char_res$thres <= 0.5] <- 1
#mean(char_res$drop) # drop 23 % of obs
#char_res_d <- char_res[which(char_res$drop == 0),]

##############################
#Alternative rule: 
#randomly censor individuals
#via a  complex, highly nonlinear fcn  of votes 4 bush in state,
#

ps.fcn <- function(v,c,pg,t){
  #v_t <- (v-.25)/.5
  v_t <- v
  ihs_pg <- log(pg + sqrt(pg ^ 2 + 1))/5
  
  #p<- ((2*c*(t*acos(v_t)) + (1-t)*atan(v_t^2))  - .5*exp(v_t) + t*((ihs_i)^4)/4 + (1-t)*(i/10000))/4
  p<- (c*(acos(v_t))*atan(v_t^2)  - .5*exp(v_t))/4 + (t*((ihs_pg)) + (1-t))/2
  p<- pmin(pmax(0,p),1)
  return(p)
}
#story to accompany this fcn: ACLU wants to help those in trouble in "red states" but do not 
#feel they can make a difference in really, really red states so target donors less often
plot(seq(0,1,.001),ps.fcn(seq(0,1,.001),2,800,1),ylim=c(0,1))#a plot of the function
lines(seq(0,1,.001),ps.fcn(seq(0,1,.001),1,800,0))
lines(seq(0,1,.001),ps.fcn(seq(0,1,.001),3,200,1))

#char$mibush=char$perbush==-999
#char$perbush[char$mibush]=.5

# Input from highly non-linear function
char.res$ps.true <- ps.fcn(char.res$perbush,char.res$cases,char.res$hpa,char.res$treatment) # hpa is highest previous contribution. cases is court cases from state which organization was involved.

# Plot CDF of this nonlinear function
ggplot(char.res,aes(x=ps.true))+ stat_ecdf()

# Set seed
set.seed(21) 

# Selection rule (=1 of uniform random [0,1] is lower, so those with higher ps.true more likely to be selected)
selection <- runif(nrow(char.res)) <= char.res$ps.true

char.censored <- char.res[selection,] #remove observations via propensity score rule

ggplot(char.res,aes(x=perbush)) + geom_histogram()+xlim(c(0,1))

ggplot(char.censored,aes(x=ps.true)) + geom_histogram() +xlim(c(0,1))

#overlap in true propensity score
ggplot(char.censored,aes(x=ps.true,colour=factor(treatment))) + stat_ecdf()

ggplot(char.censored,aes(x=ps.true,y=hpa,colour=factor(treatment))) + geom_point()
#there is clear overlap, but clearly assymetries going on with hpa as well

#################################

# New regression results with dropping

#jack's threshold rule
#reg_ols_drop <- lm(out_amountgive ~ treatment, data = char_res_d) 
#summary(reg_ols_drop) 

#Luis' PS generating rule
reg.censored <- lm(out_amountgive ~ treatment, data = char.censored) 
summary(reg.censored) 

# Old regression results (remember to drop missings to make comparable sample)
reg.ols.comp <- lm(out_amountgive ~ treatment, data = char.res) 
summary(reg.ols.comp) 
ate.true <- reg.ols.comp$coefficients[2]

# Check overlap (propensity score)

# estimate propensity score via logit regression
ps.mod <- glm(treatment ~ page18_39 + perbush + pwhite + pblack + median_hhincome + red0 + hpa + ltmedmra + freq + years + year5 + dormant + female + couple + nonlit,
              family = binomial(), data = char.censored)
summary(ps.mod)

# Put propensity scores into dataframe with actual treatment
ps.df <- data.frame(pr.score = predict(ps.mod),
                     treatment = ps.mod$model$treatment)
head(ps.df)
ps.df$pr.score <- pmin(ps.df$pr.score,1)

# Generate graph (see https://stanford.edu/~ejdemyr/r-tutorials-archive/tutorial8.html#propensity-score-estimation)
labs <- paste("Actual treatment:", c("Treated", "Control"))
ps.df %>%
  mutate(treatment = ifelse(treatment == 1, labs[1], labs[2])) %>% 
  ggplot(aes(x = pr.score)) +
  geom_histogram(color = "white") +
  facet_wrap(~treatment) +
  xlab("Probability of treatment") +
  theme_bw()

#################################
#Bias function under luis' PS rule
#since we have some continuous covariates,
#use the Mahalanobis distance to get conditional means
#in a 'neighborhood' of each set of X's

bias.fcn <- function(ps,treat,y, covars){
  x<- covars
  #stdize covariates to be z-scores mean 0 sd 1
  for (j in 1:ncol(covars)){
    x[,j] = (x[,j] - mean(x[,j]) )/sd(x[,j])
  }
  covx<- cov(x)
  bias <- matrix(NaN,nrow=length(y),ncol=1)
  mu.t <- mean(y[treat])
  mu.c <- mean(y[!treat])
  p <- mean(treat)
  for (i in 1:length(y)){
    if (i%%1000==0){
      print(i)
    }
    maxdist<- .1 #the max Mahalanobis distance
    #to compute conditional means for bias function
    distances <- mahalanobis(x,center=x[i],cov<-covx)
    while (length(y[distances <= maxdist & treat]) == 0 | 
           length(y[distances <= maxdist & !treat]) == 0){
        maxdist <- maxdist+.1
    }
    mu.t.X = mean(y[distances <= maxdist & treat])
    mu.c.X = mean(y[distances <= maxdist & !treat])
    bias[i] <- (ps[i]-p)*( p*(mu.c.X-mu.c) + (1-p)*(mu.t.X-mu.t)  )
  }
  
  return(bias)
}
covars.ps <- cbind(char.censored$hpa,char.censored$cases,char.censored$perbush) #Xs relevant for p-score
char.censored$bias <- bias.fcn(char.censored$ps.true,char.censored$treatment,
                 char.censored$out_amountgive, covars.ps)
ggplot(char.censored,aes(x=bias)) +geom_histogram(fill=I("white"),col=I("black"))
E.bias = mean(char.censored$bias)/mean(char.censored$treatment)*(1-mean(char.censored$treatment))
print(E.bias)

############################################################
### 3.


# propensity score weighting ATE
#first, estimate the propensity score with a logit regression using all covars
#since in this exercise we should not know the "ground truth" propensity score
covars.all <- char.censored[,c(14:22,23:63)] #skip the state indicator used for summ stats
ps.formula <- paste('treatment ~ ', paste(names(covars.all),collapse='+'),sep='')
m.ps <- glm(ps.formula,
              family = binomial(), data = char.censored)
char.censored$ps.est <- predict(m.ps,type='response')
summary(char.censored$ps.est)
#compare estimated p-score w/ real p-score
ggplot(melt(char.censored[,c('ps.true','ps.est')]),aes(x=value,colour=variable)) + geom_density(alpha=.2)

char.censored$w.ate[char.censored$treatment == 1] <-  1/char.censored$ps.est[char.censored$treatment == 1]
char.censored$w.ate[char.censored$treatment == 0] <-  ( 1 / (1 - char.censored$ps.est[char.censored$treatment == 0]))

#regular propensity score weighting
ate.ps <- mean(char.censored$out_amountgive[char.censored$treatment==1]*
                 char.censored$w.ate[char.censored$treatment == 1]) - 
          mean(char.censored$out_amountgive[char.censored$treatment==0]*
                 char.censored$w.ate[char.censored$treatment == 0]) 
print(ate.ps)
#gives a negative score!

# direct regression analysis ATE;
# control for Xs linearly
ols.formula <- paste('out_amountgive ~ treatment +', paste(names(covars.all),collapse='+'),sep='')
reg.ols <- lm(ols.formula, data=char.censored)
print(reg.ols$coefficients['treatment'])
#does much better; pretty close to true ATE


# traditional double robust analysis weighting using inverse propensity score weighting; the lm command in R has a weights option.
#changed this to sqrt(w) since lm automatically squares them
pweight.reg <- lm(ols.formula, weights = w.ate, data = char.censored)
summary(pweight.reg)
print(pweight.reg$coefficients['treatment'])
#gives a similar answer to non-weighted reg


# lasso or regularized logistic regression (optionally try CART or randomforest --classification trees, or classification forests), to estimate the propensity scoreand re-estimate the ATE using the methods above.
#re-estimate p-score using cross validation
ps.m.cv <- cv.glmnet(as.matrix(covars.all),char.censored$treatment, family='binomial')
coef(ps.m.cv,s='lambda.min')
char.censored$ps.lasso <- predict(ps.m.cv,as.matrix(covars.all),s='lambda.min',type='response')
#compare the method's generated p-scores
ggplot(melt(char.censored[,c('ps.true','ps.est','ps.lasso')]),aes(x=value,colour=variable)) + geom_density(alpha=.2)
#we get a worse p-score since it is regularized/shrunken so more concentrated

#redo above methods (weighted mean, DR weights)
char.censored$w.ate.lasso[char.censored$treatment == 1] <-  1/char.censored$ps.lasso[char.censored$treatment == 1]
char.censored$w.ate.lasso[char.censored$treatment == 0] <-  ( 1 / (1 - char.censored$ps.lasso[char.censored$treatment == 0]))

#propensity score weighting
ate.ps.lasso <- mean(char.censored$out_amountgive[char.censored$treatment==1]*
                 char.censored$w.ate.lasso[char.censored$treatment == 1]) - 
  mean(char.censored$out_amountgive[char.censored$treatment==0]*
         char.censored$w.ate.lasso[char.censored$treatment == 0]) 
print(ate.ps.lasso)
#still bad

#DR IPR weighting
pweight.lasso.reg <- lm(ols.formula, weights = w.ate.lasso, data = char.censored)
summary(pweight.reg)
print(pweight.lasso.reg$coefficients['treatment'])
#no gain from DR method


# a single-equation lasso of Y on X and W to estimate the ATE. Note that there is an option to not penalize the treatment indicator. You may
#wish to use that option anyway so that the treatment effect estimate is not
#shrunk. See http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html for
#the syntax for setting penalties for some coefficients to 0.

#straight from the website
#set penalty coef for treatment to zero
p.fac = rep(1, ncol(covars.all)+1)
p.fac[1]=0
lasso.reg <- cv.glmnet(as.matrix(cbind(char.censored$treatment,covars.all)),char.censored$out_amountgive,penalty.factor = p.fac)
lasso.coef<-coef(lasso.reg,s='lambda.min')
print(lasso.coef['char.censored$treatment',])
#closer to zero, but actually the non-lassoed OLS estimate performed better




#Next try using the Belloni-Chernozhukov-Hansen method, where you use the lasso to select variables with non-zero coefficients from the propensity equation and the
#outcome equation, take the union of the two sets, and finally run OLS. You can either
#use cross-validation to choose lambda in each case, or you can follow the approaches
#suggested by BCH (those are more complicated but probably doesn't make a
#                  difference).

#use CV to get union of two sets
#already have p-score regularized estimation
psreg.vars<-rownames(coef(ps.m.cv,s='lambda.min'))
psreg.vars<-psreg.vars[as.logical(coef(ps.m.cv,s='lambda.min')!=0)]
#reduced form outcome reg
yreg.rf.cv<-cv.glmnet(as.matrix(covars.all),char.censored$out_amountgive)
yreg.vars<-rownames(coef(yreg.rf.cv,s='lambda.min'))
yreg.vars<-yreg.vars[as.logical(coef(yreg.rf.cv,s='lambda.min')!=0)]
ds.vars <- union(yreg.vars,psreg.vars)
ds.vars <- ds.vars[-1]#remove intercept
paste('treatment ~ ', paste(names(covars.all),collapse='+'),sep='')
doubleselect.formula <- paste('out_amountgive ~ treatment + ',paste(ds.vars,collapse='+'),sep='')
reg.DoubleSelection <- lm(doubleselect.formula,data=char.censored)
summary(reg.DoubleSelection)
print(reg.DoubleSelection$coefficients['treatment'])
#fairly close to true ATE; still only about as good as regular OLS


# Look at how your ATE coefficient changes with regularization
# Consider the single-equation LASSO of Y on X and W. Similar to the plot
# included in the handout, plot how the coefficient on the treatment indicator
# changes with lambda. Interpret your findings. (See http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html for some sample code on plotting.)

# Lasso of outcome on treatment and covars


p.fac = rep(1, ncol(covars.all)+1)
p.fac[1]=0
lasso.reg.pen <- glmnet(as.matrix(cbind(char.censored$treatment,covars.all)),char.censored$out_amountgive,penalty.factor = p.fac)

# Plot using default plot
plot(lasso.reg.pen, label = T, xvar = "lambda")

# Plot using plotmo (more options)
plot_glmnet(lasso.reg.pen,
            xvar = c("lambda"),
            label = 10, nresponse = NA, grid.col = NA, s = NA)

# No option to plot path of single coefficient, so make plot ourselves

lambdas <- seq(from = 0, to = 1, by = 0.01) # Range of lambdas to try
t.coef <- coef(lasso.reg.pen,lambdas)[2,] # Pick out treatment coefficient
ggplot(data = NULL,aes(x=lambdas, y=t.coef)) +
  geom_line(colour="blue", size=1.5) +
  xlab("Lambda") + ylab("Treatment coefficient")

# As lambda increases, treatment effect increases as expected. Treatment effect ranges from that of full OLS 
# controlling for covariates to simple OLS with no controls.

############################################################
### 4.

# double machine learning to estimate the ATE. Specifically, use a LASSO or random forest to estimate regressions of Y on X and separately Y on W. Then, run a residual on
# residual regression.

# JB: Athey has made a typo here? Should be Y on X and W on X

# LASSO of Y on X

lasso.YX <- cv.glmnet(as.matrix(covars.all),char.censored$out_amountgive)
coef(lasso.YX, s = "lambda.min")
lasso.YX.res <- predict(lasso.YX ,as.matrix(covars.all),s=lasso.YX$lambda.min) - char.censored$out_amountgive # residuals

summary(lasso.YX.res) # very skewed

# LASSO of W on X 

lasso.WX <- cv.glmnet(as.matrix(covars.all),char.censored$treatment)
coef(lasso.WX, s = "lambda.min")
lasso.WX.res <- predict(lasso.WX ,as.matrix(covars.all),s=lasso.WX$lambda.min) - char.censored$treatment # residuals
summary(lasso.WX.res) # very skewed

# Residual on residual regression

reg.res <- lm(lasso.YX.res ~ lasso.WX.res)
summary(reg.res)


# Use approximate residual balancing (package: http://github.com/swager/balanceHD) to estimate ATE
# JB: Issue with POGS installation. Haven't managed to run. http://foges.github.io/pogs/stp/r
# POGS: https://stanford.edu/class/ee364b/projects/2014projects/reports/fougner_report.pdf

tau.hat = residualBalance.ate(as.matrix(covars.all), char.censored$out_amountgive, char.censored$treatment, estimate.se = TRUE,  optimizer = "pogs")
print(paste("true tau:", tau))
print(paste("point estimate:", round(tau.hat[1], 2)))
print(paste0("95% CI for tau: (", round(tau.hat[1] - 1.96 * tau.hat[2], 2), ", ", round(tau.hat[1] + 1.96 * tau.hat[2], 2), ")"))

# Compare and interpret your results