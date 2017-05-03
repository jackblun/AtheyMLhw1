############################################################
### ECON 293 Homework 1

# Luis Armona and Jack Blundell
# Stanford University
# Spring 2017


# Section numbers correspond to assignment page

############################################################

# set your working directory

setwd("C:\\Users\\Jack\\Documents\\Git\\Athey ML homework 1\\AtheyMLhw1") # Jack
#setwd('/home/luis/Downloads/AtheyMLhw1') #Luis
# clear things in RStudio

rm(list = ls())


# Call packages

library(ggplot2)
library(dplyr)

# set seed
set.seed(12345)

############################################################

# Load data
#let's use / instead of \\ since it should be compatible with both systems (Unix/Windows) #JB: Noted! (I'm new to PC..!)
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
hist(out_amountgive[out_amountgive<=10])

############################################################
### 1. Regression for average treatment effect

reg_ols <- lm(out_amountgive ~ treatment) 
summary(reg_ols) # show results, significant at 90% but not 95% level
# Consistent with Table 4 of paper
confint(reg_ols, level=0.95) # CI

#probit regression
gave_probit <- glm(out_gavedum ~ treatment,family=binomial(link='probit'))
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

#LA: added that we remove ppl missing income info
char_res <- char[ which(page18_39!=-999
                         & perbush!=-999
                         & median_hhincome!=-999), ] # drop all those with missings of key variables
detach(char)
attach(char_res) # attach so don't have to call each time

char_res$drop <- 0 # variable telling us to drop or not

# Make threshold rule for dropping (alternatively do with random variable)
char_res$thres <- perbush #+ 0.1*(1+perbush)^2 - 0.1*page18_39*perbush #- page18_39 - page18_39^2 + perbush^2
summary(char_res$thres)
char_res$drop[char_res$thres <= 0.5] <- 1
mean(char_res$drop) # drop 23 % of obs
char_res_d <- char_res[which(char_res$drop == 0),]

##############################
#Alternative rule: 
#randomly censor individuals
#via a  complex, highly nonlinear fcn  of votes 4 bush in state,
#

ps.fcn <- function(v,c,pg,t){
  v_t <- (v-.25)/.5
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
char_res$ps.true <- ps.fcn(char_res$perbush,char_res$cases,char_res$hpa,char_res$treatment)
ggplot(char_res,aes(x=ps.true))+ stat_ecdf()
set.seed(21) 
selection <- runif(nrow(char_res)) <= char_res$ps.true
char.censored <- char_res[selection,] #remove observations via propensity score rule
ggplot(char_res,aes(x=perbush)) + geom_histogram()+xlim(c(0,1))
ggplot(char.censored,aes(x=ps.true)) + geom_histogram() +xlim(c(0,1))

#overlap in true propensity score
ggplot(char.censored,aes(x=ps.true,colour=factor(treatment))) + stat_ecdf()
ggplot(char.censored,aes(x=ps.true,y=hpa,colour=factor(treatment))) + geom_point()
#there is clear overlap, but clearly assymetries going on with hpa as well
#End of Added PS Section by LA
#################################

# New regression results with dropping

#jack's threshold rule
reg_ols_drop <- lm(out_amountgive ~ treatment, data = char_res_d) 
summary(reg_ols_drop) 

#Luis' PS generating rule
reg_censored <- lm(out_amountgive ~ treatment, data = char.censored) 
summary(reg_censored) 

# Old regression results (remember to drop missings to make comparable sample)
reg_ols_comp <- lm(out_amountgive ~ treatment, data = char_res) 
summary(reg_ols_comp) 


# Check overlap (propensity score)
# estimate propensity score via logit regression
ps_mod <- glm(treatment ~ page18_39 + perbush,
              family = binomial(), data = char_res_d)
summary(ps_mod)

# Put propensity scores into dataframe with actual treatment

ps_df <- data.frame(pr_score = predict(ps_mod),
                     treatment = ps_mod$model$treatment)
head(ps_df)

# Generate graph (see https://stanford.edu/~ejdemyr/r-tutorials-archive/tutorial8.html#propensity-score-estimation)
labs <- paste("Actual treatment:", c("Treated", "Control"))
ps_df %>%
  mutate(treatment = ifelse(treatment == 1, labs[1], labs[2])) %>% 
  ggplot(aes(x = pr_score), data = ps_df) +
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
covars <- cbind(char.censored$hpa,char.censored$cases,char.censored$perbush)
char.censored$bias <- bias.fcn(char.censored$ps.true,char.censored$treatment,
                 char.censored$out_amountgive, covars)
ggplot(char.censored,aes(x=bias)) +geom_histogram(fill=I("white"),col=I("black"))
E.bias = mean(char.censored$bias)/mean(char.censored$treatment)*(1-mean(char.censored$treatment))
print(E.bias)

############################################################
### 3.


# propensity score weighting ATE


char.censored$w.ate[char.censored$treatment == 1] <-  1/char.censored$ps.true[char.censored$treatment == 1]
char.censored$w.ate[char.censored$treatment == 0] <-  ( 1 / (1 - char.censored$ps.true[char.censored$treatment == 0]))

  
pweight.reg <- lm(out_amountgive ~ treatment, weights = w.ate, data = char.censored)
summary(pweight.reg)

# direct regression analysis ATE;

# traditional double robust analysis weighting using inverse propensity score weighting; the lm command in R has a weights option.

# lasso or regularized logistic regression (optionally try CART or randomforest --classification trees, or classification forests), to estimate the propensity scoreand re-estimate the ATE using the methods above.


# a single-equation lasso of Y on X and W to estimate the ATE. Note that there is an option to not penalize the treatment indicator. You may
#wish to use that option anyway so that the treatment effect estimate is not
#shrunk. See http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html for
#the syntax for setting penalties for some coefficients to 0.


#Next try using the Belloni-Chernozhukov-Hansen method, where you use the lasso to select variables with non-zero coefficients from the propensity equation and the
#outcome equation, take the union of the two sets, and finally run OLS. You can either
#use cross-validation to choose lambda in each case, or you can follow the approaches
#suggested by BCH (those are more complicated but probably doesn't make a
#                  difference).


# Look at how your ATE coefficient changes with regularization
# Consider the single-equation LASSO of Y on X and W. Similar to the plot
# included in the handout, plot how the coefficient on the treatment indicator
# changes with lambda. Interpret your findings. (See http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html for some sample code on plotting.)

############################################################
### 4.

# double machine learning to estimate the ATE. Specifically, use a LASSO or random forest to estimate regressions of Y on X and separately Y on W. Then, run a residual on
# residual regression.


# Use approximate residual balancing (package: http://github.com/swager/balanceHD) to estimate ATE

# Compare and interpret your results