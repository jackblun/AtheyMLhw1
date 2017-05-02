############################################################
### ECON 293 Homework 1

# Luis Armona and Jack Blundell
# Stanford University
# Spring 2017


# Section numbers correspond to assignment page

############################################################

# set your working directory

setwd("C:\\Users\\Jack\\Documents\\Git\\Athey ML homework 1\\AtheyMLhw1") # Jack

# clear things in RStudio

rm(list = ls())


# Call packages

library(ggplot2)
library(dplyr)

# set seed
set.seed(12345)

############################################################

# Load data

char <- read.csv("analysis\\input\\charitable_withdummyvariables.csv")
attach(char) # attach so don't have to call each time


### Exploratory analysis

dim(char) # 50,083 obs, 63 vars
names(char)
head(char) # Look at first few entries of each var

# Treatment

summary(treatment) # This var seems to pool treatments
mean(treatment) # 67% treated


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

############################################################
### 2. Dropping some observations

# Try to find some drop of observations based on observables that changes the treatment effect

# look at some potential variables

summary(page18_39)
hist(page18_39[page18_39_missing!=1])

summary(perbush)
hist(perbush[perbush_missing!=1])


char_res <- char[ which(page18_39_missing!=1 
                         & perbush_missing!=1), ] # drop all those with missings of key variables
#attach(char_res) # attach so don't have to call each time

char_res$drop <- 0 # variable telling us to drop or not

# Make threshold rule for dropping (alternatively do with random variable)
char_res$thres <- perbush #+ 0.1*(1+perbush)^2 - 0.1*page18_39*perbush #- page18_39 - page18_39^2 + perbush^2

summary(char_res$thres)

char_res$drop[char_res$thres <= 0.5] <- 1

mean(char_res$drop) # drop 23 % of obs

char_res_d <- char_res[which(char_res$drop == 0),]


# New regression results with dropping
reg_ols_drop <- lm(out_amountgive ~ treatment, data = char_res_d) 
summary(reg_ols_drop) 

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

############################################################
### 3.
