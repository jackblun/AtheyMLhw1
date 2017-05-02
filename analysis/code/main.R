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

#library(stats)
#library(MASS)

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

### 1. Regression for average treatment effect

reg_ols <- lm(out_amountgive ~ treatment) 
summary(reg_ols) # show results, significant at 90% but not 95% level

# Consistent with Table 4 of paper

confint(reg_ols, level=0.95) # CI

### 

