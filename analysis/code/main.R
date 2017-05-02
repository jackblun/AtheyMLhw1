############################################################
### ECON 293 Homework 1

# Luis Armona and Jack Blundell
# Stanford University
# Spring 2017

############################################################

# set your working directory

setwd("C:\\Users\\Jack\\Documents\\Git\\Athey ML homework 1\\AtheyMLhw1")

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
