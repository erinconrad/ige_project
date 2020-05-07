# This code is part of the IGE project (github.com/erinconrad/ige_project). It takes data
# generated from the Matlab function ige_stats.m and performs a survival analysis to test
# if drug resistant patients have earlier occurrence of PST than drug responsive patients.

# To run this, you will need to have installed the survival package in R
# https://cran.r-project.org/web/packages/survival/index.html

# Parameters
file_path<-"../data/cmh_table.csv" # the file path for the table of time to first PST
#file_path<-"../data/cmh_table.csv" # file path for time to first occurrence of various features

library(samplesizeCMH)


# Load csv file
data<-read.csv(file_path,header=FALSE,sep=",")

# Make 3 dimensional table
m <- data.matrix(data)
table3d <- array(m,dim=c(2,2,2),dimnames=NULL)

mantelhaen.test(table3d,correct=FALSE)

