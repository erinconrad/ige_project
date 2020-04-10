# This code is part of the IGE project (github.com/erinconrad/ige_project). It takes data
# generated from the Matlab function ige_stats.m and performs a survival analysis to test
# if drug resistant patients have earlier occurrence of PST than drug responsive patients.

# To run this, you will need to have installed the survival package in R
# https://cran.r-project.org/web/packages/survival/index.html

# Parameters
file_path<-"../data/r_table_pst.csv" # the file path for the table of time to first PST

# Import packages
library(survival)

# Load csv file
data<-read.csv(file_path,header=TRUE,sep=",")

# perform the log-rank test to test for a difference in "survival time" (time to first PST) between drug resistant and responsive patients 
survdiff(Surv(survtime,observed)~resistant,data=data,rho=0)

# Kaplan-Meier plot (should give the same plot as in Matlab)
plot(survfit(Surv(survtime,observed)~resistant,data=data),mark.time =TRUE)
