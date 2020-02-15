# Parameters
file_path<-"/Users/erinconrad/Desktop/residency stuff/R25/ige project/data/data_for_r.csv"
output_stats<-"/Users/erinconrad/Desktop/residency stuff/R25/ige project/data/r_stats.csv"

# Import packages
library(logistf)

# Load csv file
dr<-read.csv(file_path,header=TRUE,sep=",")


# Tell it which variables are categorical?????

# Do regular model
fit1<-logistf(data=dr,drug_resistant~pst+gpfa+duration_minutes,firth=FALSE,pl=FALSE)
summary(fit1)
exp(coef(fit1))

# Do penalized model
fit2<-logistf(data=dr,drug_resistant~pst+gpfa+duration_minutes,firth=TRUE,pl=TRUE)
summary(fit2)
exp(coef(fit2))


# Export statistics
