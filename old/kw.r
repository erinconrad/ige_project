# Parameters
file_path<-"../data/r_table_kw.csv"

# Import packages

# Load csv file
data<-read.csv(file_path,header=FALSE,sep=",")


# Survdiff
kruskal.test(V1 ~ V2,data = data)