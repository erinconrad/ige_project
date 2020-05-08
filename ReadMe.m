%{
This codebase contains the 3 scripts (2 Matlab, one R) and de-identified
data needed to replicate the analysis in the paper "Using generalized 
polyspike train to predict drug-resistant idiopathic generalized epilepsy".

The codebase and data are available at:
https://github.com/erinconrad/ige_project

Matlab scripts were run on Matlab R2019a. The R script requires the
Survival package to run, which can be downloaded at:
https://cran.r-project.org/web/packages/survival/index.html


To run the code, download the github code base and navigate to the folder
containing ige_stats.m. Edit the file paths in the "Parameters" section at
the top of the ige_stats function to indicate where the de-identified data
csv file is and where you would like to output results. Also, make sure the
variable "doing_from_github" is set to 1 so that it does not attempt to
call inaccessible functions.

Then run the code:

>> ige_stats

This will run the code and save all tables and figures to the results
folder you indicate. It will also save a table of data to be used as the
input for survival.r, which is the script to perform the survival analysis
to test for a difference in the time to first PST between drug-resistant
and drug-responsive patients.


Erin Conrad, 2020, University of Pennsylvania


%}