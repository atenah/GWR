library("GWmodel")

### load the functions
source("C:/Users/atenah/Desktop/Paper_2/GWR_R_Scripts/gwr.predict.mod.R")
source("C:/Users/atenah/Desktop/Paper_2/GWR_R_Scripts/gwr.predict.mod.quick.R")
source("C:/Users/atenah/Desktop/Paper_2/GWR_R_Scripts/gwr.bootstrap.R")

### import data and convert to proper format
input = read.delim("C:/Users/atenah/Desktop/Paper_2/GWR_R_Scripts/Test_time_series.csv", sep=",")
input = SpatialPointsDataFrame(coords = input[,1:2], data = input[,-c(1:2)])

###################################################################################
###################################################################################
### SET UP YOUR REGRESSION MODEL AND PARAMETERS HEREE!!! ###

### regression model
reg.model = "Yield~March2_a_1"

### output directory (NOTE!!! DATA WILL BE OVERWRITTEN!!!)
### use a new directory, just to be safe
output.folder = "C:/Users/atenah/Desktop/Paper_2/GWR_R_Scripts/Output_GWR"

### kernal type
kernal.type = "gaussian"

### specify training and test row numbers
training.rows = 1:800
test.rows = 801:900

### Do you want to perform bootstrap?
perform.bootstrap = T  ### answer T or F
num.bootstrap = 1000  ### e.g. run the bootstrap 1000 times
percent.test = 10     ### e.g. within each run, 10% will be test data and 90% will be training data

###################################################################################
###################################################################################
###################################################################################
###################################################################################


###################################################################################
###################################################################################
############# RUN EVERYTHING FROM HERE ONWARDS ####################################
###################################################################################
###################################################################################

dir.create(output.folder)

### calculate GW distance matrix
DM = gw.dist(dp.locat=coordinates(input))

### identify best bw
bw.best = bw.gwr(formula(reg.model), data=input, kernel = "gaussian", dMat = DM)

### typical regression using all
### gwr.res = gwr.basic(formula(reg.model), data=input, bw=bw.best, kernel = kernal.type, dMat=DM)
gwr.res = gwr.predict.mod(formula(reg.model), input[training.rows,], input[test.rows,], bw=bw.best,  kernel = kernal.type, generate.prediction.var=T)

### print gwr results to file
t.ssq = sum((gwr.res$actual-mean(gwr.res$actual))^2)
r.ssq = sum((gwr.res$actual-gwr.res$prediction)^2)
R2 = 1-r.ssq/t.ssq
write(paste("R2: ", R2, sep=""), paste(output.folder, "/GW_regression_results.txt", sep=""))

sink(paste(output.folder, "/GW_regression_results_details.txt", sep=""))
gwr.res
sink()


if (perform.bootstrap) {
  gwr.bootstrap(formula(reg.model), input, output.folder, num.bootstrap, percent.test, kernal.type)
}
###################################################################################
###################################################################################
###################################################################################
###################################################################################
############# YOU WILL FIND THE RESULTS IN OUTPUT FOLDER ##########################
###################################################################################
###################################################################################

