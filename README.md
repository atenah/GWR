# GWR
Geographically Weighted Regression.
The input data: Firt column should be Longitude, The second column : Lattitude.
Run the "GWR cross validation" script, you SET UP YOUR REGRESSION MODEL AND PARAMETERS HERE in this scrip. change the input data accordingly.
Kernel type can be set to either: "bisquare","gaussian".
an example:
kernal.type = "bisquare",
adaptive.type = T,
longlat.type = F,
approach.type = "AICc",
p = 2,
theta = 0
