# Computing codes for the simulation studies and application in the manuscript entitled "Periodicity Learning in Non-stationary Functional Time Series"

Codes for simulation are in the folder named "simulation".

simu_mtfts.R: main procedure for calculating Table~1 in simulation 1
simu_mtfts2.R: main procedure for calculating Table~2 in simulation 2
my.bibasis1.R: calculate the estimators and residual sum of squares
gcvbic.R: find optimal number of knots

Codes and data for real data applications are in the folder named "real data".

The CO_NOX file contains the data of CO and NOX  named hkair.csv and the codes.
The global temperature file contains the data of global temperature  named global1850.csv and the codes.
The sunspot file contains the data of sunspot  named sunspots.csv and the codes.

hkco_s.R: real data of CO 
hknox_s.R: real data of NOX 

globaltem.R: real data of global temperature

sunspots.R: real data of sunspots numbers 

Two common Functions:
my.bibasis1.R: calculate the estimators and residual sum of squares
gcvbic.R: find optimal number of knots
