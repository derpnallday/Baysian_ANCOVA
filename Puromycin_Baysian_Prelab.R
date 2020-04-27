#Puromycin Bayeisan Pre-lab

#below is code from the Fischerian lab

# set up
require (car) # for normal quantile plots (SVAs)

#load Puromycin dataset
?Puromycin
library(datasets)
data(Puromycin)

#get summary of dataset
summary(Puromycin)

#look at variable distributions for each state
lattice::densityplot(Puromycin$conc, 
                     groups=Puromycin$state,
                     auto.key = TRUE,
                     xlab="concentration")

#slight left skew for both -> suggests a log transform

lattice::densityplot(log(Puromycin$conc), 
                     groups=Puromycin$state,
                     auto.key = TRUE,
                     xlab="log(concentration)")
#now looks like normal distribution

lattice::densityplot(Puromycin$rate, 
                     groups=Puromycin$state,
                     auto.key = TRUE,
                     xlab="rate")
#looks like a normal distribution but treated is skewed to
#right compared to untreated

#See how the data fits with regression
lattice::xyplot(Puromycin$rate~log(Puromycin$conc), 
                type = c("p","r"),
                auto.key = TRUE)

#does it look better with groups?
lattice::xyplot(Puromycin$rate~log(Puromycin$conc), 
                groups = Puromycin$state,
                type = c("p","r"),
                auto.key = TRUE)




###########################################################
#
# Model Testing
#
###########################################################


#predict the reaction rate based on ppm concentration with/without state

# simple linear regression
m1 = lm(Puromycin$rate~Puromycin$conc)
summary(m1)
confint(m1)

#check SVAs
par(mfrow=c(1,2))
qqPlot(resid(m1))
plot(resid(m1)~Puromycin$conc)
abline(h=0)

###################
# ancova model
m2 = lm( Puromycin$rate~Puromycin$conc * Puromycin$state)
summary(m2)
confint(m2)

#check SVAs
par(mfrow=c(1,2))
qqPlot(resid(m2))
plot(resid(m2)~Puromycin$conc)
abline(h=0)

####################
# simple linear regression with log transform
m3 = lm(Puromycin$rate~log(Puromycin$conc))
summary(m3)
confint(m3)

#check SVAs 
par(mfrow=c(1,2))
qqPlot(resid(m3))
plot(resid(m3)~log(Puromycin$conc))
abline(h=0)

####################
# ancova model with log transform
m4 = lm(Puromycin$rate~log(Puromycin$conc) * Puromycin$state)
summary(m4)
confint(m4)

#check SVAs 
par(mfrow=c(1,2))
qqPlot(resid(m4))
plot(resid(m4)~log(Puromycin$conc))
abline(h=0)

####################

#test models against eachother
AIC(m1,m2,m3,m4)
#non log models
anova(m1,m2)
#log transform models
anova(m3,m4)

###################################################
#Conclusion

#There is a clea relationship between the reaction rate and the log transform of concentration. Also
#We can see that the addition of the state group variable further makes the model more accurate. This means
#We may want to use a Bayesian ANCOVA to include the groups and we might see that including individual sigma parameters for
#the groups will improve our understanding of the data

# Load The data file 

myDataFrame = Puromycin
attach(myDataFrame)
myDataFrame = transform(Puromycin, logConc = log10(conc))



# Specify the column names in the data file relevant to the analysis:
yName="rate" 
xNomName="state" 
xMetName="logConc"             # the covariate

# Specify desired contrasts of slopes.
# Each main-effect contrast is a list of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL):
contrasts = list( 
  list( c("treated") , c("untreated") , compVal=0.0 , ROPE=c(-0,5,0.5) ) 
)
# Specify filename root and graphical format for saving output.
# Otherwise specify as NULL or leave saveName and saveType arguments 
# out of function calls.
fileNameRoot = "Puro-ANCOVA-" 
graphFileType = "eps" 

#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("ANCOVA.R")
#source("Jags-Ymet-Xnom1met1-MnormalHom.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , 
                    yName=yName , xNomName=xNomName , xMetName=xMetName ,
                    numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , datFrm=myDataFrame , xNomName=xNomName , 
                        xMetName=xMetName , contrasts=contrasts , 
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , datFrm=myDataFrame , yName=yName , xNomName=xNomName , 
          xMetName=xMetName , contrasts=contrasts , 
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
# treated group equation = a[1] + aMet[1] * log(conc)
# untreated group equationm = a[2] + aMeta[2] * log(conc)
# a sigma = std for both a;s
# y sigma = std for t-dist

