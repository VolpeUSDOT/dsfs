# Generating test data for meta analysis
# To generate: n, mean, and sd for each study, with three time points and three spatial locations.
# Time: before, during, and after
# Location: upstream, adjacent, downstream

# Coefficients to estimate:
# - sign activation effect: difference between before and during
# - sign deactivation effect: difference between during and after
# - upstream location effect: reduction from upstream to adjacent
# - downstream location effect: increase from adjacent to downstream

# Lots of covariates to include in final model, including distance from sign, posted speed, vehicle type, type of sign... not included now. 


#### Set up. Set working directory, read in needed packages

rm(list=ls()) # start with a blank slate! Remove everything from the workspace.

setwd("~/git/dsfs")

library(lme4)
library(sjPlot)

# Set up number of publications and number of samples per publication.

n_pub = 100
nactive = 2
ndeactive = 2
nupstream = 2
ndownstream = 2

reps = ceiling(rlnorm(n_pub, mean = 6, sd = 1)) # number of samples per treatment group per publication. Let it vary by publication, on lognormal distribution.

ntot = nactive*ndeactive*nupstream*ndownstream*reps

###### Set up differences for each level, overall. These will have some publication-specific variation, in addition to an overall difference in speed intercept by publication.
activecoef = -6 # strength of the activation effect in reducing speeds
deactivecoef = -2 # strength of the deactivation effect in the resulting increase in speeds. By itself it is zero, only comes into play in the interaction with activation.
upstreamcoef = 2 # difference between upstream and adjacent measurements, absent any treatment
downstreamcoef = 1 # difference between downstream and adjacent, absent any treament

activedeactivecoef = 5 # interaction between activation and deactivation; this is the amount that speeds bounce back after deactivation.

activeup = 5 # interaction bewteen activation and location upstream. How much higher is speed upstream during activation?
activedown = 1 # interaction between activation and location downstream
deactiveup = -2 # interaction between deactivation and location upstream
deactivedown = -2 # interaction between deactivation and location downstream

######## SD for each treatment, can change these as necessary. 
activecoef.sd = 1
deactivecoef.sd = 0.01 
upstreamcoef.sd = 0.1
downstreamcoeff.sd = 1
activedeactivecoef.sd = 1
activeup.sd = 0.5
activedown.sd = 0.5
deactiveup.sd = 0.5
deactivedown.sd = 0.5

overallspeed = 50 # general intercept; the overall mean speed across all studies
overallspeed.sd = 10

speedsd.mean = 6 # for pulling out values of sd for each study
speedsd.sd = 1

dataraw <- datax <- vector()

# small function for summarizing
summariz <- function(x){
  n <- length(x[!is.na(x)])
  mean <- mean(x, na.rm=T)
  sd <- sd(x, na.rm=T)
  return(c(n, mean, sd))
}


# now loop over publication to generate coefficients. Since length varies by publication, include the generating of levels in the loop

for(i in 1:n_pub){ # i = 1

  activation = gl(nactive, reps[i], length = ntot[i], labels = (1:nactive)-1)
  deactivation = gl(ndeactive, reps[i]*nactive, length = ntot[i], labels = (1:ndeactive)-1)
  upstream = gl(nupstream, reps[i]*nactive*ndeactive, length = ntot[i], labels = (1:nupstream)-1)
  downstream = gl(ndownstream, reps[i]*nactive*ndeactive*nupstream, length = ntot[i], labels = (1:ndownstream)-1)
  
  # Make model matrix. Limit to two-way interactions with + between predictors, and ( )^2 surrounding the term
  
  mm <- model.matrix(~(activation+deactivation+upstream+downstream+
                         activation:deactivation +
                         activation:upstream + activation:downstream +
                         deactivation:upstream + deactivation:downstream
  ), data.frame(activation, deactivation, upstream, downstream))
  
  
  # Coefficients need to match the order of the colums in the model matrix, mm. First one is the intercept
  # each publication has its own effect size for activation, deactivation, etc. 
  coeff <- c(rnorm(1, overallspeed, overallspeed.sd), 
             rnorm(1, activecoef, activecoef.sd),
             rnorm(1, deactivecoef, deactivecoef.sd),
             rnorm(1, upstreamcoef, upstreamcoef.sd),
             rnorm(1, downstreamcoef, downstreamcoeff.sd), 
             rnorm(1, activedeactivecoef, activedeactivecoef.sd), 
             rnorm(1, activeup, activeup.sd), 
             rnorm(1, deactiveup, deactiveup.sd), 
             rnorm(1, activedown, activedown.sd), 
             rnorm(1, deactivedown, deactivedown.sd)
  )

  # Now generate actual speeds for each publication.  
  speeds <- rnorm(n = ntot[i], mean = mm %*% coeff, sd = rnorm(1, speedsd.mean, speedsd.sd))
  
  testx <- data.frame(speeds, pub = i, 
                      activation, deactivation, upstream, downstream)
  
  testx$treat = paste(activation, deactivation, upstream, downstream)

  testx$treat <- as.factor(testx$treat)
  levels(testx$treat) = c("before.adj",
                          "before.down",
                          "before.up",
                          "null",
                          "null2",
                          "null3",
                          "null4",
                          "null5",
                          "during.adj", # 9
                          "during.down",
                          "during.up",
                          "null6",
                          "after.adj",
                          "after.down",
                          "after.up",
                          "null7")

  # ditch unused levels (e.g., deactivation without activation)
  testx <- testx[grep("null", testx$treat, invert = T),]
  testx$treat <- as.factor(as.character(testx$treat))
  
  # and summarize to the organization of data actually colllected: n, mean, and sd, before, during and after, upstream, adjacent, and downstream.
  
  xx <- tapply(testx$speeds, testx$treat, summariz)
  n <- unlist(lapply(xx, function(x) x[1]))
  m <- unlist(lapply(xx, function(x) x[2]))
  s <- unlist(lapply(xx, function(x) x[3]))
  
  xx <- data.frame(n, mean=m, sd=s)
  
  # organize as for data collection
  
  xx <- cbind(xx[4:6,], xx[7:9,], xx[1:3,])
  
  dimnames(xx) <- list(c("adjacent", "downstream", "upstream"),
                       paste(rep(c("before","during","after"), each = 3), c("n", "mean", "sd"), sep = "."))
  
  
  dataraw <- rbind(dataraw, testx)  
  
  datax <- rbind(datax, data.frame(pub=i, location=rownames(xx), xx))
  
  }
rownames(datax) = 1:nrow(datax)

head(datax)
tail(datax)

save(list=c("datax","dataraw"), file = "Fake_DSFS.RData")

# should all be similar before.
# effect should only appear during (adjacent less than upstream.
# effect should disappear somewhat after
