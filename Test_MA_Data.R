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

n_pub = 10
ntot = 500 # number of sample 

baseinter = 50
pubint <- baseinter + c(1:n_pub)-mean(1:n_pub) # different intercepts by species

test <- vector()

for(i in 1:n_pub){ # loop over publications. i = 1
  
# Build up the data frame. If we had factors (as in, an experimental design with definied treatment levels), would use gl to generate levels, and subtract 1 to make them start at 0 (default is 1) if there are multiple levels. For continuous predictors here, use rnorm. Easiest way to see if the input coefficients are working is to use mean 0, sd 1 in the predictors (as in, scaled predictors).
  
  activation = rnorm(ntot, 0, 1) # continuous predictor. Uncorrelated with each other
  deactivation = rnorm(ntot, 0, 1)
  upstream = rnorm(ntot, 0, 1)
  downstream = rnorm(ntot, 0, 1)
  
  ###### Set up differences for each level
  activecoef = -6 # strength of the activation effect in reducing speeds
  deactivecoef = 5 # strength of the deactivation effect in the resulting increase in speeds
  upstreamcoef = 6 # difference between upstream and adjacent measurements
  downstreamcoef = 4 # difference between downstream and adacent
  activeup = 0.5 # 
  activedown = 0.3
  deactiveup = -0.3 # 
  deactivedown = 0.1

  ######## SD for each treatment, can change these as necessary. 
  activecoef.sd = 1
  deactivecoef.sd = 0.5 
  upstreamcoef.sd = 0.1
  downstreamcoeff.sd = 1
  activeup.sd = 0.5
  activedown.sd = 0.5
  deactiveup.sd = 0.5
  deactivedown.sd = 0.5

  # Make model matrix. Limit to two-way interactions with + between predictors, and ( )^2 surrounding the term
  
  mm <- model.matrix(~(activation+deactivation+upstream+downstream+
                         activation:upstream + activation:downstream +
                         deactivation:upstream + deactivation:downstream
                         ), data.frame(activation, deactivation, upstream, downstream))


  # Coefficients need to match the order of the colums in the model matrix, mm. First one is the intercept
  coeff <- c(pubint[i], 
             rnorm(1, activecoef, activecoef.sd),
             rnorm(1, deactivecoef, deactivecoef.sd),
             rnorm(1, upstreamcoef, upstreamcoef.sd),
             rnorm(1, downstreamcoef, downstreamcoeff.sd), 
             rnorm(1, activeup, activeup.sd), 
             rnorm(1, deactiveup, deactiveup.sd), 
             rnorm(1, activedown, activedown.sd), 
             rnorm(1, deactivedown, deactivedown.sd)
  )
  
  speeds <- rnorm(n = ntot, mean = mm %*% coeff, sd = 0.1)
  
  testx <- data.frame(speeds, pub = i, 
                      activation, deactivation, upstream, downstream)
  
  test <- rbind(test, testx)  

}


library(lme4)
library(sjPlot)

m1 <- lmer(speeds ~ activation + activation:upstream + activation:downstream + (1|pub), data = test)

sjp.lmer(m1, type = "fe")
