
rm(list=ls()) # start with a blank slate! Remove everything from the workspace.

setwd("~/git/dsfs")

library(lme4)
library(sjPlot)
library(car)
library(metafor)
library(multcomp)

# metafor package is very useful for us
# http://www.metafor-project.org/doku.php/tips:multiple_factors_interactions
# rma() is the random mixed effects meta-analysis function in metafor. 

# Load in test data
load("Fake_DSFS.RData")

# make sure columns are in the right format

datax$pub <- as.factor(datax$pub)

# Rename and order the levels.
levels(datax$location) <- c("2 adjacent", "3 downstream", "1 upstream")
datax$location <- factor(as.character(datax$location))


# # If we actually had raw values, we could do the following direct test of activation hypothesis:
# m1 <- lmer(speeds ~ activation + activation:upstream + activation:downstream + (1|pub), data = dataraw)
 
# summary(m1) # should be very similar to coefficients put in initially  (see Test_MA_Data.R)
 
# # Random effects: the different groups
# sjp.lmer(m1)
 
# # Fixed effects: the treatments
# sjp.lmer(m1, type = 'fe', p.kr = FALSE)
 
# # Summary table, formatted nicely
# sjt.lmer(m1, p.kr = FALSE)

# But instead we have summarized meta-analysis data. So here we can do several different tests:

# Simple activation: (Before - During)
# True effect: (Upstream - Adjacent)_before - (Upstream - Adjacent)_during
# Simple deactivation: (During - After)
# True effect, deactivation: (Upstream - Adjacent)_during - (Upstream - Adjacent)_after
# downstream "bounce back" test also possible

# Start with simple activation (H1B). Could try rma.glmm in metafor package
# Alternatively: log response ratio within study
# use escalc() for calculating effect sizes.
# MD: raw mean difference
# SMD: standardized mean difference; Hedges' g
# SMDH: standardized mean difference with heteroscedastic population variance
# ROM: log transformed ratio of means

# Raw difference at same site
with(datax, tapply(during.mean-before.mean, location, mean))

h1b.smd <- escalc("SMD", 
              n1i = before.n,
              n2i = during.n,
              m1i = before.mean,
              m2i = during.mean,
              sd1i = before.sd,
              sd2i = during.sd,
              data = datax)


h1b.md <- escalc("MD", 
              n1i = before.n,
              n2i = during.n,
              m1i = before.mean,
              m2i = during.mean,
              sd1i = before.sd,
              sd2i = during.sd,
              data = datax)

h1b.rom <- escalc("ROM", 
                 n1i = before.n,
                 n2i = during.n,
                 m1i = before.mean,
                 m2i = during.mean,
                 sd1i = before.sd,
                 sd2i = during.sd,
                 data = datax)

# see these by location, 
tapply(h1b.smd$yi, h1b.smd$location, mean) # Use SMD for standardized mean difference.               

tapply(h1b.md$yi, h1b.md$location, mean) # Use MD for mean difference.               

# Model with main effects. yi: vector of the effect sizes; vi: vector of sampling variances.

h1b.mod <- rma(yi, vi, mods = ~ location,
               data = h1b.md, method = "HE")

summary(h1b.mod)

# Direct test specifically of difference between location adjacent to DSFS sign and upstream location
linearHypothesis(h1b.mod, c(1,0,-1,0))

# Direct test specifically of difference between location adjacent to DSFS sign and downstream location
linearHypothesis(h1b.mod, c(1,-1,0,0))

# Can also combine these in one test

summary(glht(
  h1b.mod,
  linfct = rbind(c(1,0,-1,0),
                 c(1,-1,-0,0))
        ))


plot(coef(h1b.mod)[1:3], type="o", pch=19, 
#     xlim=c(.8,3.2), ylim=c(-.1,.5), 
     xlab="Decrease in speeds as a result of DSFS activation", 
     ylab="Standardized Mean Difference", xaxt="n", bty="l")
axis(side=1, at=1:3, labels=c("Upstream","Adjacent","Downstream"))

forest(h1b.mod) 
       #xlim=c(-7,5), alim=c(-3,1), cex=.8)
text(-7, 11, "Study/Source",          pos=4, cex=.8)
text( 5, 11, "Observed SMD [95% CI]", pos=2, cex=.8)


# Try analyzing effect sizes in standard mixed effect model analysis


summary(lm1 <- lmer(yi ~ location + (1|pub), data = h1b.md))
sjp.lmer(lm1, type = "fe")
# nearly identical to rma() approach, and allows for multiple random effects and inclusion of additional continuous variables.


## Testing H3Aprime - 2018-03-05


  time1 = "after" 
  time2 = "during" 
  loc1 = "1 upstream" 
  loc2 = "2 adjacent" 
  data = dat
  measure = "MD"

  
  data$ID <- rownames(data)
  
  dx <- data[data$position == loc1 | data$position == loc2,
             c("ID",
               paste("Mean", time1, sep="."), paste("Mean", time2, sep="."),
               paste("SD", time1, sep="."), paste("SD", time2, sep="."),
               paste("N", time1, sep="."), paste("N", time2, sep=".") )
             ]
  keep <- apply(dx[2:length(dx)], 1, function(x) all(!is.na(x)))
  keeper <- names(keep)[keep]
  
  dat2 <- data[match(keeper, data$ID),]
  
  # also now make sure there are matches between both cases, reference and treatment. Everything that might vary, except for position
  dat2$compareID <- with(dat2, paste(PublicationID, StudyID, site, vehicle.type, sign.type, days.installed, distance.from.sign))
  
  lengths <- tapply(dat2$compareID, dat2$position, length)
  
  if(all(!is.na(lengths[loc1]), !is.na(lengths[loc2]))){ # Paired comparison check #1. proceed as long as there are some values at both locations. 
    # find which pairs to compare
    comparepairs <- unique(dat2$compareID)
    
    yix <- vix <- comparepairscomplete <- vector()
    
    for(i in comparepairs){ # i = comparepairs[3] # for kayemb, why doesn't it work for this publication?
      
      m1i.1 <- mean(dat2[dat2$compareID == i & dat2$position == loc1, paste("Mean", time1, sep=".")])
      m1i.2 <- mean(dat2[dat2$compareID == i & dat2$position == loc2, paste("Mean", time1, sep=".")])
      
      m2i.1 <- mean(dat2[dat2$compareID == i & dat2$position == loc1, paste("Mean", time2, sep=".")]) 
      m2i.2 <- mean(dat2[dat2$compareID == i & dat2$position == loc2, paste("Mean", time2, sep=".")])
      
      n1i.1 <- mean(dat2[dat2$compareID == i & dat2$position == loc1, paste("N", time1, sep=".")]) 
      n1i.2 <- mean(dat2[dat2$compareID == i & dat2$position == loc2, paste("N", time1, sep=".")])
      
      n2i.1 <- mean(dat2[dat2$compareID == i & dat2$position == loc1, paste("N", time2, sep=".")]) 
      n2i.2 <- mean(dat2[dat2$compareID == i & dat2$position == loc2, paste("N", time2, sep=".")])
      
      sd1i.1 <- mean(dat2[dat2$compareID == i & dat2$position == loc1, paste("SD", time1, sep=".")]) 
      sd1i.2 <- mean(dat2[dat2$compareID == i & dat2$position == loc2, paste("SD", time1, sep=".")])
      
      sd2i.1 <- mean(dat2[dat2$compareID == i & dat2$position == loc1, paste("SD", time2, sep=".")]) 
      sd2i.2 <- mean(dat2[dat2$compareID == i & dat2$position == loc2, paste("SD", time2, sep=".")])
      
      testvals <- c(m1i.1, m1i.2, m2i.1, m2i.2,  n1i.1, n1i.2, n2i.1, n2i.2, sd1i.1, sd1i.2, sd2i.1, sd2i.2)
      
      if(all(!is.na(testvals), length(testvals) == 12)){ # Paired comparison check #2. Make sure for this comparison that all values are available
        
        mi <- sum(n1i.1, n1i.2, n2i.1, n2i.2, na.rm=TRUE)-4
        
        sdpi <- sqrt( ((n1i.1 - 1) * sd1i.1^2 + (n1i.2 - 1) * sd1i.2^2 + (n2i.1 - 1) * sd2i.1^2 + (n2i.2 - 1) * sd2i.2^2 ) / mi )
        
        di <- ( (m1i.1-m1i.2) - (m2i.1 - m2i.2) )/sdpi
        
        if (measure == "MD") {
          yi <- (m1i.1-m1i.2) - (m2i.1 - m2i.2)
          
          #Use "LS" type sampling variances, large sample approximation
          vi <- sd1i.1^2/n1i.1 + sd1i.2^2/n1i.2 + sd2i.1^2/n2i.1 + sd2i.2^2/n2i.2
        }
        
        if (measure == "SMD") {
          # .cmicalc, hidden function in misc.func.hidden.r. Bias correction of SMDs
          cmi <- ifelse(mi <= 1, NA, exp(lgamma(mi/2) - log(sqrt(mi/2)) - lgamma((mi-1)/2)))
          
          yi <- cmi * di
          
          ni = sum(n1i.1, n1i.2, n2i.1, n2i.2)
          
          vi <- 1/n1i.1 + 1/n2i.1 + 1/n1i.1 + 1/n2i.2 + yi^2 / (4 * ni)
        }
        
        yix <- c(yix, yi)
        vix <- c(vix, vi)
        comparepairscomplete <- c(comparepairscomplete, i)
        
      } # end paired comparison check # 2
      
    } # end comparepairs loop
    
    result <- dat2[match(comparepairscomplete, dat2$compareID),c("PublicationID", "StudyID", "site", "vehicle.type2","safety.focus2", "sign.type", "days.installed", "distance.from.sign","posted.speed")]
    
  testout <-  data.frame(result, yi = yix, vi = vix)
  } # end paired com
