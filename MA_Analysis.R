

# Analysis script for meta-analysis of Dynamics Speed Feedback Signs
# Based on studies compiled by Volpe group 2015-2016
# Dan Flynn | daniel.flynn.ctr@dot.gov

setwd("~/git/dsfs")

# For MA_Report.Rmd -- first run MA_Analysis_Prep.R

source("MA_Analysis_Prep.R")

# Model with main effects. yi: vector of the effect sizes; vi: vector of sampling variances.

h1b.mod <- rma(yi, vi, mods = ~ location,
               data = h1b.md)

summary(h1b.mod)

plot(coef(h1b.mod)[1:3], type="o", pch=19, 
#     xlim=c(.8,3.2), ylim=c(-.1,.5), 
     xlab="Decrease in speeds as a result of DSFS activation", 
     ylab="Standardized Mean Difference", xaxt="n", bty="l")
axis(side=1, at=1:3, labels=c("Upstream","Adjacent","Downstream"))

forest(h1b.mod) 
       #xlim=c(-7,5), alim=c(-3,1), cex=.8)


# Try analyzing effect sizes in standard mixed effect model analysis

summary(lm1 <- lmer(yi ~ location + (1|PublicationID), data = h1b.md))
sjp.lmer(lm1, type = "fe", show.intercept = TRUE)



aggregate(h1b.md[,c("Mean.before","Mean.during","yi")],
                    by = list(h1b.md$PublicationID, h1b.md$location),
                    FUN = mean, na.rm=T)  
          
# nearly identical to rma() approach, and allows for multiple random effects and inclusion of additional continuous variables.

###############################


summary(lm2 <- lmer(yi ~ location + (1 | vehicle.type) + (1 | PublicationID/StudyID), 
                    data = h1b.md[h1b.md$sign.type != "DVMS",]))


sjp.lmer(lm2, type = "fe", show.intercept = T)

sjp.lmer(lm2, type = "re", show.intercept = T)



AIC(lm1, lm2) # improved model with study ID nested in publication, as it should be. With 

# now add posted speed

summary(lm3 <- lmer(yi ~ location + posted.speed + 
                      (1|safety.focus) +
                      (1|vehicle.type) +
                      (1|PublicationID/StudyID), data = h1b.md))

sjp.lmer(lm3, type = "fe", show.intercept = TRUE)

sjp.lmer(lm3, type = "re", show.intercept = TRUE)


# add sign type

summary(lm4 <- lmer(yi ~ location + posted.speed + (1|PublicationID/StudyID), data = h1b.md))

sjp.lmer(lm3, type = "fe", show.intercept = TRUE)



# Test of Simple Deactivation hypothesis (H3B)

h3b.md <- escalc("MD", 
                 n1i = N.during,
                 n2i = N.after,
                 m1i = Mean.during,
                 m2i = Mean.after,
                 sd1i = SD.during,
                 sd2i = SD.after,
                 data = dat)

# how many values can we actually calculate deactivation effects?
summary(is.na(h3b.md$yi)) # only 42 total values here.

tapply(h3b.md$yi, h3b.md$location, mean, na.rm=T) # Use MD for mean difference.               

summary(lm2.deactive <- lmer(yi ~ location + (1|PublicationID/StudyID), data = h3b.md))

# compare to the right data set

h3b.1b.md <- escalc("MD", 
                    n1i = N.before,
                    n2i = N.during,
                    m1i = Mean.before,
                    m2i = Mean.during,
                    sd1i = SD.before,
                    sd2i = SD.during,
                    data = dat[!is.na(h3b.md$yi),])
summary(lm2.active <- lmer(yi ~ location + (1|PublicationID/StudyID), data = h3b.1b.md))



sjp.lmer(lm2.deactive, type = "fe", show.intercept = TRUE)


hist(dat$Mean.before, xlim = c(20, 70))

hist(dat$Mean.during, add = T, col = alpha("blue", 0.2))

hist(dat$Mean.after, add = T, col = alpha("red", 0.2))

aggregate(dat[,c("Mean.before","Mean.during","Mean.after")],
          by = list(dat$location),
          FUN = mean, na.rm=T)


