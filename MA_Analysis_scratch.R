# Analysis script for meta-analysis of Dynamics Speed Feedback Signs
# Based on studies compiled by Volpe group 2015-2016
# Dan Flynn | daniel.flynn.ctr@dot.gov

setwd("~/git/dsfs")

library(lme4)
library(sjPlot)
library(scales) # for alpha()

# For MA_Report.Rmd -- first run MA_Analysis_Prep.R

source("MA_Analysis_Prep.R")

# Model with main effects. yi: vector of the effect sizes; vi: vector of sampling variances.

h1b.mod <- rma(yi, vi, mods = ~ location,
               data = h1b.md) # positive number is reduction in speed from before to during

summary(h1b.mod)

plot(coef(h1b.mod)[1:3], type="o", pch=19, 
#     xlim=c(.8,3.2), ylim=c(-.1,.5), 
     xlab="Decrease in speeds as a result of DSFS activation", 
     ylab="Standardized Mean Difference", xaxt="n", bty="l")
axis(side=1, at=1:3, labels=c("Upstream","Adjacent","Downstream"))


#forest(h1b.mod) 
#xlim=c(-7,5), alim=c(-3,1), cex=.8)


# Try analyzing effect sizes in standard mixed effect model analysis
# 
# summary(lm1 <- lmer(yi ~ location + (1|PublicationID), data = h1b.md))
# 
# sjp.lmer(lm1, type = "fe", show.intercept = TRUE)
# 
# sjt.lmer(lm1)

aggregate(h1b.md[,c("Mean.before","Mean.during","yi")],
                    by = list(h1b.md$PublicationID, h1b.md$location),
                    FUN = mean, na.rm=T)  
          
# nearly identical to rma() approach, and allows for multiple random effects and inclusion of additional continuous variables.

###############################

# summary(lm2 <- lmer(yi ~ location + (1 | vehicle.type) + (1 | PublicationID/StudyID), 
#                     data = h1b.md[h1b.md$sign.type != "DVMS",]))
# 
# sjp.lmer(lm2, type = "fe", show.intercept = T)
# 
# sjp.lmer(lm2, type = "re", show.intercept = T)
# 
# AIC(lm1, lm2) # improved model with study ID nested in publication, as it should be. 

# now add posted speed
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# ** Use this model ***

# H1B
summary(lm3 <- lmer(yi ~ location + posted.speed + 
                      (1|safety.focus2) +
                      (1|vehicle.type2) +
                      (1|PublicationID/StudyID/site), data = h1b.md))

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

sjp.lmer(lm3, type = "fe", show.intercept = TRUE)

# sjp.lmer(lm3, type = "re", show.intercept = TRUE)

sjt.lmer(lm3)

# # add sign type
# 
# summary(lm4 <- lmer(yi ~ location + posted.speed + (1|PublicationID/StudyID), data = h1b.md))
# 
# sjp.lmer(lm4, type = "fe", show.intercept = TRUE)

# H1A
summary(h1a.mod <- lmer(yi ~ posted.speed + 
                      (1|safety.focus2) +
                      (1|vehicle.type2) +
                      (1|PublicationID/StudyID/site), data = h1a.md))

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

sjp.lmer(lm3, type = "fe", show.intercept = TRUE)

# sjp.lmer(lm3, type = "re", show.intercept = TRUE)

sjt.lmer(lm3)
