# Analysis script for meta-analysis of Dynamics Speed Feedback Signs
# Based on studies compiled by Volpe group 2015-2016
# Dan Flynn | daniel.flynn.ctr@dot.gov

setwd("~/git/dsfs")

library(lme4)
library(sjPlot)
library(scales) # for alpha()

# For MA_Report.Rmd -- first run MA_Analysis_Prep.R

source("MA_Analysis_Prep.R") # or load the resulting .RData

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# safety.focus2 and vehicle.type2 are cleaned versions of these variables 
# H1B
summary(h1b.mod <- lmer(yi ~ location + posted.speed + 
                      (1|safety.focus2) +
                      (1|vehicle.type2) +
                      (1|PublicationID/StudyID/site), data = h1b.md))

summary(h1b.mod2 <- lmer(yi ~ location + posted.speed + 
                          (location|safety.focus2) +
                          (location|vehicle.type2) +
                          (1|PublicationID/StudyID/site), data = h1b.md))


sjp.lmer(h1b.mod, type = "fe", show.intercept = TRUE, p.kr = F)

sjp.lmer(h1b.mod2, type = "re",p.kr = F)

sjp.lmer(h1b.mod2, type = "eff.ri", p.kr = F)

sjp.lmer(h1b.mod2, type = "rs.ri")

sjp.lmer(h1b.mod2, type = "ri.slope")


sjt.lmer(h1b.mod2,p.kr = F)



summary(h1b.smod <- lmer(yi ~ location + posted.speed + 
                          (1|safety.focus2) +
                          (1|vehicle.type2) +
                          (1|PublicationID/StudyID/site), data = h1b.smd))

sjp.lmer(h1b.smod, type = "fe", show.intercept = TRUE)

sjt.lmer(h1b.smod)


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# H1A. Remove 'location', normalized already
summary(h1a.mod <- lmer(yi ~ posted.speed + 
                      (1|safety.focus2) +
                      (1|vehicle.type2) +
                      (1|PublicationID/StudyID/site), data = h1a.md))

sjp.lmer(h1a.mod, type = "fe", show.intercept = TRUE)

sjt.lmer(h1a.mod)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# H1C. Remove 'location', normalized already
summary(h1c.mod <- lmer(yi ~ posted.speed + 
                          (1|safety.focus2) +
                          (1|vehicle.type2) +
                          (1|PublicationID/StudyID/site), data = h1c.md))

sjp.lmer(h1c.mod, type = "fe", show.intercept = TRUE)

sjt.lmer(h1c.mod)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# H2A. Remove 'location', normalized already
summary(h2a.mod <- lmer(yi ~ posted.speed + 
                          (1|safety.focus2) +
                          (1|vehicle.type2) +
                          (1|PublicationID/StudyID/site), data = h2a.md))

sjp.lmer(h2a.mod, type = "fe", show.intercept = TRUE)

sjt.lmer(h2a.mod)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# H2B: see H1B, downstream


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# H3B. Only at adjacent, don't need location 
summary(h3b.mod <- lmer(yi ~ location + posted.speed + 
                          (location|safety.focus2) +
                          (location|vehicle.type2) +
                          (1|PublicationID/StudyID/site),
                        data = h3b.md))

sjp.lmer(h3b.mod, type = "fe", show.intercept = TRUE)
sjp.lmer(h3b.mod, type = "re")

sjt.lmer(h3b.mod)
# halo effect: 1.26 mph compared to upstream

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# H3B'. Only at adjacent, don't need location 
summary(h3bprime.mod <- lmer(yi ~ location + posted.speed + 
                          (location|safety.focus2) +
                          (location|vehicle.type2) +
                          (1|PublicationID/StudyID/site),
                        data = h3bprime.md))

sjp.lmer(h3bprime.mod, type = "fe", show.intercept = TRUE)

sjp.lmer(h3bprime.mod, type = "re")


sjt.lmer(h3b.mod)
