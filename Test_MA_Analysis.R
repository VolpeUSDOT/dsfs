
rm(list=ls()) # start with a blank slate! Remove everything from the workspace.

setwd("~/git/dsfs")

library(lme4)
library(sjPlot)
library(metafor)

# Load in test data
load("Fake_DSFS.RData")


# If we actually had raw values, we could do the following direct test of activation hypothesis:
m1 <- lmer(speeds ~ activation + activation:upstream + activation:downstream + (1|pub), data = dataraw)

summary(m1) # should be very similar to coefficients put in above.

# Random effects: the different groups
sjp.lmer(m1)

sjp.lmer(m1, type = 'fe', p.kr = FALSE)

sjt.lmer(m1, p.kr = FALSE)

# But instead we have summarized meta-analysis data. So here we can do several different tests:

# Simple activation: (Before - During)
# True effect: (Upstream - Adjacent)_before - (Upstream - Adjacent)_during
# Simple deactivation: (During - After)
# True effect, deactivation: (Upstream - Adjacent)_during - (Upstream - Adjacent)_after
# downstream "bounce back" test also possible


# Start with simple activation (H1B). Could try rma.glmm in metafor package
# Alternatively: log response ratio within study
