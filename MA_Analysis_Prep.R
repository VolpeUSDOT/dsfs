
# Prep script for meta-analysis of Dynamics Speed Feedback Signs
# Based on studies compiled by Volpe group 2015-2016
# Dan Flynn | daniel.flynn.ctr@dot.gov

# Set working directory and read in necessary scripts
setwd("~/git/dsfs")

library(lme4)
library(sjPlot)
library(car)
library(metafor)
library(multcomp)
library(scales)
library(xtable)

# metafor package is very useful for us
# http://www.metafor-project.org/doku.php/tips:multiple_factors_interactions
# rma() is the random mixed effects meta-analysis function in metafor. 

# Load in data. Map network drive to your computer if necessary
dat <- read.csv("Y:/Meta-Analysis/Meta-Analyses/DSFS Meta-Analysis Data.csv")

# make sure columns are in the right format

dat$pub <- as.factor(dat$PublicationID)

# Rename and order the levels.
levels(dat$position) <- c("2 adjacent", "3 downstream", "1 upstream")
dat$location <- factor(as.character(dat$position))

# Simple activation: (Before - During). Combines at DSFS (H1B), downstream (H2B)

# True effect: (Upstream - Adjacent)_before - (Upstream - Adjacent)_during. This is H1A
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
# with(dat, tapply(Mean.during-Mean.before, location, mean, na.rm=T))

h1b.smd <- escalc("SMD", 
              n1i = N.before,
              n2i = N.during,
              m1i = Mean.before,
              m2i = Mean.during,
              sd1i = SD.before,
              sd2i = SD.during,
              data = dat)


h1b.md <- escalc("MD", 
                 n1i = N.before,
                 n2i = N.during,
                 m1i = Mean.before,
                 m2i = Mean.during,
                 sd1i = SD.before,
                 sd2i = SD.during,
                 data = dat)

h1b.rom <- escalc("ROM", 
                  n1i = N.before,
                  n2i = N.during,
                  m1i = Mean.before,
                  m2i = Mean.during,
                  sd1i = SD.before,
                  sd2i = SD.during,
                  data = dat)

# see these by location, 
# tapply(h1b.smd$yi, h1b.smd$location, mean, na.rm=T) # Use SMD for standardized mean difference.               
# tapply(h1b.md$yi, h1b.md$location, mean, na.rm=T) # Use MD for mean difference.               
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
# summary(is.na(h3b.md$yi)) # only 42 total values here.

# tapply(h3b.md$yi, h3b.md$location, mean, na.rm=T) # Use MD for mean difference.        
aggregate(dat[,c("Mean.before","Mean.during","Mean.after")],
          by = list(dat$location),
          FUN = mean, na.rm=T)


