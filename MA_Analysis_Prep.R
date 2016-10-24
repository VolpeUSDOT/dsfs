
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
dat <- read.csv("Y:/Meta-Analysis/Meta-Analyses/DSFS Meta-Analysis Data For R.csv")

# delete any studies which are incomplete
# dat$PublicationID

incompletes <- dat$Mean.during=="" & dat$position==""

dat <- dat[!is.na(!incompletes),]

# make sure columns are in the right format

dat$pub <- as.factor(dat$PublicationID)

# Rename and order the levels.
dat$position <- factor(as.character(dat$position))
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

# use time1: before, during, after
# loc1: upstream, adjacent, downstream


# in the case of calculating hypothesis H1A, need to have a separate escalc ("Effect Size Calculation") function to normalize between two time points. E.g. bertini06
TEST = TRUE#FALSE
if(TEST){
  data = dat[dat$PublicationID == "gambatese15",]
  time1='during'
  time2='before'
  loc1='2 adjacent'
  loc2='1 upstream'
  }

escalc.normalize.time <- function(time1, time2, loc1, loc2, data, measure = c("SMD", "MD")){
  # pull data using appropriate time and location
  # first need to line up data sets for which all data are available
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
  
  for(i in comparepairs){ # i = comparepairs[1]
  
  m1i.1 <- dat2[dat2$compareID == i & dat2$position == loc1, paste("Mean", time1, sep=".")] 
  m1i.2 <- dat2[dat2$compareID == i & dat2$position == loc1, paste("Mean", time2, sep=".")]
  
  m2i.1 <- dat2[dat2$compareID == i & dat2$position == loc2, paste("Mean", time1, sep=".")] 
  m2i.2 <- dat2[dat2$compareID == i & dat2$position == loc2, paste("Mean", time2, sep=".")]
  
  n1i.1 <- dat2[dat2$compareID == i & dat2$position == loc1, paste("N", time1, sep=".")] 
  n1i.2 <- dat2[dat2$compareID == i & dat2$position == loc1, paste("N", time2, sep=".")]

  n2i.1 <- dat2[dat2$compareID == i & dat2$position == loc2, paste("N", time1, sep=".")] 
  n2i.2 <- dat2[dat2$compareID == i & dat2$position == loc2, paste("N", time2, sep=".")]

  sd1i.1 <- dat2[dat2$compareID == i & dat2$position == loc1, paste("SD", time1, sep=".")] 
  sd1i.2 <- dat2[dat2$compareID == i & dat2$position == loc1, paste("SD", time2, sep=".")]
  
  sd2i.1 <- dat2[dat2$compareID == i & dat2$position == loc2, paste("SD", time1, sep=".")] 
  sd2i.2 <- dat2[dat2$compareID == i & dat2$position == loc2, paste("SD", time2, sep=".")]
  
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

  result <- dat2[match(comparepairscomplete, dat2$compareID),c("PublicationID", "StudyID", "site", "vehicle.type", "sign.type", "days.installed", "distance.from.sign")]
  
  data.frame(result, yi = yix, vi = vix)
  } # end paired comparison check # 1
  
}

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# 1. Activation hypothesis
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

h1a.md <- escalc.normalize.time(
            time1 = "during", 
            time2 = "before", 
            loc1 = "2 adjacent", 
            loc2 = "1 upstream", 
            data = dat, 
            measure = "MD")

h1a.smd <- escalc.normalize.time(
  time1 = "during", 
  time2 = "before", 
  loc1 = "2 adjacent", 
  loc2 = "1 upstream", 
  data = dat, 
  measure = "SMD")


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

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# 2. Downstream
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

h1a.md <- escalc.normalize.time(
  time1 = "during", 
  time2 = "before", 
  loc1 = "2 adjacent", 
  loc2 = "1 upstream", 
  data = dat, 
  measure = "MD")

h1a.smd <- escalc.normalize.time(
  time1 = "during", 
  time2 = "before", 
  loc1 = "2 adjacent", 
  loc2 = "1 upstream", 
  data = dat, 
  measure = "SMD")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# 3. Deactivation
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

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

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Summary statistics
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Raw difference at same site
with(dat, tapply(Mean.during-Mean.before, location, mean, na.rm=T))

# see these by location, 
tapply(h1b.smd$yi, h1b.smd$location, mean, na.rm=T) # Use SMD for standardized mean difference.               

tapply(h1b.md$yi, h1b.md$location, mean, na.rm=T) # Use MD for mean difference. Will be identical to raw diff calculation above
tapply(h3b.md$yi, h3b.md$location, mean, na.rm=T) # Using MD for mean difference.        

aggregate(dat[,c("Mean.before","Mean.during","Mean.after")],
          by = list(dat$location),
          FUN = mean, na.rm=T)

# Summary about dataset in general: Which hypotheses can we address?

id <- paste(dat$PublicationID, dat$StudyID, dat$ID)

aggregate(dat[,c("Mean.before", "Mean.during", "Mean.after")],
          by = list(id),
          FUN = function(x) length(x[!is.na(x)]))

means.loc <- 
  aggregate(dat[,c("Mean.before", "Mean.during", "Mean.after")],
            by = list(id, dat$location),
            FUN = function(x) length(x[!is.na(x)])>0)

sd.loc <- 
  aggregate(dat[,c("SD.before", "SD.during", "SD.after")],
            by = list(id, dat$location),
            FUN = function(x) length(x[!is.na(x)])>0)

# H1A: need mean and sd (and N, assume if have sd also have N). During and before, upstream and adjacent.
down <- sum(apply(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.before", "Mean.during")],
          1, all))
adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.before", "Mean.during")],
                1, all))

h1a <- min(down, adj)

# H1B
adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.before", "Mean.during")],
                 1, all))

h1b <- adj

# H1C 
down <- sum(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.during")])
adj <- sum(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.during")])


h1c <- min(down, adj)

# H2A 
down <- sum(apply(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.before", "Mean.during")],
                1, all))
down <- sum(apply(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.before", "Mean.during")],
                 1, all))

h2a <- min(down, adj)


# H2B

down <- sum(apply(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.before", "Mean.during")],
                 1, all))

h2b <- down

# H2C 
down <- sum(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.during")])
down <- sum(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.during")])


h2c <- min(down, down)

# H2Aprime

adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.before", "Mean.during")],
                1, all))
down <- sum(apply(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.before", "Mean.during")],
                  1, all))

h2aprime <- min(down, adj)


# H2Cprime

adj <- sum(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.during")])
down <- sum(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.during")])

h2cprime <- min(adj, down)

# H3A: need mean and sd (and N, assume if have sd also have N). During and before, upstream and adjacent.
down <- sum(apply(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.before", "Mean.after")],
                1, all))
adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.before", "Mean.after")],
                 1, all))

h3a <- min(down, adj)

# H3B
adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.before", "Mean.after")],
                 1, all))
h3b <- adj

# H3C 
up <- sum(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.after")])
adj <- sum(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.after")])


h3c <- min(up, adj)

# H3Aprime

adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.after", "Mean.during")],
                 1, all))
up <- sum(apply(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.after", "Mean.during")],
                  1, all))

h3aprime <- min(up, adj)


# H3Bprime

adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.during", "Mean.after")],
                 1, all))

h3bprime <- adj


studies <- data.frame(h1a, h1b, h1c, h2a, h2b, h2c, h2aprime, 
                      h2cprime, h3a, h3b, h3c, h3aprime, h3bprime)

means.loc <- 
  aggregate(dat[,c("Mean.before", "Mean.during", "Mean.after")],
            by = list(dat$PublicationID, dat$location),
            FUN = function(x) length(x[!is.na(x)])>0)

sd.loc <- 
  aggregate(dat[,c("SD.before", "SD.during", "SD.after")],
            by = list(dat$PublicationID, dat$location),
            FUN = function(x) length(x[!is.na(x)])>0)

# H1A: need mean and sd (and N, assume if have sd also have N). During and before, upstream and adjacent.
down <- sum(apply(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.before", "Mean.during")],
                  1, all))
adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.before", "Mean.during")],
                 1, all))

h1a <- min(down, adj)

# H1B
adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.before", "Mean.during")],
                 1, all))

h1b <- adj

# H1C 
down <- sum(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.during")])
adj <- sum(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.during")])


h1c <- min(down, adj)

# H2A 
down <- sum(apply(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.before", "Mean.during")],
                  1, all))
down <- sum(apply(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.before", "Mean.during")],
                  1, all))

h2a <- min(down, adj)


# H2B

down <- sum(apply(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.before", "Mean.during")],
                  1, all))

h2b <- down

# H2C 
down <- sum(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.during")])
down <- sum(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.during")])


h2c <- min(down, down)

# H2Aprime

adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.before", "Mean.during")],
                 1, all))
down <- sum(apply(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.before", "Mean.during")],
                  1, all))

h2aprime <- min(down, adj)


# H2Cprime

adj <- sum(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.during")])
down <- sum(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.during")])

h2cprime <- min(adj, down)

# H3A: need mean and sd (and N, assume if have sd also have N). During and before, upstream and adjacent.
down <- sum(apply(means.loc[means.loc$Group.2 == "3 downstream",c("Mean.before", "Mean.after")],
                  1, all))
adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.before", "Mean.after")],
                 1, all))

h3a <- min(down, adj)

# H3B
adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.before", "Mean.after")],
                 1, all))
h3b <- adj

# H3C 
up <- sum(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.after")])
adj <- sum(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.after")])


h3c <- min(up, adj)

# H3Aprime

adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.after", "Mean.during")],
                 1, all))
up <- sum(apply(means.loc[means.loc$Group.2 == "1 upstream",c("Mean.after", "Mean.during")],
                1, all))

h3aprime <- min(up, adj)


# H3Bprime

adj <- sum(apply(means.loc[means.loc$Group.2 == "2 adjacent",c("Mean.during", "Mean.after")],
                 1, all))

h3bprime <- adj


publications <- data.frame(h1a, h1b, h1c, h2a, h2b, h2c, h2aprime, 
                      h2cprime, h3a, h3b, h3c, h3aprime, h3bprime)


summarytable <- data.frame(t(rbind(publications, studies)))


names(summarytable) <- c("Publications", "Studies")

print(summarytable)




