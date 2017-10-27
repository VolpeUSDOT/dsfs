# Re-creating Don's vote count ..
library(scales) # for alpha()

setwd("//vntscex/dfs/3BC-Share/NHTSA/FY15/Meta-Analysis/Statistical Summary Tables")


d <- read.csv("dsfs __ nhtsa __ vote count __ 3 23 17.csv")

# For each safety focus: 
# Make a table of Safety Focuse, Mean Speed, 85th Percentile, and Over Speed Limit
# Rows: Overall proportion, and then count of studies within each hypothesis H1, H2, and H3.


table(d$Safety.Focus)

# Now make the proportion tables

h1 <- d[c("H1.Decrease","H1.NC","H1.Increase")]
h1.85 <- d[c("H1.85.Decrease","H1.85.Increase")]
h1.pct <- d[c("H1.Pct.Decrease","H1.Pct.Increase")]

h2 <- d[c("H2.Decrease","H2.No.Change","H2.Increase")]
h2.85 <- d[c("H2.85.Decrease","H2.85.Increase")]
h2.pct <- d[c("H2.Pct.Decrease","H2.Pct.Increase")]

h3 <- d[c("H3.Decrease","H3.Increase")]
h3.85 <- d[c("H3.85.Decrease","H3.85.Increase")]
h3.pct <- d[c("H3.Pct.Decrease","H3.Pct.Increase")]

makesummarytable <- function(x){
  
  x$Total <- rowSums(x, na.rm = T)
  
  # Aggregate to totals across all publications
  x <- aggregate(x, by=list(Safety.Focus = d$Safety.Focus), sum, na.rm=T)
  
  # find the 'decrease' column
  Prop <- x[,grep("Decrease", names(x))]/x$Total
  No.Decrease = x$Total - x[,grep("Decrease", names(x))]
  
  data.frame(Safety.Focus = x[,1], 
             Decrease = x[,grep("Decrease", names(x))],
             No.Decrease = No.Decrease, Prop = Prop)
  }


result <- rbind(data.frame(Hypothesis = "H1", makesummarytable(h1)),
      data.frame(Hypothesis = "H1 85", makesummarytable(h1.85)),
      data.frame(Hypothesis = "H1 Pct Speed", makesummarytable(h1.pct)),
      data.frame(Hypothesis = "H2", makesummarytable(h2)),
      data.frame(Hypothesis = "H2 85", makesummarytable(h2.85)),
      data.frame(Hypothesis = "H2 Pct Speed", makesummarytable(h2.pct)),
      data.frame(Hypothesis = "H3", makesummarytable(h3)),
      data.frame(Hypothesis = "H3 85", makesummarytable(h3.85)),
      data.frame(Hypothesis = "H3 Pct Speed", makesummarytable(h3.pct))
      )


# Reorder by safety focus, delete blank focus lines
result <- result[result$Safety.Focus != "",]

result.saf <- result[order(result$Safety.Focus),]

write.csv(result, file = "Vote Count Summary.csv", row.names=F)
write.csv(result.saf, file = "Vote Count Summary by Safety Focus.csv", row.names=F)

### Re-creating the figures ----
# Three-part barcharts

# For first three, overall hypotheses, sum across safety foci


makesummarytable2 <- function(x){
  
  x$Total <- rowSums(x, na.rm = T)
  
  # Aggregate to totals across all publications
  x <- aggregate(x, by=list(Safety.Focus = d$Safety.Focus), sum, na.rm=T)
  
  if(length(grep("[NC]", names(x))) > 0) { nc <- x[,grep("[NC]", names(x))] } else nc <- rep(0, nrow(x))
   
  data.frame(Safety.Focus = x[,1], 
             Decrease = x[,grep("Decrease", names(x))],
             No.Change = nc,
             Increase =  x[,grep("Increase", names(x))]
             )
}


result2 <- rbind(data.frame(Hypothesis = "H1", makesummarytable2(h1)),
                data.frame(Hypothesis = "H1 85", makesummarytable2(h1.85)),
                data.frame(Hypothesis = "H1 Pct Speed", makesummarytable2(h1.pct)),
                data.frame(Hypothesis = "H2", makesummarytable2(h2)),
                data.frame(Hypothesis = "H2 85", makesummarytable2(h2.85)),
                data.frame(Hypothesis = "H2 Pct Speed", makesummarytable2(h2.pct)),
                data.frame(Hypothesis = "H3", makesummarytable2(h3)),
                data.frame(Hypothesis = "H3 85", makesummarytable2(h3.85)),
                data.frame(Hypothesis = "H3 Pct Speed", makesummarytable2(h3.pct))
)


r2 <- aggregate(result2[3:5], by = list(result2$Hypothesis), FUN = sum)
rownames(r2) = r2[,1]; r2 <- r2[,-1]

pdf(width = 5, height = 3, file = "Vote Count Figs.pdf")

# Figure 10: Activation Hypothesis
par(xpd = T, cex = 0.7, mar = c(2.5, 2, 2, 1))
barcols = alpha(c("dodgerblue3", "lightblue3", "thistle3"), 0.8)

dx <- as.matrix(r2[c("H1","H1 85", "H1 Pct Speed"),])

b1 <- barplot(dx,
        beside = T,
        space = c(0.1,1),
        legend = T,
        col = barcols,
        names.arg = c("Decrease", "No Change", "Increase"),
        main = "H1: Activation Hypothesis",
        args.legend = list(bty = "n",
                           legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
                            )
        
        )

offset = 5
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
      labels = as.vector(dx),
     cex = 0.8)

# Figure 11: Downstream

dx <- as.matrix(r2[c("H2","H2 85", "H2 Pct Speed"),])

b1 <- barplot(dx,
              beside = T,
              space = c(0.1,1),
              legend = T,
              col = barcols,
              names.arg = c("Decrease", "No Change", "Increase"),
              main = "H2: Downstream Hypothesis",
              args.legend = list(bty = "n",
                                 legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
              )
              
)

offset = 2
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
     labels = as.vector(dx),
     cex = 0.8)


# No deactivation hypothesis fig overall.

result2$overallhyp <- substr(result2$Hypothesis, 1, 2)

# Fig 12: workzones H1

dx <- as.matrix(result2[result2$overallhyp == "H1" & result2$Safety.Focus == "Work zone", 3:5])

b1 <- barplot(dx,
              beside = T,
              space = c(0.1,1),
              legend = T,
              col = barcols,
              names.arg = c("Decrease", "No Change", "Increase"),
              main = "Work Zones: H1",
              args.legend = list(bty = "n",
                                 legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
              )
              
)

offset = 2
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
     labels = as.vector(dx),
     cex = 0.8)

# Fig 13: Work zones H2

dx <- as.matrix(result2[result2$overallhyp == "H2" & result2$Safety.Focus == "Work zone", 3:5])

b1 <- barplot(dx,
              beside = T,
              space = c(0.1,1),
              legend = T,
              col = barcols,
              names.arg = c("Decrease", "No Change", "Increase"),
              main = "Work Zones: H2",
              args.legend = list(bty = "n",
                                 legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
              )
              
)

offset = 2
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
     labels = as.vector(dx),
     cex = 0.8)

# Fig 14: Work zones H3
dx <- as.matrix(result2[result2$overallhyp == "H3" & result2$Safety.Focus == "Work zone", 3:5])

b1 <- barplot(dx,
              yaxp = c(0, 3, 3),
              beside = T,
              space = c(0.1,1),
              legend = T,
              col = barcols,
              names.arg = c("Decrease", "No Change", "Increase"),
              main = "Work Zones: H3",
              args.legend = list(bty = "n",
                                 legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
              )
              
)

offset = 0.1
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
     labels = as.vector(dx),
     cex = 0.8)

# Fig 16: School Zone H1
dx <- as.matrix(result2[result2$overallhyp == "H1" & result2$Safety.Focus == "School zone", 3:5])

b1 <- barplot(dx,
              beside = T,
              ylim = c(0, 25),
              space = c(0.1,1),
              legend = T,
              col = barcols,
              names.arg = c("Decrease", "No Change", "Increase"),
              main = "School Zone: H1",
              args.legend = list(bty = "n",
                                 legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
              )
              
)

offset = 1
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
     labels = as.vector(dx),
     cex = 0.8)

# Fig 17: School zone H2
dx <- as.matrix(result2[result2$overallhyp == "H2" & result2$Safety.Focus == "School zone", 3:5])

b1 <- barplot(dx,
              beside = T,
              yaxp = c(0, 3, 3),
              ylim = c(0, 3.5),
              space = c(0.1,1),
              legend = T,
              col = barcols,
              names.arg = c("Decrease", "No Change", "Increase"),
              main = "School Zone: H2",
              args.legend = list(bty = "n",
                                 legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
              )
              
)

offset = 0.1
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
     labels = as.vector(dx),
     cex = 0.8)

# Fig 19: Transition Zone: H1
dx <- as.matrix(result2[result2$overallhyp == "H1" & result2$Safety.Focus == "Transition zone", 3:5])

b1 <- barplot(dx,
              beside = T,
              ylim = c(0, 30),
              space = c(0.1,1),
              legend = T,
              col = barcols,
              names.arg = c("Decrease", "No Change", "Increase"),
              main = "Transition Zone: H1",
              args.legend = list(bty = "n",
                                 legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
              )
              
)

offset = 1
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
     labels = as.vector(dx),
     cex = 0.8)

# Fig 20: Transition Zone: H2
dx <- as.matrix(result2[result2$overallhyp == "H2" & result2$Safety.Focus == "Transition zone", 3:5])

b1 <- barplot(dx,
              beside = T,
              yaxp = c(0, 2, 2),
              ylim = c(0, 2),
              space = c(0.1,1),
              legend = T,
              col = barcols,
              names.arg = c("Decrease", "No Change", "Increase"),
              main = "Transition Zone: H2",
              args.legend = list(bty = "n",
                                 legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
              )
              
)

offset = 0.1
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
     labels = as.vector(dx),
     cex = 0.8)

# Fig 22: Horizontal curve H1
dx <- as.matrix(result2[result2$overallhyp == "H1" & result2$Safety.Focus == "Horizontal curve", 3:5])

b1 <- barplot(dx,
              beside = T,
            #  yaxp = c(0, 2, 2),
             ylim = c(0, 30),
              space = c(0.1,1),
              legend = T,
              col = barcols,
              names.arg = c("Decrease", "No Change", "Increase"),
              main = "Horizontal Curve: H1",
              args.legend = list(bty = "n",
                                 legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
              )
              
)

offset = 1
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
     labels = as.vector(dx),
     cex = 0.8)

# Fig 23: Horizontal curve H2
dx <- as.matrix(result2[result2$overallhyp == "H2" & result2$Safety.Focus == "Horizontal curve", 3:5])

b1 <- barplot(dx,
              beside = T,
              yaxp = c(0, 10, 2),
              ylim = c(0, 11),
              space = c(0.1,1),
              legend = T,
              col = barcols,
              names.arg = c("Decrease", "No Change", "Increase"),
              main = "Horizontal Curve: H2",
              args.legend = list(bty = "n",
                                 legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
              )
              
)

offset = 0.5
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
     labels = as.vector(dx),
     cex = 0.8)

# Fig 25: Straight section H1
dx <- as.matrix(result2[result2$overallhyp == "H1" & result2$Safety.Focus == "Straight zone", 3:5])

b1 <- barplot(dx,
              beside = T,
         #     yaxp = c(0, 10, 2),
         #     ylim = c(0, 11),
              space = c(0.1,1),
              legend = T,
              col = barcols,
              names.arg = c("Decrease", "No Change", "Increase"),
              main = "Straight Section: H1",
              args.legend = list(bty = "n",
                                 legend = c("Mean speed", "85th Percentile Speed", "Percentage speeding")
              )
              
)

offset = 1
text(x = as.vector(b1),
     y = as.vector(dx)+ offset,
     labels = as.vector(dx),
     cex = 0.8)

#  Fig 26: Straight section H2
# All zero -- no sig results from chang nolan nihan in this version.

# Fig 27: Straight section H3
# Also all zero.

dev.off()

# Looking at school and transition zones

dx <- as.matrix(result2[result2$overallhyp == "H1" & result2$Safety.Focus == "Transition zone", 3:5])


result2[result2$Safety.Focus == "Transition zone",]

d[d$Safety.Focus == "School zone",]


d[d$Safety.Focus == "Transition zone",]

d[d$Safety.Focus == "Horizontal curve",]

d[d$Safety.Focus == "Straight zone",]

