# Re-creating Don's vote count ..
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


result <- rbind(data.frame(Hypotheisis = "H1", makesummarytable(h1)),
      data.frame(Hypotheisis = "H1 85", makesummarytable(h1.85)),
      data.frame(Hypotheisis = "H1 Pct Speed", makesummarytable(h1.pct)),
      data.frame(Hypotheisis = "H2", makesummarytable(h2)),
      data.frame(Hypotheisis = "H2 85", makesummarytable(h2.85)),
      data.frame(Hypotheisis = "H2 Pct Speed", makesummarytable(h2.pct)),
      data.frame(Hypotheisis = "H3", makesummarytable(h3)),
      data.frame(Hypotheisis = "H3 85", makesummarytable(h3.85)),
      data.frame(Hypotheisis = "H3 Pct Speed", makesummarytable(h3.pct))
      )


# Reorder by safety focus, delete blank focus lines
result <- result[result$Safety.Focus != "",]

result.saf <- result[order(result$Safety.Focus),]

write.csv(result, file = "Vote Count Summary.csv", row.names=F)
write.csv(result.saf, file = "Vote Count Summary by Safety Focus.csv", row.names=F)
