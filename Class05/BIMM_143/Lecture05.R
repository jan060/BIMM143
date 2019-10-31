---
#' title: "Class05 Data Exploration and Visualization in R"
#' author: "Julie Nguyen"
#' output: github_document
---


###Section 2A: Line Plot###
weight <- read.table("Lecture05/bimm143_05_rstats/weight_chart.txt", header = TRUE)
head(weight)
plot(weight$Age, weight$Weight,type = "o", pch = 15,cex = 1.5, lwd = 2,ylim = c(2, 10), xlab = "Age (months)", ylab = "Weight (kg)", main = "Baby Weight with Age",col = "blue")


###Section 2B: Barplot###
mouse <- read.delim("Lecture05/bimm143_05_rstats/feature_counts.txt")
#Set the margin sizes in a vector c(bot, left, top, right)
par(mar = c(5, 12, 4, 5)) 
barplot(mouse$Count,
        horiz = TRUE, #Horizontal bar plot
        names.arg = mouse$Feature, #Names to be plotted below each bar
        main = "Features of Mouse GRCm38 Genome",
        las = 1, #Labels are always horizontal; Can pass as argument in par()
        xlim = c(0, 80000),
        xlab = "Amount")


###Section 2C: Histograms###
x <- c(rnorm(1000), rnorm(1000)+4)
hist(x)


###Section 3A: Providing color vectors###
mf <- read.delim("Lecture05/bimm143_05_rstats/male_female_counts.txt")
par(mar = c(7, 5, 4, 2))
#Create rainbow bar plot
barplot(mf$Count, 
        col = rainbow(nrow(mf)), 
        names.arg = mf$Sample, 
        las = 2, 
        ylab = "Counts")
#Create bar plot with colors corresponding to male and female
barplot(mf$Count, 
        col = c("blue3", "red3"), 
        names.arg = mf$Sample, 
        las = 2, 
        ylab = "Counts")


###Section 3B: Coloring by value###
express <- read.delim("Lecture05/bimm143_05_rstats/up_down_expression.txt")
table(express$State)
plot(express$Condition1,
     express$Condition2, 
     col = express$State)
palette() #Allows you to see the colors chosen for each State
levels(express$State) #Allows you to see how the colors were chosen for each State

palette(c("red3", "blue2", "grey"))
plot(express$Condition1, express$Condition2,
     xlab = "Condition 1", ylab = "Condition 2",
     main = "Gene Expression Comparison", 
     col = express$State)
legend(8, 5, legend = c("Up", "Down", "Unchanging"),
       pch = 1,
       col = c("red3", "blue2", "grey"))


###Section 3C: Dynamic use of color###
meth <- read.delim("Lecture05/bimm143_05_rstats/expression_methylation.txt")
dcols <- densCols(meth$gene.meth, meth$expression) #Changes the colors by density
plot(meth$gene.meth, meth$expression, col = dcols, pch = 20)
#Restrict to genes that have > 0 expression levels
inds <- meth$expression > 0
#Create plot with color gradient corresponding to density
dcols <- densCols(meth$gene.meth[inds], meth$expression[inds], 
                  colramp = colorRampPalette(c("blue2", "green1", "red2", "yellow1")))
plot(meth$gene.meth[inds], meth$expression[inds],
     xlab = "Gene Methylation", ylab = "Expression", 
     main = "Effect of Methylation on Expression",
     col = dcols, pch = 20)