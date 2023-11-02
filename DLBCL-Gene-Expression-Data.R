# Step 1: Data Acquisition and Preprocessing

# Assign the path to the eisen.txt file and read it into a data frame
data <- read.table("eisen.txt", header = TRUE, na.strings = "NA", blank.lines.skip = FALSE, row.names = 1)

# Assign the path to eisenClasses.txt and read it into a data frame for class labels
class_labels <- read.table("eisenClasses.txt", header = TRUE)

# Check the first few rows of the data frame and class_labels data frame
head(data)
head(class_labels)


# Step 2: Data Subset and Exploration

# Subset the data frame with class labels
cl <- as.character(class_labels[, 2])  
data <- data[, cl]

# Rearrange the columns
dimnames(data)[[2]]

# Split up classes and identify the gene of interest 
GC_group <- 1:19 
ACT_group <- 20:39

# Splitting up the classes and look for the gene of interest
x <- as.numeric(data[8000, GC_group])
y <- as.numeric(data[8000, ACT_group])

# Remove all NA's
x <- x[!is.na(x)]; y <- y[!is.na(y)]


# Step 3: Data Visualization

# Plot a boxplot
boxplot(list(GC_group = x, ACT_group = y), col = c("red", "blue"), 
        main = 'Example gene from DLBCL cDNA 2-channel dataset', 
        ylab = "log2(ratio intensity)",
        names = c("GC_group", "ACT_group"))

# Plot individual histograms for both x and y
hist(x, xlab = 'Selected gene', ylab = 'GC group', 
     main = "Gene selected from DLBCL cDNA 2-channel dataset", 
     col = c("red"))

hist(y, xlab = 'Selected gene', ylab = 'ACT group', 
     main = "Gene selected from DLBCL cDNA 2-channel dataset", 
     col = c("blue"))


# Step 4: Pooled Variance and Sample Size Calculation

# Calculate pooled variance
nx <- length(x)
ny <- length(y)
pool.var <- (((nx - 1) * var(x)) + ((ny - 1) * var(y))) / (nx + ny - 2)

# Step 5: Sample Size Calculation

# Load the pwr library for power analysis
library(pwr)

# Calculate minimum sample size necessary to detect a 1.5-fold difference
dif.1.5fold <- log2(1.5) / sqrt(pool.var)
pl.ss3 <- pwr::pwr.t.test(d = dif.1.5fold, sig.level = 0.01, power = 0.8, type = "two.sample")

# Print the pooled variance and minimum sample size
cat("Pooled Variance:", pool.var, "\n")
cat("Minimum Sample Size to Detect 1.5-Fold Difference:", pl.ss3$n, "\n")

# Step 6: Sample Size Calculation with Empirical Data

# Calculate the difference (dif) based on your data
dif <- abs(mean(x) - mean(y)) / sqrt(pool.var)

# Define your desired significance level, power, and type of test
sig.level <- 0.01
power <- 0.8
test_type <- "two.sample"

# Calculate the sample size for the new dataset
pl.ss31 <- pwr::pwr.t.test(d = dif, sig.level = sig.level, power = power, type = test_type)

# Print the calculated difference and sample size
cat("Difference (dif):", dif, "\n")
cat("Sample Size (pl.ss31$n):", pl.ss31$n, "\n")


# Step 7: Standard Deviation and Sample Size Calculation

# Load the ssize library for sample size calculations
library(ssize)

# Load the gdata library for data manipulation
library(gdata)

# Read the data from eisen.txt again
data_new <- read.table("C:/Users/sweth/OneDrive/Documents/R learning/Lab3/eisen.txt", 
          header = TRUE, na.strings = "NA", blank.lines.skip = FALSE, row.names = 1)

# Calculate standard deviation for each gene in the matrix
sd_df <- apply(data_new, 1, function(row) sd(row, na.rm = TRUE))

# Plot a histogram of the standard deviations
hist(sd_df, n = length(sd_df), col = "red", border = "pink", main = "", xlab = "Standard Deviation (for data on the log2 scale)")


# Step 8: Sample Size vs. Proportion of Genes

fold.change <- 3.0
power <- 0.8
sig.level <- 0.05

# Calculate standard deviation for each gene in new_data
exp.sd <- apply(data_new, 1, sd, na.rm = TRUE)

# Calculate sample size based on the new_data
all.size <- ssize(sd = exp.sd, delta = log2(fold.change), sig.level = sig.level, power = power)

# Plot the proportion of genes vs. sample size graph
ssize.plot(all.size, lwd = 2, col = "blue", xlim = c(1, 20))
xmax <- par("usr")[2] - 1
ymin <- par("usr")[3] + 0.05

# Add a legend to the graph
legend(
  x = xmax,
  y = ymin,
  legend = strsplit(paste("fold change =", fold.change, ",", "alpha =", sig.level, ",", "power =", power, ",", "# genes =", nrow(data_new), sep = ','), ",")[[1]],
  xjust = 1,
  yjust = 0,
  cex = 1.0
)

# Set the title for the graph
title("Sample Size to Detect Fold Change")
