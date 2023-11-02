#Step 1: Data Retrieval
#We begin by downloading the DLBCL dataset and its associated class labels. 
#The dataset is loaded into R, and missing values are accounted for using the read.table function.

#Assign the path to the eisen.txt file and read it into a data frame
data <- read.table("eisen.txt", 
                   header = TRUE, na.strings = "NA", blank.lines.skip = FALSE, row.names = 1)

#Step 2: Data Preprocessing
#We process the dataset by reordering the columns to match the class labels. 
#This allows us to separate data related to different DLBCL classes for further analysis.

#Subset the data frame with class labels
cl <- as.character(class_labels[, 2])  
data <- data[, cl]


#Step 3: Data Exploration
#To better understand the data, we select a gene of interest, remove any missing values, 
#and visualize the gene expression differences between two DLBCL classes using boxplots and histograms. 
#This helps us identify the gene's characteristics within the dataset.

#Splitting up the classes and look for the gene of interest
x <- as.numeric(data[8000, GC_group])
y <- as.numeric(data[8000, ACT_group])

#Remove all NA's
x <- x[!is.na(x)]; y <- y[!is.na(y)]

#Plotting a boxplot
boxplot(list(GC_group = x, ACT_group = y), col = c("red", "blue"), 
        main = 'Example gene from DLBCL cDNA 2-channel dataset', 
        ylab = "log2(ratio intensity)",
        names = c("GC_group", "ACT_group"))

#Plot individual histograms for both x and y
hist(x, xlab = 'Selected gene', ylab = 'GC group', 
     main = "Gene selected from DLBCL cDNA 2-channel dataset", 
     col = c("red"))

hist(y, xlab = 'Selected gene', ylab = 'ACT group', 
     main = "Gene selected from DLBCL cDNA 2-channel dataset", 
     col = c("blue"))

#Step 4: Power and Sample Size Calculation - Single Gene
#We calculate the pooled variance for the selected gene, which is a crucial parameter for power and sample 
#size calculations. We then compute the minimum sample size needed to detect a 1.5-fold difference with 80% 
#power and 99% confidence. This is accomplished using the pwr library for power analysis.

#Calculate pooled variance
nx <- length(x)
ny <- length(y)
pool.var <- (((nx - 1) * var(x)) + ((ny - 1) * var(y))) / (nx + ny - 2)


#Step 5: Power and Sample Size Calculation - Multivariate Analysis
#We calculate the sample size required to detect a difference empirically for the 
#same gene using the determined delta between the two classes, assuming 99% confidence and 80% power. This further 
#ensures that our analysis is robust and generalizable.

#Load the pwr library for power analysis
library(pwr)

#Calculate minimum sample size necessary to detect a 1.5-fold difference for the dataset
dif.1.5fold <- log2(1.5) / sqrt(pool.var)
pl.ss3 <- pwr::pwr.t.test(d = dif.1.5fold, sig.level = 0.01, power = 0.8, type = "two.sample")

#Calculate the difference (dif) based on your data
dif <- abs(mean(x) - mean(y)) / sqrt(pool.var)

#Define your desired significance level, power, and type of test
sig.level <- 0.01
power <- 0.8
test_type <- "two.sample"

# Calculate the sample size for the new dataset
pl.ss31 <- pwr::pwr.t.test(d = dif, sig.level = sig.level, power = power, type = test_type)

#Step 6: Standard Deviation Calculation
#To prepare for a proportion of genes vs. sample size analysis, we calculate the standard deviation 
#for each gene in the dataset, considering missing values. We use the ssize and gdata libraries 
#for this purpose.

#Calculate standard deviation for each gene in new_data
exp.sd <- apply(data_new, 1, sd, na.rm = TRUE)

#Plot histogram of the variable
hist(exp.sd, n = length(exp.sd), col = "red", border = "pink", main = "", xlab = "
