# Date: 22 July 2021
# Quiz 01: Gene Expression Analysis
# Name: # TODO

# Instructions: You are to follow the analysis detailed below in code. From time-to-time, you will encounter questions. Please be sure to add all your responses to these questions in the report.md file of your writing/ directory in your repository. Please see details of the assignment sheet.



rm(list = ls()) # remove all variables in memory.

# load the libraries
library(tidyverse)
library(dplyr)

#install.packages("mosaic") # you will likely have to install this new library to load the Mosaic library
library(mosaic)


# Select file NCI60.rda from your repository directory data/ when prompted. You will see this file as a variable the Global Environment panel of rStudio.

load(file.choose()) # select file data/NCI60.rda

# Select file NCI60cells.rda from your repository directory data/ when prompted. You will see this file as a variable the Global Environment panel of rStudio.

load(file.choose()) # select data/NCI60cells.rda



#The below code block creates new and smaller datasets from the larger ones which came in the files that were just loaded.

Narrow <- NCI60 %>% tidyr::gather(cellLine, expression, -Probe)
CellTypes <- NCI60cells %>% select(cellLine, tissue) %>% mutate(cellLine=gsub("\\:",".", as.character(cellLine)))


# The below code block is reorganizing the code block by combining ("joining") the "CellTypes" data with the "Narrow" set. 
Narrow <- Narrow %>% inner_join(CellTypes)

##############################################
# Question 1: What new column has been added as a result of applying the join function to the "Narrow" dataset and CellTypes data sets in the above line of code?


# The below code block is isolating the TOP3A types of probe samples
Probe_TOP3A <- Narrow %>% filter(Probe == "TOP3A")

# The below code block is finding the mean of each type of each type of cancer tissue
SummaryStats <- Probe_TOP3A %>% group_by(tissue) %>% summarise(mn_exp = exp(mean(expression, na.rm=TRUE)))




##############################################
# Question 2: Plot 1 is created below. In terms of the x-axis and y-axes (tissue and gene expression, respectively), offer an interpretation of what the below plot is describing about tissue types. 

plot(SummaryStats, main ="Plot 1 ") 


##############################################
# Question 3: Plot 2 (created below) is a bar-chart view of the same information from Plot 1. How is this plot more informative than the previous Plot 1? 

#A new plot
SummaryStats %>% ggplot(aes(x = tissue, y = mn_exp)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle =30, hjust=1)) + ggtitle("Plot 2")


# The following code creates a new data set from existing datasets.

SummaryStats <- SummaryStats %>% mutate(tissue = reorder(tissue, mn_exp))

SummaryStats <- Probe_TOP3A %>% group_by(tissue) %>% summarise(mn=mean(expression, na.rm=TRUE), se = (sd(expression, na.rm=TRUE) / sqrt(n())))

##############################################
# Question 4: What is the function "mutate()" doing in the above code block?



##############################################
# Question 5: Inside the above code block, we note that there is "na.rm=TRUE". What does this code do?



##############################################
# Question 6: In terms of the x-axis and y-axes (tissue and gene expression, respectively), interpret the results of Plot 3, created below.

SummaryStats %>% ggplot(aes(x = tissue, y = exp(mn))) +  geom_point(data = Probe_TOP3A, aes(x = tissue, y = exp(expression))) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ggtitle("Plot 3")


##############################################
# Question 7: Interpret the results of Plot 4, created by the code block below.


# plot a histogram
SummaryStats %>% ggplot(aes(x = tissue, y = exp(mn))) + geom_bar(stat = "identity", fill = "grey", color = NA) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ggtitle("Plot 4")


##############################################
# Question 8: In Plot 5, both of the previous plots have been combined together on the same canvas. How does this, either help or confuse the output of these lines of code?

# The below code places both plots together on same canvas
SummaryStats %>% ggplot(aes(x = tissue, y = exp(mn))) + geom_bar(stat = "identity", fill = "grey", color = NA) + geom_point(data = Probe_TOP3A, aes(x = tissue, y = exp(expression))) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ggtitle("Plot 5")


# The below code adds some error bars for the next plot
SummaryStats <- SummaryStats %>% mutate(top = mn * 2 * se , bottom = mn - 2 * se)




##############################################
# Question 9: What do the whiskers (i.e., the big "I" shaped structures above each histogram bar) convey in this plot? Why might we need to know this information?

SummaryStats %>% ggplot(aes(x = tissue, y = exp(mn),color= tissue)) + geom_bar(stat = "identity", alpha = 0.2) + geom_errorbar(aes(x = tissue, ymin = exp(bottom), ymax = exp(top)), width = .5) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ggtitle("Plot 7")


# Define a function to create models to aquire R-squared values. This function uses the "mosaic" library.

##############################################
# Question 10: The word, "function" is below in the code block. What does the "lm()" part of the code block do? 

r2 <- function(data){
  mosaic::rsquared(lm(data$expression ~ data$tissue))
}


# In the following code block creates studies cases to describe relationships between anti-cancer drugs and probes (tests). The do() function is analogous to summarise() and the unlist() functions, which facilitate plotting by integrating the data as vectors.

# Note: The next line may take some time to complete due to making models out of all groups which were created from the data.
ProbeR2 <- Narrow %>% group_by(Probe) %>% do(r2 = r2(.)) %>% mutate(r2 = unlist(r2))


##############################################
# Next, we would like to pull out the 30 probes with the largest R-squared values to be plotted. In terms of the test, these values indicate that the indepentant variables were able to successfully predict the dependent variables. In this test, these probes responded well (statistically) to tests.

#  Plot 8 is showing R-squared values which have been extracted from a linear model of each Probe. Although we have not yet discussed this form of analysis (regression analysis) in class, it is sufficient to say that R-squared values provide a statistic to determine a fitness of a model. For instance, the R-squared statistic provides a measurement of how much of the variation in a variable, known as a response variable, is accounted for in other explanatory variables. R-squared values extend from 0 to 1 where a zero-value represents a poor prediction by the explanatory variables, where-as a score of 1 indicates a perfect prediction by the explanatory variables. 


# The following code block creates a plot to show R-squared ("r2") values of each of the probes (the individual tests).
Actual <- ProbeR2 %>% arrange(desc(r2)) %>% head(30) %>% mutate(Probe = reorder(Probe, desc(r2)))

Actual %>% ggplot(aes(x = Probe, y = r2))  + geom_point()  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Plot 8")


##############################################
# Question 11:  With respect to the above discussion of R-squared variables, which Probe in Plot 8 appears to be able to best explain its model?



# Next, we would like to determine whether the top 30 Probes that we plotted in Plot 8 were actually not relevant or were false positives. Now, we will build a Null Hypothesis model out of a set of R-squared values to be used for a comparison. In Plot 8, we argued that those points described relationships between anti-cancer drugs and Probes from tissues. Now, we are going to try to build a model from our data in which there are NO relationships described by the anti-cancer drugs and the Probes in tissues. Then we will test that all the relationships that we plotted in Plot 8, are not included in a plot showing that NO relationships exist.


# The below code creates null hypotheses cases from the Narrow data set to study.  

# Note: The next line may take some time to complete due to making models out of all groups which were created from the data.

NullR2 <- Narrow %>% group_by(Probe) %>% mutate(expression = mosaic::shuffle(expression)) %>% group_by(Probe) %>% do(r2 = r2(.)) %>% mutate(r2 = unlist(r2))



##############################################
# Question 12: In Plot 9 we see two curves. What can you infer about the relationships between the Actual cases of relationships (blue curve) and those cases for which no relationship exists (red curve)? Remember, the blue curve corresponds to the variable "r2" and the red curve corresponds to the Null values.

ProbeR2 %>% ggplot(aes(x = r2)) + geom_density(fill = "blue", color = "blue4", alpha = .50) + geom_density(data = NullR2, aes(x = r2), fill = "pink", alpha = 0.45, color = "red") + ggtitle("Plot 9")
#Note: Actual cases are blue, Null cases are red.



Null <- NullR2 %>% arrange(desc(r2)) %>% head(30) 
Actual$null <- Null$r2


##############################################
# Question 13: Plot 10 created below, is very important to determine a result concerning the Actual (blue) data sets and the Null (red) data sets. In terms of the distribution of the red and blue points, explain what this plot, created using the below code block, describes about the two groups. 

Actual %>% ggplot(aes(x = Probe, y = r2)) + geom_point(color = "blue") + geom_point(aes(y = null), color = "red") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + ggtitle("Plot 10")

#Note: Actual cases are blue, Null cases are red.

