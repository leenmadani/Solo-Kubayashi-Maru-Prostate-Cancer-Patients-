
setwd("C:/Users/Leen/Desktop/Mbiotech/BTC1877")

# Read dataset into a dataframe using read.csv function 
data <- read.csv("study2_n1.csv")

# Take a deeper look at the dataset's structure to understand its components using glimpse()
library(dplyr)
glimpse(data)

# Select the variables that we are interested in and filter for cohort A  
selected_data <- data %>%
  select(AGE, COHORT, DXSTAGE, T1CHEMO, T1DIABET, T1EQTOT, T1KIDNEY, T1RP, T1RT, T1SPBONE, T1VIAGRA, ptid, t1ADT, t1arthritis, t1heart, t1sporg, t2ADT, t2arthritis, t2rp, t2rt, t3ADT, t3arthritis, t3rp, t3rt, T2EQTOT, T3EQTOT) %>% 
  filter(COHORT == "A")



# Dealing with missing values 
any(is.na(selected_data)) # is there any missing values portrayed as NA?
# Let's visualize our missing data using nania package
#install.packages("naniar")
library(naniar)
vis_miss(selected_data) # shows 0% missingess
# there isn't but we need to double check


# Maybe na's are found as character strings and not as the natural encoding of missingness in R as NA 
na_check <- apply(selected_data, 1, function(x) any(grepl("na", x, ignore.case = TRUE)))
which(na_check)
na_observations <- selected_data[na_check, ] # 3 are na as character 

# Check emptiness 
empty_entries <- which(apply(selected_data, 1, function(x) any(x == "")))
# Subset the original dataset to display the entire observations with empty entries
rows_empty <- selected_data[empty_entries, ]
# 18 observations have empty entries 


#Replace both character strings "na" and empty entries with NA that is encoded by R 
data_na <- selected_data %>%
  mutate_all(~ifelse(. %in% c("na", ""), NA, .))

vis_miss(data_na) 

# 15% of DXSTAGE is missing (0.6% of the total dataset missingness)

glimpse(data_na)

table(data_na$DXSTAGE)

# Replace "MO" with "M0" in the DXSTAGE column
data_na$DXSTAGE <- gsub("MO", "M0", data_na$DXSTAGE)

# Now 2 observations have T3BNXM0 because previously the MO was 0. 

#Load necessary packages for descriptive stats and EDA
library(funModeling)
library(tidyverse)
library(Hmisc)


# Let's encode categorical variables as factors: 

data_categ <- data_na %>%
  mutate(
    COHORT = as.factor(COHORT),
    DXSTAGE = as.factor(DXSTAGE),
    T1CHEMO = as.factor(T1CHEMO),
    T1DIABET = as.factor(T1DIABET),
    T1KIDNEY = as.factor(T1KIDNEY),
    T1RP = as.factor(T1RP),
    T1RT = as.factor(T1RT),
    T1SPBONE = as.factor(T1SPBONE),
    T1VIAGRA = as.factor(T1VIAGRA),
    t1ADT = as.factor(t1ADT),
    t1arthritis = as.factor(t1arthritis),
    t1heart = as.factor(t1heart),
    t1sporg = as.factor(t1sporg),
    t2ADT = as.factor(t2ADT),
    t2arthritis = as.factor(t2arthritis),
    t2rp = as.factor(t2rp),
    t2rt = as.factor(t2rt),
    t3ADT = as.factor(t3ADT),
    t3arthritis = as.factor(t3arthritis),
    t3rp = as.factor(t3rp),
    t3rt = as.factor(t3rt)
  )

# The following basic_eda function is adapted from datascienceheroes.com
basic_eda <- function(data)
{
  # Remove patient id column
  data <- data[, colnames(data) != "ptid"]
  
  glimpse(data)
  print(status(data))
  print(profiling_num(data))
  plot_num(data, path_out = ".") # to export plot as jpeg to current directory, add path_out = "."
  describe(data)
}

basic_eda(data_categ)

# CHecking normality of for continuous variables: 

# Q-Q plot and Shapiro-Wilk test for AGE
qqnorm(data_categ$AGE)
qqline(data_categ$AGE)
shapiro.test(data_categ$AGE)

# Q-Q plot and Shapiro-Wilk test for T1EQTOT
qqnorm(data_categ$T1EQTOT)
qqline(data_categ$T1EQTOT)
shapiro.test(data_categ$T1EQTOT)

# Recall for saphiro wilk test null hypothesis is that data follows normality
#Based on these results, both "AGE" and "T1EQTOT" are not normally distributed as p-value is lower than 0.05 (t1eqtot signficantly lower) so reject null hypothesis. 

anova_t1eqtot <- aov(T1EQTOT ~ DXSTAGE, data = data_categ)

# Summarize the ANOVA results
summary(anova_t1eqtot)

#Create data frames for missing and complete DXSTAGE data
data_missing_dxstage <- data_categ %>% filter(is.na(DXSTAGE))
data_complete_dxstage <- data_categ %>% filter(!is.na(DXSTAGE))

#Perform a Wilcoxon rank-sum test to compare the distribution of numerical variables (e.g., AGE) between the groups with and without missing DXSTAGE data.
wilcox.test(data_missing_dxstage$AGE, data_complete_dxstage$AGE)
wilcox.test(data_missing_dxstage$T1EQTOT, data_complete_dxstage$T1EQTOT)

# it appears that age is not significantly different between patients with missing and complete DXSTAGE data 


######### Let's check if theres relationship between missing dxstage dataset and comorbid conditions. Literature says that 
# those with  higher comorbidity index are more likely to havig missing staged data. 

# Create a new variable for each comorbid condition where 1 and 2 indicate "yes" (1)
comorbid_conditions <- c("T1DIABET", "T1KIDNEY", "t1arthritis", "t1heart", "t2arthritis", "t3arthritis")

# Calculate the proportions of each comorbid condition in the missing and complete groups
proportions_missing <- sapply(comorbid_conditions, function(condition) {
  sum(data_missing_dxstage[[condition]] %in% c(1, 2)) / nrow(data_missing_dxstage)
})

proportions_complete <- sapply(comorbid_conditions, function(condition) {
  sum(data_complete_dxstage[[condition]] %in% c(1, 2)) / nrow(data_complete_dxstage)
})

# Combine the proportions into a data frame for easy comparison
comparison <- data.frame(Comorbid_Condition = comorbid_conditions, 
                         Proportion_Missing = proportions_missing, 
                         Proportion_Complete = proportions_complete)

# Save the comparison data frame as a CSV file
write.csv(comparison, file = "comorbid_condition_comparison.csv", row.names = FALSE)



###################

###### HANDLING MISSING DATA!!!!!

# Extract the first two letters from each stage label
T_stage <- substr(data_categ$DXSTAGE, 1, 2)

# Count the occurrences of each unique first two letters
Tstage_counts <- table(T_stage)


# Extract the first two letters from each stage label
T_stage_impute <- substr(data_impute$DXSTAGE, 1, 2)

# Count the occurrences of each unique first two letters
Tstage_count_impute <- table(T_stage_impute)

# Impute DXSTAGE based on T1SPBONE
data_impute <- data_categ %>%
  mutate(
    DXSTAGE = case_when(
      T1SPBONE == 1 & is.na(DXSTAGE) ~ "T4",
      t1sporg == 1 & is.na(DXSTAGE) ~ "T4",
      t1sporg == 0 & T1RP == 0 & T1RT == 0 & is.na(DXSTAGE) ~ "T1",
      t1sporg == 0 & (T1RP == 1 | T1RT == 1) & is.na(DXSTAGE) ~ "T2",
      TRUE ~ DXSTAGE  # Keep existing DXSTAGE values for other cases
    )
  )

vis_miss(data_categ)
vis_miss(data_impute)

na_left <- data_impute[is.na(data_impute$DXSTAGE) == T,]

data_impute <- na.omit(data_impute)
basic_eda(data_impute)

library(knitr)
##### Visualize cont. variables 
table1 <- profiling_num(data_impute)
sd(data_impute$AGE)
sd(data_impute$T1EQTOT)

##### Visualize the categorical variables
# Define a vector of categorical variables
categorical_vars <- c("T1DIABET", "T1KIDNEY", "T1RP", "T1RT", "T1SPBONE", "T1VIAGRA", "t1ADT", "t1arthritis", "t1heart", "t1sporg",
                      "t2ADT", "t2arthritis", "t2rp", "t2rt", "t3ADT", "t3arthritis", "t3rp", "t3rt")

# Loop through each categorical variable and create bar plots 
for (var in categorical_vars) {
  p <- ggplot(data_impute, aes(x = factor(data_impute[[var]]))) +
    geom_bar(aes(fill = factor(data_impute[[var]])),  show.legend = FALSE) +
    labs(x = var, y = "Frequency") +
      scale_fill_brewer(palette = "Set2") + 
    theme_minimal()
  print(p)  # Display the plot
}


##### Now lets analyze the impact of the different treatments!
t1 <- aov(T1EQTOT ~ T1RP + T1RT + t1ADT, data = data_impute)
summary(t1)



# Q-Q plot and Shapiro-Wilk test for transformed T1EQTOT

shapiro.test(log(data_impute$T1EQTOT))

hist(log(data_impute$T1EQTOT))


# Reflect & Take the square root of the reflected variable
T1EQTOT_reflect_sqrt <- sqrt(max(data_impute$T1EQTOT) + 1 - data_impute$T1EQTOT)
qqnorm(sqrt(max(data_impute$T1EQTOT) + 1 - data_impute$T1EQTOT))
qqline(sqrt(max(data_impute$T1EQTOT) + 1 - data_impute$T1EQTOT))

# Set custom plot margins
par(mar = c(4, 4, 2, 2)) 
# Check if the transformed variable is more normally distributed
hist(T1EQTOT_reflect_sqrt, main = "Histogram of sqrt(reflect(T1EQTOT))")

# Square root
hist(sqrt(data_impute$T1EQTOT))
qqnorm(sqrt(data_impute$T1EQTOT))
qqline(sqrt(data_impute$T1EQTOT))

#Reflect then natural log
hist(log(max(data_impute$T1EQTOT) + 1 - data_impute$T1EQTOT))
qqnorm(log(max(data_impute$T1EQTOT) + 1 - data_impute$T1EQTOT))
qqline(log(max(data_impute$T1EQTOT) + 1 - data_impute$T1EQTOT))


hist(log10(max(data_impute$T1EQTOT) + 1 - data_impute$T1EQTOT))
qqnorm(log10(max(data_impute$T1EQTOT) + 1 - data_impute$T1EQTOT))
qqline(log10(max(data_impute$T1EQTOT) + 1 - data_impute$T1EQTOT))


# Kruskal-Wallis test for T1, T2, and T3
kruskal.test(T1EQTOT ~ T1RP, data = data_impute)
kruskal.test(T1EQTOT ~ T1RT, data = data_impute)
kruskal.test(T1EQTOT ~ t1ADT, data = data_impute)
kruskal.test(T1EQTOT ~ t1arthritis, data = data_impute)
kruskal.test(T1EQTOT ~ t1arthritis, data = data_impute)

# Kruskal-Wallis test for T2
kruskal.test(T2EQTOT ~ t2rp, data = data_impute)
kruskal.test(T2EQTOT ~ t2rt, data = data_impute)
kruskal.test(T2EQTOT ~ t2ADT, data = data_impute)
pairwise.wilcox.test(data_impute$t3ADT, data_impute$t2eqtot, data_impute$t2rp, p.adjust.method = "BH")

# Kruskal-Wallis test for T3
kruskal.test(T3EQTOT ~ t3rp, data = data_impute)
kruskal.test(T3EQTOT ~ t3rt, data = data_impute)
kruskal.test(T3EQTOT ~ t3ADT, data = data_impute)

# Perform pairwise Wilcoxon rank-sum tests for t3ADT and t3eqtot by t3rp
pairwise.wilcox.test(data_impute$t3ADT, data_impute$t3eqtot, data_impute$t3rp, p.adjust.method = "BH")

# Fit a linear regression model
lm_modelt1 <- lm(T1EQTOT ~ T1RP + T1RT + t1ADT + T1DIABET + T1KIDNEY + t1arthritis + t1heart + t1sporg + T1VIAGRA + T1SPBONE, data = data_impute)


plot(lm_modelt1$fitted.values, residuals(lm_modelt1))
abline(h = 0, col = "red")

plot(lm_modelt1$fitted.values, residuals(lm_modelt1))

# View the summary of the regression model
summary(lm_modelt1)

write.csv(lmt1, file = "lm_results_T1.csv")


# Fit a linear regression model for T2
lm_modelt2 <- lm(T2EQTOT ~ t2rp + t2rt + t2ADT + t2arthritis, data = data_impute)


plot(lm_modelt2$fitted.values, residuals(lm_modelt2 ))
abline(h = 0, col = "red")

qqnorm(residuals(lm_modelt1))
qqline(residuals(lm_modelt1))


# View the summary of the regression model
summary(lm_modelt2)


# Fit a linear regression model for T3
lm_modelt3 <- lm(T3EQTOT ~ t3rp + t3rt + t3ADT + t3arthritis, data = data_impute)

plot(lm_modelt3$fitted.values, residuals(lm_modelt3))
abline(h = 0, col = "red")
# View the summary of the regression model
summary(lm_modelt3)


qqnorm(residuals(lm_modelt2))
qqline(residuals(lm_modelt2))

qqnorm(residuals(lm_modelt3))
qqline(residuals(lm_modelt3))


########### DIAGNOSIS STAGING 

# Extract the first two characters from DXSTAGE
data_impute$DXSTAGE_NEW <- substr(data_impute$DXSTAGE, 1, 2)

# Convert DXSTAGE_NEW to a factor
data_impute$DXSTAGE_NEW <- factor(data_impute$DXSTAGE_NEW)


# Fit a linear regression model with the new DXSTAGE_NEW variable
lm_model <- lm(T1EQTOT ~ T1RP + T1RT + t1ADT + DXSTAGE_NEW, data = data_impute)

# View the summary of the regression model
summary(lm_model)


# Create a summary table of frequencies for DXSTAGE_NEW
frequency_table <- table(data_impute$DXSTAGE_NEW)

# Convert the table to a data frame
frequency_df <- as.data.frame(frequency_table)
colnames(frequency_df) <- c("DXSTAGE_NEW", "Frequency")

# Create the bar plot
ggplot(data_impute, aes(x = factor(DXSTAGE_NEW))) +
  geom_bar(fill = "lightblue") +  # You can change the fill color
  labs(x = "Diagnosis Stage (DXSTAGE_NEW)", y = "Frequency") +
  theme_minimal()


dev.off()