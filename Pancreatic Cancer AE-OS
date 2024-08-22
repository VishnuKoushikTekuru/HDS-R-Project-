# Load necessary libraries
library(survival)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggfortify)
library(survminer)
library(grid)
library(gridExtra)

# Set the working directory
setwd("/Users/koush/Desktop/restart")

# Load data
osdata <- read.csv("OSdata.csv")
aedata <- read.csv("aeData.csv")

# Clean and prepare data
aedata <- aedata %>%
  select(PersonId, AE, AEgrade, Date) %>%
  mutate(Date = as.Date(Date, "%d/%m/%Y")) %>%
  filter(AE != "")

osdata <- osdata %>%
  mutate(REGDate = as.Date(REGDate, "%d/%m/%Y"),
         OSdt = as.Date(OSdt, "%d/%m/%Y"))

# Filter to keep only AE's with words/alphabets
aedata <- aedata %>%
  filter(!grepl("^[0-9]+$", AE))

# Remove rows with missing or NA values
aedata <- aedata %>%
  filter(!is.na(AE) & !is.na(AEgrade) & !is.na(Date))

# Removing AE's that don't occur very often
tb <- table(aedata$AE)
ae.rem <- names(which(tb < 20))
aedata <- aedata[!aedata$AE %in% ae.rem, ]

# Merging with OSdata
aedata <- merge(aedata, osdata, by = "PersonId")

# Filtering dates to keep only those AEs that occurred between the registration date and last follow-up or death date
aedata <- subset(aedata, Date >= REGDate & Date <= OSdt)

# Prepare data for Cox model
aedata <- aedata %>%
  mutate(Event = OScen,
         SurvivalTime = SurvTime)
# Basic Toxicity Analysis:
# Tables of toxicity between treatment groups
toxicity_table <- aedata %>%
  group_by(Arm, AE, AEgrade) %>%
  summarize(Count=n())

print(toxicity_table)


# Identify the date of the first AE for each person
first_ae_date <- aedata %>%
  group_by(PersonId) %>%
  summarise(First_AE_Date = min(Date, na.rm = TRUE)) %>%
  ungroup()

# View the first few rows of the result
head(first_ae_date)

# Calculating frequency counts for the identified columns
ae_counts <- table(aedata$AE)
cgem_counts <- table(aedata$Cgem)
cvan_counts <- table(aedata$Cvan)
sae_counts <- table(aedata$sae)
ctccat_counts <- table(aedata$CTCcat)

# Print the frequency counts
# Create a list of tables
tables_list <- list(AE = ae_counts, 
                    Cgem = cgem_counts, 
                    Cvan = cvan_counts, 
                    SAE = sae_counts, 
                    CTCcat = ctccat_counts)

# Display the result
print(tables_list)

# Filter out relevant columns for clustering
clustering_data <- osdata[, c('PatAge', 'DEMECOG', 'DASCA19_9', 'SurvTime')]

# Handle missing values using median imputation
clustering_data <- as.data.frame(lapply(clustering_data, function(x) {
  ifelse(is.na(x), median(x, na.rm = TRUE), x)
}))

# Normalize the data
clustering_data_scaled <- scale(clustering_data)

# Determine the optimal number of clusters using the Elbow method
wcss <- numeric()
for (i in 1:10) {
  kmeans_result <- kmeans(clustering_data_scaled, centers=i)
  wcss[i] <- kmeans_result$tot.withinss
}

# Plotting the Elbow method graph
plot(1:10, wcss, type="b", main="Elbow Method", xlab="Number of clusters", ylab="WCSS")

# Apply K-means clustering with 3 clusters
set.seed(42)
kmeans_result <- kmeans(clustering_data_scaled, centers=3)
osdata$Cluster <- kmeans_result$cluster

# Summarize each cluster's characteristics
cluster_summary <- osdata %>%
  group_by(Cluster) %>%
  summarise(
    Avg_Age = mean(PatAge, na.rm=TRUE),
    Age_Std_Dev = sd(PatAge, na.rm=TRUE),
    Avg_DEMECOG = mean(DEMECOG, na.rm=TRUE),
    DEMECOG_Std_Dev = sd(DEMECOG, na.rm=TRUE),
    Avg_DASCA19_9 = mean(DASCA19_9, na.rm=TRUE),
    DASCA19_9_Std_Dev = sd(DASCA19_9, na.rm=TRUE),
    Avg_SurvTime = mean(SurvTime, na.rm=TRUE),
    SurvTime_Std_Dev = sd(SurvTime, na.rm=TRUE)
  )

print(cluster_summary)

write.csv(cluster_summary,"cluster_summary.csv",row.names = FALSE)

# Create the Kaplan-Meier survival object by Treatment Arm
km_fit_treatment <- survfit(Surv(SurvTime, OScen) ~ Arm, data = osdata)

# Plot the Kaplan-Meier survival curves
treatment_km_plot <- ggsurvplot(
  km_fit_treatment,
  data = osdata,
  pval = TRUE,                  # Add p-value to the plot
  conf.int = TRUE,              # Add confidence interval
  risk.table = TRUE,            # Show risk table
  ggtheme = theme_minimal(),    # Set the theme
  title = "Kaplan-Meier Survival Curves by Treatment Arm",
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  palette = c("#E7B800", "#2E9FDF") # Custom colors for the curves
)

# Save the plot to a PDF file
pdf("Kaplan-Meier_Survival_Curves_by_Treatment_Arm.pdf", width = 8, height = 6)
print(treatment_km_plot)
dev.off()

# Extended Cox Model Analysis
# Create a data frame to store results
results <- data.frame(AE = character(),
                      HR = numeric(),
                      CI_lower = numeric(),
                      CI_upper = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each unique AE
for (ae in unique(aedata$AE)) {
  ae_data <- aedata %>%
    mutate(AE_presence = ifelse(AE == ae, 1, 0))
  
  # Fit the Cox model
  cox_model <- coxph(Surv(SurvivalTime, Event) ~ AE_presence, data = ae_data)
  cox_summary <- summary(cox_model)
  
  # Extract the hazard ratio, confidence interval, and p-value
  hr <- round(cox_summary$coefficients[1, "exp(coef)"], 2)
  ci_lower <- round(cox_summary$conf.int[1, "lower .95"], 2)
  ci_upper <- round(cox_summary$conf.int[1, "upper .95"], 2)
  p_value <- round(cox_summary$coefficients[1, "Pr(>|z|)"], 2)
  
  # Append results to the data frame
  results <- rbind(results, data.frame(AE = ae,
                                       HR = hr,
                                       CI_lower = ci_lower,
                                       CI_upper = ci_upper,
                                       p_value = p_value,
                                       stringsAsFactors = FALSE))
}

# Display results
print(results)

# Filter significant results (p < 0.05)
significant_results <- results %>%
  filter(p_value < 0.05)

# Save results to a CSV file
write.csv(results, "cox_model_results.csv", row.names = FALSE)
write.csv(significant_results, "significant_cox_model_results.csv", row.names = FALSE)

# Plot results for all adverse events
p_all <- ggplot(results, aes(x = reorder(AE, -HR), y = HR)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  coord_flip() +
  labs(title = "Hazard Ratios of All Adverse Events",
       x = "Adverse Event",
       y = "Hazard Ratio (95% CI)") +
  theme_minimal()
# Save the plot to a PDF file
pdf("Hazard Ratios of All Adverse Events.pdf", width = 8, height = 6)
print(p_all)
dev.off()

# Plot significant results
p_significant <- ggplot(significant_results, aes(x = reorder(AE, -HR), y = HR)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  coord_flip() +
  labs(title = "Hazard Ratios of Significant Adverse Events",
       x = "Adverse Event",
       y = "Hazard Ratio (95% CI)") +
  theme_minimal()
# Save the plot to a PDF file
pdf("Hazard Ratios of Significant Adverse Events.pdf", width = 8, height = 6)
print(p_significant)
dev.off()

# Kaplan-Meier plots for significant adverse events using OS dataset
pdf(file = "significant_ae_km_plots.pdf", width = 8, height = 6)

for (ae in significant_results$AE) {
  # Determine which patients experienced the adverse event
  ae_patients <- aedata %>%
    filter(AE == ae) %>%
    pull(PersonId) %>%
    unique()
  
  # Add the AE_presence column to osdata
  os_data <- osdata %>%
    mutate(AE_presence = ifelse(PersonId %in% ae_patients, 1, 0))
  
  # Fit the Kaplan-Meier survival model
  km_fit <- survfit(Surv(SurvTime, OScen) ~ AE_presence, data = os_data)
  
  # Create the Kaplan-Meier plot
  p <- ggsurvplot(km_fit,
                  data = os_data,
                  pval = TRUE,
                  conf.int = TRUE,
                  risk.table = TRUE,
                  ggtheme = theme_minimal(),
                  title = paste("Kaplan-Meier Survival Curve for", ae))
  
  # Print the plot
  print(p)
  
  # Add a new page for each plot
  grid.newpage()
}

dev.off()

# Additional analyses for placebo and vandetanib groups
results_existing <- data.frame(AE = character(),
                               HR_All = numeric(),
                               CI_Lower_All = numeric(),
                               CI_Upper_All = numeric(),
                               P_Value_All = numeric(),
                               HR_Placebo = numeric(),
                               CI_Lower_Placebo = numeric(),
                               CI_Upper_Placebo = numeric(),
                               P_Value_Placebo = numeric(),
                               HR_Vandetanib = numeric(),
                               CI_Lower_Vandetanib = numeric(),
                               CI_Upper_Vandetanib = numeric(),
                               P_Value_Vandetanib = numeric(),
                               stringsAsFactors = FALSE)

for (cand.ae in unique(aedata$AE)) {
  cand.id <- which(aedata$AE == cand.ae)
  cand.date <- aggregate(list("aeDate" = aedata$Date[cand.id]), by = list("PersonId" = aedata$PersonId[cand.id]), FUN = min)
  cand.date$AE <- cand.ae
  
  data <- merge(osdata, cand.date, by = "PersonId", all.x = TRUE)
  
  data$OSdt <- as.Date(data$OSdt, format = "%d/%m/%Y")
  data$REGDate <- as.Date(data$REGDate, format = "%d/%m/%Y")
  data$aeDate <- as.Date(data$aeDate)
  
  anal_data <- data.frame()
  unique_patients <- unique(data$PersonId)
  
  for (pid in unique_patients) {
    p_dat <- subset(data, PersonId == pid)
    p_dat$start <- 0
    p_dat$end <- p_dat$SurvTime
    p_dat$ae <- 0
    p_dat$status <- p_dat$OScen
    newDat <- p_dat
    
    if (!is.na(p_dat$aeDate[1])) {
      newDat <- rbind(newDat, p_dat)
      newDat$end[1] <- as.numeric(difftime(p_dat$aeDate[1], p_dat$REGDate[1], units = "weeks")) / 4.35
      newDat$start[2] <- as.numeric(difftime(p_dat$aeDate[1], p_dat$REGDate[1], units = "weeks")) / 4.35
      newDat$ae[2] <- 1  # Corrected to set the AE occurrence properly
    }
    
    anal_data <- rbind(anal_data, newDat)
  }
  
  # Debugging step: Check for negative or NA end times
  print(head(anal_data))  # Check the first few rows of the data
  invalid_times <- anal_data[which(anal_data$end <= anal_data$start | is.na(anal_data$end) | is.na(anal_data$start)), ]
  if (nrow(invalid_times) > 0) {
    print("Invalid time intervals found:")
    print(invalid_times)
  }
  
  anal_data <- unique(anal_data)
  
  # Only proceed if there are valid times
  if (nrow(anal_data) > 0 && all(anal_data$end > anal_data$start, na.rm = TRUE)) {
    all.cph <- coxph(Surv(start, end, status) ~ ae, data = anal_data)
    all.cph.sum <- summary(all.cph)
    hr.all <- round(all.cph.sum$coefficients[1, 2], 2)
    ci.l.all <- round(all.cph.sum$conf.int[1, 3], 2)
    ci.u.all <- round(all.cph.sum$conf.int[1, 4], 2)
    p.all <- round(all.cph.sum$coefficients[1, 5], 2)
    
    placebo.cph <- coxph(Surv(start, end, status) ~ ae, data = anal_data, subset = (Arm == "Placebo"))
    placebo.cph.sum <- summary(placebo.cph)
    hr.placebo <- round(placebo.cph.sum$coefficients[1, 2], 2)
    ci.l.placebo <- round(placebo.cph.sum$conf.int[1, 3], 2)
    ci.u.placebo <- round(placebo.cph.sum$conf.int[1, 4], 2)
    p.placebo <- round(placebo.cph.sum$coefficients[1, 5], 2)
    
    van.cph <- coxph(Surv(start, end, status) ~ ae, data = anal_data, subset = (Arm == "Vandetanib"))
    van.cph.sum <- summary(van.cph)
    hr.van <- round(van.cph.sum$coefficients[1, 2], 2)
    ci.l.van <- round(van.cph.sum$conf.int[1, 3], 2)
    ci.u.van <- round(van.cph.sum$conf.int[1, 4], 2)
    p.van <- round(van.cph.sum$coefficients[1, 5], 2)
    
    temp <- data.frame(AE = cand.ae,
                       HR_All = hr.all,
                       CI_Lower_All = ci.l.all,
                       CI_Upper_All = ci.u.all,
                       P_Value_All = p.all,
                       HR_Placebo = hr.placebo,
                       CI_Lower_Placebo = ci.l.placebo,
                       CI_Upper_Placebo = ci.u.placebo,
                       P_Value_Placebo = p.placebo,
                       HR_Vandetanib = hr.van,
                       CI_Lower_Vandetanib = ci.l.van,
                       CI_Upper_Vandetanib = ci.u.van,
                       P_Value_Vandetanib = p.van,
                       stringsAsFactors = FALSE)
    
    results_existing <- rbind(results_existing, temp)
  } else {
    print(paste("Skipping analysis for AE:", cand.ae, "due to invalid times."))
  }
}

# Display the results
print(results_existing)

# Save results to CSV files
write.csv(results_existing, "results_existing.csv", row.names = FALSE)

# Plot the results
plot_all <- ggplot(results_existing, aes(x = reorder(AE, -HR_All), y = HR_All)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_Lower_All, ymax = CI_Upper_All), width = 0.2) +
  coord_flip() +
  labs(title = "Hazard Ratios for All Adverse Events",
       x = "Adverse Event",
       y = "Hazard Ratio (95% CI)") +
  theme_minimal()
# Save the plot to a PDF file
pdf("All Adverse Events.pdf", width = 8, height = 6)
print(plot_all)
dev.off()

# Filter significant results (p < 0.05) and plot
significant_existing <- results_existing %>%
  filter(P_Value_All < 0.05)

# Save results to CSV files
write.csv(significant_existing, "significant_results_existing.csv", row.names = FALSE)

plot_significant <- ggplot(significant_existing, aes(x = reorder(AE, -HR_All), y = HR_All)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_Lower_All, ymax = CI_Upper_All), width = 0.2) +
  coord_flip() +
  labs(title = "Hazard Ratios for Significant Adverse Events",
       x = "Adverse Event",
       y = "Hazard Ratio (95% CI)") +
  theme_minimal()
# Save the plot to a PDF file
pdf("Significant Adverse Events.pdf", width = 8, height = 6)
print(plot_significant)
dev.off()

#-------------------------------------------------------------------------

# Categorize AEs based on their severity using the Grade column
aedata <- aedata %>%
  mutate(AE_severity = AEgrade)

# Define the Extended Cox proportional hazards model
extended_cox_model <- coxph(
  Surv(SurvTime, OScen) ~ AE + AE_severity ,
  data = aedata
)

# Summarize the model
summary(extended_cox_model)

# Visualize the results with Kaplan-Meier plots for significant AEs
pdf(file = "extended_cox_km_plots.pdf", width = 8, height = 6)
for (ae in unique(aedata$AE)) {
  os_data <- aedata %>%
    mutate(AE_presence = ifelse(AE == ae, 1, 0))
  
  km_fit <- survfit(Surv(SurvTime, OScen) ~ AE_presence, data = os_data)
  
  p <- ggsurvplot(km_fit,
                  data = os_data,
                  pval = TRUE,
                  conf.int = TRUE,
                  risk.table = TRUE,
                  ggtheme = theme_minimal(),
                  title = paste("Kaplan-Meier Survival Curve for", ae))
  
  print(p)
  
  # # Add a new page for each plot
  # grid.newpage()
}
dev.off()

# Evaluate potential confounders and alternative models
# Example: adding age, gender, and other demecog values as potential confounders
# Assuming demecog variables include Age, Gender, and Comorbidity1, Comorbidity2, etc.
extended_cox_model_conf <- coxph(
  Surv(SurvTime, OScen) ~ AE + AE_severity + PatAge + REGSex + DEMECOG,
  data = aedata
)

# Summarize the model with confounders
summary(extended_cox_model_conf)

# Discuss results
# Extracting significant results
significant_results <- summary(extended_cox_model)$coefficients
significant_results_conf <- summary(extended_cox_model_conf)$coefficients

# Print significant results
print(significant_results)
print(significant_results_conf)
