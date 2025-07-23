# Script for cleaning and merging Zarahs soil aggregate data

# Load necessary libraries
library(dplyr)
library(readr)
library(readxl)
library(writexl)

# load data
spring24 <- read_xlsx("C:/Users/mm68tulo/Desktop/Jena Chapter/Data Analysis/1_Data/Aggregates/Results_Spring2024_WSA_all.xlsx")
august24 <- read_xlsx("C:/Users/mm68tulo/Desktop/Jena Chapter/Data Analysis/1_Data/Aggregates/WSA_August_clean.xlsx")
november24 <- read_xlsx("C:/Users/mm68tulo/Desktop/Jena Chapter/Data Analysis/1_Data/Aggregates/WSA_November_clean.xlsx")
winter25 <- read_xlsx("C:/Users/mm68tulo/Desktop/Jena Chapter/Data Analysis/1_Data/Aggregates/WSA_Winter_all.xlsx")

# Combine datasets
mean_aggregates <- spring24 %>% 
  group_by(Plotcode) %>%
  summarise(
    mean_spring = mean(agg_overall, na.rm = TRUE),
    .groups = 'drop'
  )

mean_winter <- winter25 %>% 
  group_by(Plotcode) %>%
  summarise(mean_winter = mean(agg_overall, na.rm = TRUE), .groups = 'drop')

mean_aggregates$mean_august <- august24$agg_overall
mean_aggregates$mean_november <- november24$agg_overall
mean_aggregates$mean_winter <- mean_winter$mean_winter

# calculate mean across seasons
mean_aggregates <- mean_aggregates %>%
  mutate(
    mean_all = rowMeans(select(., starts_with("mean_")), na.rm = TRUE)
  )

# Save the final dataset
write_xlsx(mean_aggregates, "C:/Users/mm68tulo/Desktop/Jena Chapter/Data Analysis/1_Data/Aggregates/mean_aggregates.xlsx")

