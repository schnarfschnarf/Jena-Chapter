rm(list = ls()) #clear everything

#install.packages("readxl")
#install.packages("tidyr")
#install.packages("tidyverse")
#install.packages("plm")
#install.packages("ggplot2")
#install.packages("car")
#install.packages("caret")
#install.packages("dplyr")
#install.packages("zoo")
#install.packages("stargazer")
#install.packages("systemfit")
#install.packages("robustsur")
#install.packages("estimatr")
#install.packages("broom")
#install.packages("gridExtra")
#install.packages("MASS")
#install.packages("repmod")
#install.packages("plotly")

library(readxl)
library(tidyr)
library(tidyverse)
library(plm)
library(ggplot2)
library(car)
library(caret)
library(dplyr)
library(zoo)
library(stargazer)
library(systemfit)
library(robustsur)
library(estimatr)
library(broom)
library(gridExtra)
library(MASS)
library(repmod)
library(plotly)

#Set Working Directory
getwd()
setwd("C:/Users/idivsh33caqo/Documents/00_PhDPES/InsuranceGrass/PES/Empirischer Teil/JenaExperiment/R")

#########################################
###DATA IMPORT ##########################
#########################################
finaldata <- read_xlsx("Data5.xlsx") #read data
finaldata <- finaldata[complete.cases(finaldata),]
finaldata <- subset(finaldata, finaldata$plot != "B1A09")
finaldata <- subset(finaldata, finaldata$plot != "B4A03")
#finaldata <- subset(finaldata, finaldata$year != "2002")
#finaldata <- subset(finaldata, finaldata$year != "2004")
#finaldata <- subset(finaldata, finaldata$year != "2006")
#finaldata <- subset(finaldata, finaldata$year != "2008")
#finaldata <- subset(finaldata, finaldata$year != "2011")
#finaldata <- subset(finaldata, finaldata$year != "2014")
#finaldata <- subset(finaldata, finaldata$year != "2017")
finaldata$sowndiv <- as.numeric(finaldata$sowndiv) 
finaldata$func.group <- as.numeric(finaldata$func.group)
finaldata$numgrass <- as.numeric(finaldata$numgrass)
finaldata$numsherb <- as.numeric(finaldata$numsherb)
finaldata$numtherb <- as.numeric(finaldata$numtherb)
finaldata$numleg <- as.numeric(finaldata$numleg)


#Create Column with absolute changes of SOC on individual plots 
finaldata <- arrange(finaldata, plot)

finaldata <- finaldata %>%
  group_by(plot) %>%
  mutate(growth_SOC = `Corg_Soil (t/kt)` - dplyr::lag(`Corg_Soil (t/kt)` , order_by = plot))

finaldata <- arrange(finaldata, year)
ungroup(finaldata)
finaldata <- finaldata[complete.cases(finaldata),]


#Normalizing SOC_growth by calculating the average yearly growth rate (because time steps are irregular)

finaldata$SOC_growth_a <- ifelse(finaldata$year %in% c(2004, 2006, 2008), finaldata$growth_SOC / 2,
                                 ifelse(finaldata$year %in% c(2011, 2014, 2017), finaldata$growth_SOC / 3, finaldata$growth_SOC))


#####################
### DATA CREATION ###
#####################
##Getting mean, variance and covariance for each plot over time 
## SOC in g/kg
#Mean
SOC_Mean <- ave(finaldata$SOC_growth_a, finaldata$plot, FUN = mean) #calculates the mean of C on each plot over time and assigns this value to each time step for the respective plot 
finaldata$SOC_Mean <- SOC_Mean  #add vector as column to data

#Variance
SOC_Var <- ave(finaldata$SOC_growth_a, finaldata$plot, FUN = var)  #calculates the variance of C on each plot over time and assigns this value to each time step for the respective plot 
finaldata$SOC_Var <- SOC_Var #add vector as column to data

#Standard Deviation
SOC_Sd <- ave(finaldata$SOC_growth_a, finaldata$plot, FUN = sd)  #calculates the variance of C on each plot over time and assigns this value to each time step for the respective plot 
finaldata$SOC_Sd<-SOC_Sd #add vector as column to data


## Biomass in g/m2
#Mean
Yield_Mean <- ave(finaldata$`Sown plant biom [t/km2]`, finaldata$plot, FUN = mean) #calculates the mean of yield on each plot over time and assigns this value to each time step for the respective plot 
finaldata$Yield_Mean<- Yield_Mean #add vector as column to data 

#Variance
Yield_Var <- ave(finaldata$`Sown plant biom [t/km2]`, finaldata$plot, FUN = var) #calculates the variance of yield on each plot over time and assigns this value to each time step for the respective plot 
finaldata$Yield_Var <- Yield_Var #add vector as column to data

#Standard Deviation 
Yield_Sd <- ave(finaldata$`Sown plant biom [t/km2]`, finaldata$plot, FUN = sd) #calculates the variance of yield on each plot over time and assigns this value to each time step for the respective plot 
finaldata$Yield_Sd <- Yield_Sd #add vector as column to data

##Covariance SOC and Yield: cov(x) = sd_y* sd_C
CovCY <- Yield_Sd*SOC_Sd
finaldata$Cov_CY <- CovCY


#Coefficient of Variation 
CoefficientVariation_Y <- Yield_Sd/Yield_Mean
finaldata$CoefficientVariation_Y <- CoefficientVariation_Y

CoefficientVariation_C <- SOC_Sd/SOC_Mean
finaldata$CoefficientVariation_C <- CoefficientVariation_C 


# Cutting the dataset to the relevant cross-sectional plot data
data80 <- finaldata[1:80, -c(6:11)]
##########################
#Creating Relevant Variables 
Yield <- data80$Yield_Mean
Yield_Var <- data80$Yield_Var
Yield_CV <- data80$CoefficientVariation_Y
SOC <- data80$SOC_Mean
SOC_Var <- data80$SOC_Var
SOC_CV <- data80$CoefficientVariation_C
Cov <- data80$Cov_CY
Sowndiv <- data80$sowndiv


###################
#Seemingly Unrelated Regression Means
eq1 <- Yield ~ log(Sowndiv)  #Yield 
eq2 <- SOC   ~ log(Sowndiv)  #SOC 
system <- list(eq1 = eq1, eq2 = eq2)
sur_fit <- surerob(system, data = data80) 
summary(sur_fit)


residuals <- sur_fit$residuals
residuals_yield <- residuals[1:80]
residuals_SOC <- residuals[81:160]
qqnorm(residuals_yield)
qqline(residuals_yield)
qqnorm(residuals_SOC)
qqline(residuals_SOC)

#Plot 
b0_Yield <- sur_fit$coefficients[1]
b1_Yield <- sur_fit$coefficients[2]
b0_SOC <- sur_fit$coefficients[3]
b1_SOC <- sur_fit$coefficients[4]


# Create a ggplot plot for MeanYield with logarithmic axes and abline
plot_yield <- ggplot(data80, aes(x = log(Sowndiv), y = Yield)) +
  geom_point() +
  labs(x = "log(Sown Diversity)", y = "Aboveground Biomass in t/km2") +
  ggtitle("Biomass and Sown Diversity") +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 6)
  )+
  geom_abline(intercept = b0_Yield, slope = b1_Yield, col = "green")


# Create a ggplot plot for MeanSOC with logarithmic axes and abline
plot_soc <- ggplot(data80, aes(x = log(Sowndiv), y = SOC)) +
  geom_point() +
  labs(x = "log(Sown Diversity)", y = "Change in SOC in t/kt") +
  ggtitle("SOC Change in Sown Diversity") +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 6)
  )+
  geom_abline(intercept = b0_SOC, slope = b1_SOC, col = "green")


#Seemingly Unrelated Regression Variances
eq3 <- Yield_Var ~ log(Sowndiv) #Variance Yield
eq4 <- SOC_Var ~ log(Sowndiv) #Variance SOC 
system2 <- list(eq3 = eq3, eq4 = eq4)
sur_fit2 <- surerob(system2, data = data80)
summary(sur_fit2)

#Plot 
b0_Yield_Var <- sur_fit2$coefficients[1]
b1_Yield_Var <- sur_fit2$coefficients[2]
b0_SOC_Var <- sur_fit2$coefficients[3]
b1_SOC_Var <- sur_fit2$coefficients[4]


# Create a ggplot plot for MeanYield with logarithmic axes and abline
plot_yield_var <- ggplot(data80, aes(x = log(Sowndiv), y = Yield_Var)) +
  geom_point() +
  labs(x = "log(Sown Diversity)", y = "Variance of Biomass in t/km2") +
  ggtitle("Variance of Biomass and Sown Diversity") +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 6)
  )+
  geom_abline(intercept = b0_Yield_Var, slope = b1_Yield_Var, col = "green")


# Create a ggplot plot for MeanSOC with logarithmic axes and abline
plot_soc_var <- ggplot(data80, aes(x = log(Sowndiv), y = SOC_Var)) +
  geom_point() +
  labs(x = "log(Sown Diversity)", y = "Variance of Change in SOC in t/kt") +
  ggtitle("Variance of SOC change and Sown Diversity") +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 6)
  )+
  geom_abline(intercept = b0_SOC_Var, slope = b1_SOC_Var, col = "green")


#Coefficients of Variation
eq5 <- Yield_CV ~ log(Sowndiv) #Coefficient of Variation Yield
eq6 <- SOC_CV ~ log(Sowndiv) #Coefficient of Variation SOC 
system3 <- list(eq5 = eq5, eq6 = eq6)
sur_fit3 <- surerob(system3, data = data80)
summary(sur_fit3)

#Plot 
b0_Yield_CV <- sur_fit3$coefficients[1]
b1_Yield_CV <- sur_fit3$coefficients[2]
b0_SOC_CV <- sur_fit3$coefficients[3]
b1_SOC_CV <- sur_fit3$coefficients[4]


# Create a ggplot plot for MeanYield with logarithmic axes and abline
plot_yield_CV <- ggplot(data80, aes(x = log(Sowndiv), y = Yield_CV)) +
  geom_point() +
  labs(x = "log(Sown Diversity)", y = "Coefficient of Variation Biomass in t/km2") +
  ggtitle("Coefficient of Variation Biomass and Sown Diversity") +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 6)
  )+
  geom_abline(intercept = b0_Yield_CV, slope = b1_Yield_CV, col = "green")


# Create a ggplot plot for MeanSOC with logarithmic axes and abline
plot_soc_CV <- ggplot(data80, aes(x = log(Sowndiv), y = SOC_CV)) +
  geom_point() +
  labs(x = "log(Sown Diversity)", y = "Coefficient of Variation SOC change in t/kt") +
  ggtitle("Coefficient of Variation SOC change and Sown Diversity") +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 6)
  )+
  geom_abline(intercept = b0_SOC_CV, slope = b1_SOC_CV, col = "green")

# Convert to an interactive plotly plot
interactive_plot_soc <- ggplotly(plot_soc_CV)
interactive_plot_soc


#Covariance
model4 <- rlm(Cov ~ log(Sowndiv))
summary(model4)
rob.pvals(model4)


#Plot 
b0_Cov <- model4$coefficients[1]
b1_Cov <- model4$coefficients[2]

# Create a ggplot plot for MeanSOC with logarithmic axes and abline
plot_Cov <- ggplot(data80, aes(x = log(Sowndiv), y = Cov)) +
  geom_point() +
  labs(x = "log(Sown Diversity)", y = "Covariance") +
  ggtitle("Covariance of Biomass and SOC change and Sown Diversity") +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 6)
  )+
  geom_abline(intercept = b0_Cov, slope = b1_Cov, col = "green")


#########################################
###DATA IMPORT ##########################
#########################################
finaldata2 <- read_xlsx("C:/Users/idivsh33caqo/Documents/00_PhDPES/InsuranceGrass/PES/Empirischer Teil/JenaExperiment/Kosten/seed_price_data.xlsx") #read data

x <- finaldata2$number_species 
y <- finaldata2$price_kg #EUR/kg oder TEUR/ton

model1 <- rlm(log(y) ~ x)
summary(model1)
rob.pvals(model1)

#Plot 
b0_Costs <- model1$coefficients[1]
b1_Costs <- model1$coefficients[2]

# Create a ggplot plot for MeanSOC with logarithmic axes and abline
plot_Costs <- ggplot(finaldata2, aes(x = x, y = log(y))) +
  geom_point() +
  ggtitle("Costs and Sown Diversity") +
  labs(x = "Sown Diversity)", y = "Costs of Seeds in TEUR/t") +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 6)) + 
  geom_abline(intercept = b0_Costs, slope = b1_Costs, col = "green")
print(plot_Costs)




grid.arrange(plot_yield, plot_yield_CV, plot_soc, plot_soc_CV, plot_Cov, plot_Costs, ncol = 3)

















