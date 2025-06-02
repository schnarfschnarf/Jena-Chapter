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
library(repmod)
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
#finaldata <- subset(finaldata, finaldata$sowndiv != "60")
#finaldata <- subset(finaldata, finaldata$plot != "B4A03")
finaldata$sowndiv <- as.numeric(finaldata$sowndiv) 
finaldata$func.group <- as.numeric(finaldata$func.group)
#finaldata$numgrass <- as.numeric(finaldata$numgrass)
#finaldata$numsherb <- as.numeric(finaldata$numsherb)
#finaldata$numtherb <- as.numeric(finaldata$numtherb)
#finaldata$numleg <- as.numeric(finaldata$numleg)

#Create Column with absolute changes of SOC on individual plots 
finaldata <- arrange(finaldata, plot)

finaldata <- finaldata %>%
  group_by(plot) %>%
  mutate(growth_SOC = `Corg_Soil (t/kt)` - dplyr::lag(`Corg_Soil (t/kt)` , order_by = plot))

finaldata <- arrange(finaldata, year)
ungroup(finaldata)
finaldata <- finaldata[complete.cases(finaldata),]

#Normalizing SOC_growth by calculating the average yearly growth rate (because time steps are irregular)

finaldata$growth_SOC <- ifelse(finaldata$year %in% c(2004, 2006, 2008), finaldata$growth_SOC / 2,
                               ifelse(finaldata$year %in% c(2011, 2014, 2017), finaldata$growth_SOC / 3, finaldata$growth_SOC))

#Creating Relevant Variables 
Yield <- finaldata$`Sown plant biom [t/km2]` 
SOC <- finaldata$growth_SOC 
Sowndiv <- finaldata$sowndiv
time <- finaldata$year


#Seemingly Unrelated Regression 
#Define Variables
eq1 <- Yield ~ log(Sowndiv)  #Yield 
eq2 <- SOC   ~ log(Sowndiv)  #SOC 
system <- list(eq1 = eq1, eq2 = eq2)
sur_fit <- surerob(system, data = finaldata) 
summary(sur_fit)

#Plot 
b0_Yield <- sur_fit$coefficients[1]
b1_Yield <- sur_fit$coefficients[2]
b0_SOC <- sur_fit$coefficients[3]
b1_SOC <- sur_fit$coefficients[4]

# Create a ggplot plot for MeanYield with logarithmic axes and abline
plot_yield <- ggplot(finaldata, aes(x = Sowndiv, y = Yield)) +
  geom_point() +
  scale_x_log10() + 
  labs(x = "Sown Plant Species Diversity", y = "Aboveground Biomass in t/km2") +
  ggtitle("Diversity and Yield") +
  geom_abline(intercept = b0_Yield, slope = b1_Yield, col = "green")+
  theme(
    plot.title = element_text(size = 6),
    axis.title = element_text(size = 6)
  )

# Create a ggplot plot for MeanSOC with logarithmic axes and abline
plot_soc <- ggplot(finaldata, aes(x = Sowndiv, y = SOC)) +
  geom_point() +
  scale_x_log10() +
  labs(x = "Sown Plant Species Diversity", y = "Change in SOC in t/kt") +
  ggtitle("Diversity and SOC change") +
  geom_abline(intercept = b0_SOC, slope = b1_SOC, col = "green")+
  theme(
    plot.title = element_text(size = 6),
    axis.title = element_text(size = 6)
  )

#Extract residuals
residuals <- residuals(sur_fit)
residuals_yield <- residuals[, 1]
Var_Yield <- residuals_yield^2
residuals_SOC <- residuals[, 2]
Var_SOC <- residuals_SOC^2
Cov <- residuals_yield * residuals_SOC 

finaldata$Var_Yield <- Var_Yield
finaldata$Var_SOC <- Var_SOC 
finaldata$Cov <- Cov



# Regress Variances and Covariance against Sowndiv 
LM_VAR_YIELD <- rlm(Var_Yield ~ log(Sowndiv)) 
summary(LM_VAR_YIELD)
rob.pvals(LM_VAR_YIELD)

#Plot 
b0_yieldvar <- LM_VAR_YIELD$coefficients[1]
b1_yieldvar <- LM_VAR_YIELD$coefficients[2]
# Create a ggplot plot for MeanSOC with logarithmic axes and abline
plot_yieldvar <- ggplot(finaldata, aes(x = Sowndiv, y = Var_Yield)) +
  geom_point() +
  scale_x_log10() + 
  labs(x = "Sown Plant Species Diversity", y = "Change in biomass in t/km2") +
  ggtitle("Diversity and Variance of Biomass") +
  geom_abline(intercept = b0_yieldvar, slope = b1_yieldvar, col = "green")+
  theme(
    plot.title = element_text(size = 6),
    axis.title = element_text(size = 6)
  )


LM_VAR_SOC <- rlm(Var_SOC ~ log(Sowndiv)) 
summary(LM_VAR_SOC) 
rob.pvals(LM_VAR_SOC)

#Plot 
b0_socvar <- LM_VAR_SOC$coefficients[1]
b1_socvar <- LM_VAR_SOC$coefficients[2]
plot_socvar <- ggplot(finaldata, aes(x = Sowndiv, y = Var_SOC)) +
  geom_point() +
  scale_x_log10() +
  labs(x = "Sown Plant Species Diversity", y = "Change in SOC in t/kt") +
  ggtitle("Diversity and Variance of SOC") +
  geom_abline(intercept = b0_socvar, slope = b1_socvar, col = "green")+
  theme(
    plot.title = element_text(size = 6),
    axis.title = element_text(size = 6)
  )

LM_COV <- rlm(Cov ~ log(Sowndiv))
summary(LM_COV)
rob.pvals(LM_COV)

#Plot 
b0_cov <- LM_COV$coefficients[1]
b1_cov <- LM_COV$coefficients[2]

plot_cov <- ggplot(finaldata, aes(x = Sowndiv, y = Cov)) +
  geom_point() +
  scale_x_log10() + 
  labs(x = "Sown Plant Species Diversity", y = "Covariance") +
  ggtitle("Diversity and Covariance") +
  geom_abline(intercept = b0_cov, slope = b1_cov, col = "green")+
  theme(
    plot.title = element_text(size = 6),
    axis.title = element_text(size = 6)
  )


#Coefficients of Variation
VarCoef_Yield <- sqrt(Var_Yield)/Yield
VarCoef_Yield[is.infinite(VarCoef_Yield)] <- 0
finaldata$VarCoef_Yield <- VarCoef_Yield
VarCoef_SOC <- sqrt(Var_SOC)/SOC
VarCoef_SOC[is.infinite(VarCoef_SOC)] <- 0
finaldata$VarCoef_SOC <- VarCoef_SOC

VarCoefYield <- rlm(VarCoef_Yield ~ log(Sowndiv))
summary(VarCoefYield)
rob.pvals(VarCoefYield)

#Plot 
b0_yieldCV <- VarCoefYield$coefficients[1]
b1_yieldCV <- VarCoefYield$coefficients[2]

plot_yieldCV <- ggplot(finaldata, aes(x = Sowndiv, y = VarCoef_Yield)) +
  geom_point() +
  scale_x_log10() + 
  labs(x = "Sown Plant Species Diversity", y = "Coefficient of Variation Biomass (t/km2)") +
  ggtitle("Diversity and Coefficient of Variation Biomass") +
  geom_abline(intercept = b0_yieldCV, slope = b1_yieldCV, col = "green")+
  theme(
    plot.title = element_text(size = 6),
    axis.title = element_text(size = 6)
  )


VarCoefSOC<- rlm(VarCoef_SOC ~ log(Sowndiv))
summary(VarCoefSOC)
rob.pvals(VarCoefSOC)
#Plot 
b0_socCV <- VarCoefSOC$coefficients[1]
b1_socCV <- VarCoefSOC$coefficients[2]
plot_socCV <- ggplot(finaldata, aes(x = Sowndiv, y = VarCoef_SOC)) +
  geom_point() +
  scale_x_log10() + 
  labs(x = "Sown Plant Species Diversity", y = "Coefficient of Variation SOC change (t/kt)") +
  ggtitle("Diversity and Coefficient of Variation SOC change") +
  geom_abline(intercept = b0_socCV, slope = b1_socCV, col = "green")+
  theme(
    plot.title = element_text(size = 6),
    axis.title = element_text(size = 6)
  )



grid.arrange(plot_yield,plot_yieldvar, plot_yieldCV, plot_soc, plot_socvar, plot_socCV, plot_cov,  ncol = 3)



































# PLOTS
Mean_Yield <- ggplot(data = finaldata, aes(x = Sowndiv, y = Yield)) +
  geom_point() +   # Add dots for the data points
  geom_smooth(method = "lm", se = TRUE) +  # Add the regression line
  labs(x = "Sowndiv", y = "Biomass (t/km2)") +  # Set axis labels
  ggtitle("Plant Species Diversity and Biomass")+  # Set the plot title
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 6)
  )
print(Mean_Yield)

Mean_SOC <- ggplot(data = finaldata, aes(x = Sowndiv, y = SOC)) +
  geom_point() +   # Add dots for the data points
  geom_smooth(method = "lm", se = TRUE) +  # Add the regression line
  labs(x = "Sowndiv", y = "SOC change  (t/kt)") +  # Set axis labels
  ggtitle("Plant Species Diversity and SOC change")+  # Set the plot title
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 6)
  )
print(Mean_SOC)

# Create a scatter plot with regression line
Variance_Yield <- ggplot(data = finaldata, aes(x = finaldata$sowndiv, y = finaldata$Var_Yield)) +
  geom_point() +   # Add dots for the data points
  geom_smooth(method = "lm", se = TRUE) +  # Add the regression line
  labs(x = "Sowndiv", y = "Variance of Biomass (t/km2)") +  # Set axis labels
  ggtitle("Variance of Biomass and Sowndiv")+  # Set the plot title
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 6)
  )
print(Variance_Yield)

# Create a scatter plot with regression line
Variance_SOC <- ggplot(data = finaldata, aes(x = finaldata$sowndiv, y = finaldata$Var_SOC)) +
  geom_point() +   # Add dots for the data points
  geom_smooth(method = "lm", se = TRUE) +  # Add the regression line
  #scale_x_log10() +  # Set x-axis to logarithmic scale
  #scale_y_log10() +  # Set y-axis to logarithmic scale
  labs(x = "Sowndiv", y = "Variance of SOC change (t/kt)") +  # Set axis labels
  ggtitle("Variance of SOC change and Sowndiv") + # Set the plot title
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 6)
  )
print(Variance_SOC)

# Create a scatter plot with regression line
Covariance <- ggplot(data = finaldata, aes(x = finaldata$sowndiv, y = finaldata$Cov)) +
  geom_point() +   # Add dots for the data points
  geom_smooth(method = "lm", se = TRUE) +  # Add the regression line
  labs(x = "Sowndiv", y = "Covariance of SOC change and Biomass") +  # Set axis labels
  ggtitle("Covariance of SOC change and Biomass and Sowndiv") + # Set the plot title
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 6)
  )
print(Covariance)


# Create a scatter plot with regression line
CVYield <- ggplot(data = finaldata, aes(x = finaldata$sowndiv, y = finaldata$VarCoef_Yield)) +
  geom_point() +   # Add dots for the data points
  geom_smooth(method = "lm", se = TRUE) +  # Add the regression line
  labs(x = "Sowndiv", y = "CV Biomass (t/km2)") +  # Set axis labels
  ggtitle("Coefficient of Variation of  Biomass and Sowndiv") +  # Set the plot title
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 6)
  )
print(Covariance)


# Create a scatter plot with regression line
CVSOC<- ggplot(data = finaldata, aes(x = finaldata$sowndiv, y = finaldata$VarCoef_SOC)) +
  geom_point() +   # Add dots for the data points
  geom_smooth(method = "lm", se = TRUE) +  # Add the regression line
  labs(x = "Sowndiv", y = "CV SOC (t/kt)") +  # Set axis labels
  ggtitle("Coefficient of Variation of  SOC and Sowndiv") +  # Set the plot title
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 6)
  )
print(Covariance)



#Grid with plots
grid.arrange(Mean_Yield,Variance_Yield, CVYield, Mean_SOC, Variance_SOC, CVSOC, Covariance,  ncol = 3)

