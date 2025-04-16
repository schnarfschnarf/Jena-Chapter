#---------------------------------------------------------------------

# rBEXIS package to read data directly form the https://jexis.idiv.de/
# plot informaiton from the main experiment
# ID: 90

#======================================================================
#==================================================================================
rm(list=ls()) 
#==========================  Install related packages  ===============================
#install devtools and add to libary
#install.packages("devtools")
#library(devtools)

#install required packages
#install.packages("httr")
#install.packages("XML")
#install.packages("jsonlite")
#install_github("BEXIS2/rBExIS", subdir = "rBExIS")

#add packages to the libary
library(httr)
library(jsonlite)
library(XML)
library(rBExIS)

library(tidyverse)

#==================================================================================
# Here start functions implemented by the package
#==================================================================================

bexis.options() # info of the url and token number

bexis.options(base_url = "https://jexis.uni-jena.de")

# Except the website url, you also need your account ID to see which dataset you have the rights to download
# This Identitity number is called "token"
# How to get the token number?
# 1. log in the https://jexis.idiv.de/
# 2. on the right upside, clike your name
# 3. there you will see the "token"; click it and copy the long string to "token_num"
token_num<-"" # add your token id here
bexis.options(token = token_num )

#get list of all dataset ids
bexis_dataset_ids <- bexis.get.datasets()

#read the data of id = 90, which is the main plot information
dp <- bexis.get.dataset_by(id = 90)


for (i in 1:nrow(dp)){
  for (v in c("grass","sherb","therb","leg")){
    v2<-sprintf("%s%s","num",v)
    dp[[v]][i]<-ifelse(dp[[v2]][i]>0,1,0)
  }
}


splitSet <- function(com) 
  sort(strsplit(com,"|",fixed=TRUE)[[1]])

species_list<-"Aju.rep"
for (i in 1:nrow(dp)){
  species_list1<-as.character(splitSet(dp$composition[i]))
  species_list<-c(species_list,species_list1)
}
sp<-unique(species_list) #60 species, exclude the blank
for (i in 1:nrow(dp)){
  for (v in sp){
    dp[[v]][i]<-ifelse(grepl(v,dp$composition[i],perl=T),1,0)
  }
}

write.csv(dp,file="dp.csv")
