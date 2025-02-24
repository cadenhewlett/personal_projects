library(here)
library(tidyverse)
library(MASS)
set.seed(1928)
load(here("OSF_WFH.RData"))
data <- OSF_WFH

# self-reported productivity and burnout propensity, as well as their willingness
# to continue working from home. higher satisfaction with home
# environment factors significantly predicts increased productivity and 
# decreased burnout propensity

# GC - general control
# Both models include a set of
# carefully selected general controls (GC) that could otherwise confound the estimators.
# Specifically, these include demographic characteristics (gender and age), job characteristics
# (company size, job suitable for working from home, type of work, income, and working hours),
# and household characteristics (household size, children at home, partner at home, and pets
# y is the predicted value of either productivity or burnout propensity for each
# participant i. 

# combine pets into single column
data$pets <- (data$`Pet - Dog` == "Yes") | (data$`Pet - Cat` == "Yes")
dims <- data$`WFH office length` * data$`WFh office Width` 
# quantile(dims)

# ------------ WFH ENVIRONMENT ------------ #
wfh = data.frame(
  work_loc      = data$`Where Work?`,
  productivity  = data$`Single-item Scale Productivity Home`,
  gender        = data$geslacht,
  age           = data$leeftijd_cat,
  partner       = data$`Partner Home during Office Hours`,
  office_room   = data$`Home Office Room`,
  office_light  = data$`Home Office Light`,
  home_light    = data$`Lighting Home`,
  home_desk     = data$`Desk Home`,
  home_wifi     = data$`WiFi Home`,
  home_stress   = data$`Stress and Irritabilitu Home`,
  work_hours    = data$werkuren,  #as.factor(ifelse(data$werkuren == "36+ hours", "Full-time", "Part-time")),
  comp_size     = as.factor(ifelse(data$`Company size` == "50+", "Large", "Small"))
)

# recast integers to factors
for (i in 1:ncol(wfh)){
  if(class(wfh[, i]) == "integer"){
    wfh[, i] = as.factor(wfh[, i])
  }
}

# # refine results to exclusively home workers

wfh <- wfh[wfh$work_loc ==  "Exclusively Home" &
             !is.na(wfh$work_loc), -1]
# wfh <- wfh[wfh$work_loc ==  "Exclusively Home" &
#              !is.na(wfh$work_loc), -1]
# random indices
indices <- sample(1:nrow(wfh), size = nrow(wfh)*0.75)
# train and test split
wfh_train <- wfh[indices,  ]
wfh_test  <- wfh[-indices, ]

