install.packages("tableone")  
install.packages("data.table")

library(tableone)
library(data.table)


df <- fread("aortic_Height_merged.tsv")


vars <- c(
  "Gender_Legal_Sex", 
  "Age", 
  "Race1", 
  "SBP_Mean", 
  "DBP_Result", 
  "smoking_status",
  "Prevalent_Hypertension", 
  "Prevalent_Hyperlipidemia", 
  "Prevalent_CAD", 
  "Incident_CAD",
  "Prevalent_COPD", 
  "ldl", 
  "hdl", 
  "bmi", 
  "tc", 
  "tg",
  "Mean_Result"
)


cat_vars <- c(
  "Gender_Legal_Sex", "Race1", "smoking_status",
  "Prevalent_Hypertension", "Prevalent_Hyperlipidemia", 
  "Prevalent_CAD", "Incident_CAD", "Prevalent_COPD",
  "Incident_Hypertension", "Incident_Hyperlipidemia"
)


table1 <- CreateTableOne(vars = vars, data = df, factorVars = cat_vars)


print(table1, showAllLevels = TRUE)

