rm(list=ls())
getwd()
library(readxl)
library(dplyr)
library(openxlsx)
library(tidyr)

#import data with xlsx format
ALA <- read_excel("ALA_ComALL.xlsx")

#select data with high FDR and IsMasterProtein
ALA_select <- ALA[ALA$FDR == 'High' & ALA$Master == 'IsMasterProtein', ]

#select data
ALA_grouped <-  ALA_select[,c(colnames(ALA_select)[grep("(Grouped)",colnames(ALA_select))],"Accession")]
ALA_grouped_1 <- ALA_grouped[,1:6]
ALA_grouped_2 <- ALA_grouped[,"Accession"]
ALA_grouped_3 <- mutate(ALA_grouped_2, ALA_grouped_1)

#change column names
colnames(ALA_grouped_1)[1] ="ALA_0µM_24hr"
colnames(ALA_grouped_1)[2] ="ALA_75µM_24hr"
colnames(ALA_grouped_1)[3] ="ALA_150µM_24hr"
colnames(ALA_grouped_1)[4] ="ALA_0µM_48hr"
colnames(ALA_grouped_1)[5] ="ALA_75µM_48hr"
colnames(ALA_grouped_1)[6] ="ALA_150µM_48hr"



#export data with selected column with multiple sheets
ALA_126 <-  ALA_select[,c(colnames(ALA_select)[grep("126",colnames(ALA_select))],"Accession","Description")]
ALA_127N <- ALA_select[,c(colnames(ALA_select)[grep("127N",colnames(ALA_select))],"Accession","Description")]
ALA_127C <- ALA_select[,c(colnames(ALA_select)[grep("127C",colnames(ALA_select))],"Accession","Description")]
ALA_128N <- ALA_select[,c(colnames(ALA_select)[grep("128N",colnames(ALA_select))],"Accession","Description")]
ALA_128C <- ALA_select[,c(colnames(ALA_select)[grep("128C",colnames(ALA_select))],"Accession","Description")]
ALA_129N <- ALA_select[,c(colnames(ALA_select)[grep("129N",colnames(ALA_select))],"Accession","Description")]

require(openxlsx)
list_of_datasets <- list("ALA_select" = ALA_select, 
                         "ALA_0µM_24hr" = ALA_126,
                         "ALA_75µM_24hr" = ALA_127N,
                         "ALA_150µM_24hr" = ALA_127C,
                         "ALA_0µM_48hr" = ALA_128N,
                         "ALA_75µM_48hr" = ALA_128C,
                         "ALA_150µM_48hr" = ALA_129N
                         )
#write.xlsx(list_of_datasets, file = "ALA_filtered.xlsx")

# Load promor
BiocManager::install("limma")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pcaMethods")
library(promor)

#data not using non numeric value and drop na value
rawdata <- 
data.frame(ID=rawdata$ID,value=apply(rawdata[,-1],1,quantile,0.75))

# Filter out proteins with high levels of missing data in either condition/group
raw_filtered <- filterbygroup_na(ALA_grouped_3)

# Impute missing data and create an imp_df object.
imp_df <- impute_na(raw_filtered)

# Normalize data and create a norm_df object
norm_df <- normalize_data(imp_df)

# Perform differential expression analysis and create a fit_df object
fit_df <- find_dep(norm_df)

volcano_plot(fit_df, text_size = 5)


