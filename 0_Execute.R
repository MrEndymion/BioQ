#========================================================================
#
# DATA SET TESTING
#
# by Lars O. Mortensen, DHI June 2023
#
#========================================================================

library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(rredlist)
library(vegan)
library(jsonlite)
library(FD)
library(igraph)
library(bipartite)
library(tidyr)
library(ggplot2)
library(terra)
library(raster)
library(fundiversity)
library(randomForest)

#=============================
# Reading data
#=============================
Species_Summary <- fread(file.path("./Indput-Output/EBM_BioQ_Indput.csv"))
Species_Summary <- Species_Summary[,c("Sample", "Species", "Abundance")]
Trait_list <- fread(file.path("./bioq-api/data/CEFAS_Trait_DataBase.csv"))

# #run the 4 scripts:
   source("./1_BioQ_descriptiv.R")
   source("./2_BioQ_traitbased.R")
   source("./3_BioQ_interactions.R")
   source("./4_Standardization.R")

#Outputs write to output file directory
write.csv(BioQ_scores, "./Indput-Output/EBM_BioQ_output.csv")

