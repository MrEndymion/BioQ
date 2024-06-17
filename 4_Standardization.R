#========================================================================
#
# STANDARDIZATION OF INDICIE
#
# by Wong Xin Huei, DHI June 2023
#
#========================================================================

#=============================================
# Reading data from previous scripts
#=============================================

Tax_data <- Taxonomic_indicies
Fun_data <- Functional_trait
Inc_data <- Interactiv_indicies

#Removing S and N from Tax
# Tax_data <- subset(Tax_data, select = -c(S, NIS,n.rare , a.rare , N)) #Remove NIS n.rare etc as well?
Tax_data <- subset(Tax_data, select = -c(S, N)) #Remove NIS n.rare etc as well?

Indicie_all <- merge(Tax_data, Fun_data, by = c("Sample"))
Indicie_all <- merge(Indicie_all, Inc_data, by = c("Sample"))

Reference_site <- c("Ref.1", "Ref.2", "Ref", "Ref.") 

#=============================================
# Calculating standard site if multiple refs
#=============================================

Reference_globe <- subset(Indicie_all, Sample %in% Reference_site)

Reference_globe$Sample <- "Ref_globe"
if(nrow(Reference_globe) > 1) {Reference_globe <- ddply(Reference_globe, .(Sample), colwise(mean, na.rm=TRUE))}

#=============================================
# Standardizing and calculating
#=============================================

Score_year <- data.frame()

Score_each<-Indicie_all
var_cal <- Indicie_all[,-1]
var_overall <- ddply(var_cal, .(), colwise(sd, na.rm=TRUE)) # calculation of the column-wise standard deviation.

#Normalizing scores
Normalized_station <- data.frame()
for (station in unique(Score_each$Sample)){
  Impact <- subset(Score_each, Sample == station)
  
  Results <- Impact[,"Sample"]
  for (s in names(Impact[,-1])){
    Indicie <- Impact[,s]
    Ref_val <- Reference_globe[,s]
    Score <- ifelse(var_overall[,s] > 0,(Indicie- Ref_val)/var_overall[,s], (Indicie- Ref_val)/Ref_val)
    
    Results <- cbind(Results, Score) 
    
  }
  colnames(Results) <- c("Sample", names(Impact[,-1]))
  Normalized_station <- rbind(Normalized_station, Results)
  
}

Score_year <- rbind(Score_year, Normalized_station)

#Change all the numeric varibales back to numeric

Score_year[,-1] <- lapply(Score_year[, -1], as.numeric)



#=============================================
# Adding BioQ calculations
#=============================================

Tax_indicies <- names(subset(Tax_data, select = -c(Sample)))
Fun_indicies <- names(subset(Fun_data, select = -c(Sample)))
Inc_indicies <- names(subset(Inc_data, select = -c(Sample)))

BioQ_scores <- data.frame()


Tax_Scores <-  Score_year[,c(Tax_indicies)]
Tax_Scores[is.na(Tax_Scores)] <- 0
Fun_Scores <-  Score_year[,c(Fun_indicies)]
Fun_Scores[is.na(Fun_Scores)] <- 0
Inc_Scores <-  Score_year[,c(Inc_indicies)]
Inc_Scores[is.na(Inc_Scores)] <- 0

Tax_Score <- rowMeans(Tax_Scores)
Fun_Score <- rowMeans(Fun_Scores)
Inc_Score <- rowMeans(Inc_Scores)

Score_year$Tax_mean <- Tax_Score
Score_year$Fun_mean <- Fun_Score
Score_year$Inc_mean <- Inc_Score


Score_year$BioQ <-  rowMeans(cbind(Tax_Score, Fun_Score, Inc_Score))

BioQ_scores<-Score_year


# Gettign S and N for the visualisation:

Abundance <- Taxonomic_indicies
Abundance <- subset(Abundance, select = c(Sample, S, N))

Abundance <- merge(Abundance, Phyla, by = c("Sample"))

Abundance$TotalN <- sum(Abundance$N)

BioQ_scores <- merge(BioQ_scores, Abundance, by = c("Sample"))

#write.csv(BioQ_scores, "./Indput-Output/Standardized_BioQ.csv")

print("4_BioQ - Done")
# print("BioQ_scores calculated")



# Perform operations on the data
summary_data <- BioQ_scores

# Return the summary data
return(BioQ_scores)
