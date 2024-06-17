#========================================================================
#
# FUNCTIONAL INDICIE CALCULATIONS
#
# by Lars O. Mortensen, DHI Feb 2023
#
#========================================================================


Platform_year_genus <- ddply(Species_Summary, .(Sample, Species), summarise, 
                             Abundance = sum(Abundance, na.rm = TRUE))
Species_to_genus <- unlist(lapply(strsplit(Platform_year_genus$Species, " "), "[", 1))

Platform_year_genus$Species <- Species_to_genus

names(Trait_list)[names(Trait_list) == '_Order'] <- 'Order'

#Trait_simplification
Trait_list_Family <- subset(Trait_list, Genus == "")
Trait_list_Order <- subset(Trait_list, Genus == "" & Family == "")
Trait_list_Class <- subset(Trait_list, Genus == "" & Family == "" & Order == "")
Trait_list_Phylum <- subset(Trait_list, Genus == "" & Family == "" & Order == "" & Class == "")

#Subsertting data

Trait_genus <- subset(Platform_year_genus, Species %in% Trait_list$Genus)
Trait_genus$Genus <- Trait_genus$Species
Family_genus <- subset(Platform_year_genus, Species %in% Trait_list_Family$Family)
Family_genus$Family <- Family_genus$Species
Order_genus <- subset(Platform_year_genus, Species %in% Trait_list_Order$Order)
Order_genus$Order <- Order_genus$Species
Class_genus <- subset(Platform_year_genus, Species %in% Trait_list_Class$Class)
Class_genus$Class <- Class_genus$Species
Phylum_genus <- subset(Platform_year_genus, Species %in% Trait_list_Phylum$Phylum)
Phylum_genus$Phylum <- Phylum_genus$Species
Leftover <- subset(Platform_year_genus, !Species %in% c(Trait_list$Genus, Trait_list_Family$Family,
                                                        Trait_list_Order$Order, Trait_list_Class$Class, Trait_list_Phylum$Phylum))

Numbers_iden <- as.data.frame(cbind(length(unique(Trait_genus$Genus)), 
                                    length(unique(Family_genus$Family)), 
                                    length(unique(Order_genus$Order)), 
                                    length(unique(Class_genus$Class)), 
                                    length(unique(Phylum_genus$Phylum)), 
                                    length(unique(Leftover$Genus))))



names(Numbers_iden) <- c("Genus", "Familiy", "Order", "Class", "Phylum", "Unident")

Trait_genus <- merge(Trait_genus, Trait_list, by = "Genus", all.x = TRUE)
Family_genus <- merge(Family_genus, Trait_list_Family, by = "Family", all.x = TRUE)
Family_genus$Genus <- Family_genus$Species
Family_genus <- Family_genus[,names(Trait_genus)]
Order_genus <- merge(Order_genus, Trait_list_Order, by = "Order", all.x = TRUE)
Order_genus$Genus <- Order_genus$Species
Order_genus <- Order_genus[,names(Trait_genus)]
Class_genus <- merge(Class_genus, Trait_list_Class, by = "Class", all.x = TRUE)
Class_genus$Genus <- Class_genus$Species
Class_genus <- Class_genus[,names(Trait_genus)]
Phylum_genus <- merge(Phylum_genus, Trait_list_Phylum, by = "Phylum", all.x = TRUE)
Phylum_genus$Genus <- Phylum_genus$Species
Phylum_genus <- Phylum_genus[,names(Trait_genus)]

Trait_platform <- rbind(Trait_genus, Family_genus, Order_genus, Class_genus,Phylum_genus)

# save amount of phyla for later use in the output file
Phyla<-Trait_platform %>%
  group_by(Sample) %>%
  summarise(numberPhyla = n_distinct(Phylum))

Phyla$TotalPhyla <- n_distinct(Trait_platform$Phylum)
Phyla$TotalSpecies <- n_distinct(Trait_platform$Species)

Trait_platform <- Trait_platform[,c(-1,-5:-9)]

#######################################################################################################
# Trait analysis
#######################################################################################################

#=============================================
# Community weighted mean
#=============================================
Abundance_matrix <- Trait_platform[,c("Sample", "Species", "Abundance")]

trait_matrix <- ddply(Trait_platform[,-1], .(Species), colwise(mean))
row.names(trait_matrix) <- trait_matrix$Species
trait_matrix <- trait_matrix[,-1:-2]

Abundance_matrix_cast <- acast(Abundance_matrix, Sample ~ Species , value.var='Abundance', fun.aggregate=sum, margins=FALSE)
Abundance_matrix_cast <- data.matrix(Abundance_matrix_cast)

resCWM <- functcomp(trait_matrix, Abundance_matrix_cast, CWM.type = "all")

resCWM$Sample <- row.names(resCWM)

CWM_frame_average <- rowMeans(resCWM[,-ncol(resCWM)], na.rm=TRUE)

CWM_cal <- as.data.frame(CWM_frame_average)
CWM_cal$Station <- row.names(CWM_cal)
names(CWM_cal) <- c("CWM", "Sample")

#=============================================
# making distance matrix https://cran.r-project.org/web/packages/fundiversity/vignettes/fundiversity.html
#=============================================

sr <- pcoa(gowdis(trait_matrix[, 1:6]) / max(gowdis(trait_matrix[, 1:6])))$vectors[,1]
Body <- pcoa(gowdis(trait_matrix[, 7:12]) / max(gowdis(trait_matrix[, 7:12])))$vectors[,1]
Size <- pcoa(gowdis(trait_matrix[, 13:16]) / max(gowdis(trait_matrix[, 13:16])))$vectors[,1]
Reproduction <- pcoa(gowdis(trait_matrix[, 17:20]) / max(gowdis(trait_matrix[, 17:20])))$vectors[,1]
zone <- pcoa(gowdis(trait_matrix[, 21:23]) / max(gowdis(trait_matrix[, 21:23])))$vectors[,1]
habit <- pcoa(gowdis(trait_matrix[, 24:29]) / max(gowdis(trait_matrix[, 24:29])))$vectors[,1]
layer <- pcoa(gowdis(trait_matrix[, 30:33]) / max(gowdis(trait_matrix[, 30:33])))$vectors[,1]
Trophic <- pcoa(gowdis(trait_matrix[, 34:39]) / max(gowdis(trait_matrix[, 34:39])))$vectors[,1]
Mobility <- pcoa(gowdis(trait_matrix[, 40:43]) / max(gowdis(trait_matrix[, 40:43])))$vectors[,1]
Turbation <- pcoa(gowdis(trait_matrix[, 44:48]) / max(gowdis(trait_matrix[, 44:48])))$vectors[,1]

all.dist <- cbind(data.frame(sr),data.frame(Body), data.frame(Size), data.frame(Reproduction),
                  data.frame(zone), data.frame(habit), data.frame(layer), data.frame(Trophic),
                  data.frame(Mobility), data.frame(Turbation))

#=============================================
# Testing fundiversity https://cran.r-project.org/web/packages/fundiversity/vignettes/fundiversity.html
#=============================================

#library(fundiversity)
library(FD)

resFD.alltraits <- suppressWarnings(dbFD(all.dist, Abundance_matrix_cast,
                                         calc.CWM = FALSE, stand.FRic = TRUE,
                                         messages = F, m = "min", calc.FRic = TRUE, calc.FDiv = FALSE))

# #Scores for FRic
FRic_frame <- data.frame(resFD.alltraits$FRic)
FRic_frame$Sample <- row.names(FRic_frame)

# Calculatin the rest of the functional indicies

FDis <- fd_fdis(all.dist, Abundance_matrix_cast)
FEve <- fd_feve(all.dist, Abundance_matrix_cast)
RaoQ <- fd_raoq(all.dist, Abundance_matrix_cast)

FDiv <- data.frame()
for(i in unique(Trait_platform$Sample)){
  data <- subset(Trait_platform, Sample == i)
  traits <- all.dist[rownames(all.dist) %in% unique(data$Species), ]  
  
  Test_matrix <- data[,c("Sample", "Species", "Abundance")]
  Test_matrix_cast <- acast(Test_matrix, Sample ~ Species , value.var='Abundance', fun.aggregate=sum, margins=FALSE)
  Test_matrix_cast <- data.matrix(Test_matrix_cast)
  
  temp <- c()
  for(t in 1:10){
    Result <- fd_fdiv(traits[,t, drop = FALSE], Test_matrix_cast)
    temp <- c(temp, Result[,2])
  }
  
  Result <- cbind(i, mean(temp))
  
  FDiv <- rbind(FDiv, Result)
}

FDiv$i <- as.numeric(FDiv$i)
FDiv$V2 <- as.numeric(FDiv$V2)

Functional_trait <- cbind(FDis, FDiv[,2], FEve[,2], FRic_frame[,1], RaoQ[,2])
names(Functional_trait) <- c("Sample", "FDis", "FDiv", "FEve", "FRic", "RaoQ")

Functional_trait <- merge(Functional_trait, CWM_cal, by = "Sample", all = TRUE)


print("2_BioQ - Done")