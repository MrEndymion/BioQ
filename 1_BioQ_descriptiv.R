#========================================================================
#
# TAXONOMIC INDICI CALCULATIONS
#
# by Lars O. Mortensen, DHI June 2022
#
#========================================================================

Taxonomic_indicies<-c()

#=============================
# Number of Individuals and Species
#=============================

#N
Species.N <- ddply(Species_Summary, .(Sample), summarise, N = sum(Abundance))

Species.N <- Species_Summary %>% 
  group_by(Sample) %>% 
  summarize(N = sum(Abundance))

# S
Species.S <- Species_Summary %>% 
  group_by(Sample, Species) %>% 
  summarize(S = sum(Abundance))

Species.S <- Species.S %>% 
  group_by(Sample) %>% 
  summarize(S = sum(S > 0, na.rm = TRUE))

#=============================================
# Species richness (d, Margaleff) d = (S - 1) / ln N
#=============================================

# S = number of species, N = total number of individuals in the sample
Richness <- merge(Species.S, Species.N, by = c("Sample"), all = TRUE )
Richness$d <- (Richness$S-1) / log(Richness$N)

#=============================================
# Diversity (Simpson index (1-lambda))
#=============================================

#Making Matrix
Div_matrix <- acast(Species_Summary, Sample ~ Species , value.var='Abundance', fun.aggregate=sum, margins=FALSE)

#Calculating simpson diveristy
Simp <- vegan::diversity(Div_matrix, "simpson")
Simp <- data.frame(Simp)
Simp$Sample <- rownames(Simp)

#=============================================
# Diversity (Shannon)
#=============================================

#Calculating shannon diveristy
Shan <- vegan::diversity(Div_matrix, "shannon")
Shan <- data.frame(Shan)
Shan$Sample <- rownames(Shan)

#=============================================
# Aggregating last diversity table
#=============================================

Diversity_table <- merge(Richness, Simp, by = "Sample")
Diversity_table <- merge(Diversity_table, Shan, by = "Sample")

#=============================================
# Calculating Pielous index
#=============================================

Diversity_table$Pielou <- Diversity_table$Shan/log(Diversity_table$S)

#=============================================
# Non-indigenous species
## @knitr NIS_overivew
#=============================================

# #Searching EASIN
# Species_names <- unique(Species_Summary$Species)
# 
# Species_names <- unique(Species_Summary$Species)
# Species_names_EOL <- Species_names
# Species_names <- gsub(" ", "%20", Species_names)
# 
# Species_status <- data.frame()
# cal <- 0
# for(i in Species_names){
#   AphiaID <- i
#   url <- sprintf("https://easin.jrc.ec.europa.eu/apixg/catxg/term/%s", AphiaID);
# 
#   #Get the actual data from the URL
#   classificationTree <- fromJSON(url)
#   
#   if(length(classificationTree) == 1){
#     Name <- i
#     IsEUConcern <- FALSE
#     needed <- as.data.frame(cbind(Name, IsEUConcern))
#     needed$Name <- gsub("%20", " ", needed$Name)
#   }
#   
#   if(length(classificationTree) > 1){
#     needed <- classificationTree[,c("Name", "IsEUConcern")]
#     needed <- needed[grepl(i, needed$Name, fixed = TRUE),]
#     if(nrow(needed) > 0){needed$Name <- i}
#     if(nrow(needed) == 0){
#       needed <- classificationTree[,c("Name", "IsEUConcern")]
#       needed <- needed[1,]
#       needed$Name <- i
#       needed$IsEUConcern <- FALSE
#       }
#     
#   }
#   
#   needed$NIS <- ifelse(needed$IsEUConcern == "FALSE", 0, 1)
#   needed <- ddply(needed, .(Name), summarise, NIS = sum(NIS))
#   needed$Name <- gsub("%20", " ", needed$Name)
#   
#   Species_status <- rbind(Species_status, needed)
#   
#   cal <- cal+1
#   print(cal/length(Species_names)*100)
# }
# 
# names(Species_status) <- c("Species", "NIS")
# Alien <- merge(Species_Summary, Species_status, by = "Species", all = TRUE)
#  
# Alien$NIS <- Alien$NIS * Alien$Abundance
# Alien_sub <- ddply(Alien, .(Sample), summarise, NIS = sum(NIS))
# 
# Diversity_table <- merge(Diversity_table, Alien_sub, by = c("Sample"), all.x = TRUE)
# 
# #=============================================
# # Estimating rarity
# #=============================================
# 
# data_species <- Species_Summary
# 
# Rarity <- data.frame()
# cal <- 0
# for(r in unique(Species_Summary$Species)){
#   data <- rl_history(r, key = "3d359a007f87594b24edd35a3d8a53bfb7a7ca81657eb97ea8f4ee2d958d82dd")
#   status <- data$result$category
#   status <- ifelse(length(status)> 0, status, "NON")
#   result <- cbind(r, status)
#   Rarity <- rbind(Rarity, result)
#   
#   cal <- cal+1
#   print(cal/length(unique(Species_Summary$Species))*100)
# }
# 
# names(Rarity) <- c("Species", "status")
# 
# data_species <- merge(data_species, Rarity, by = "Species")
# data_species$status.factor <- ifelse(data_species$status == "NON", 0, 1)
# data_species$status.abund <- ifelse(data_species$status == "NON", 0,data_species$Abundance)
# 
# Rarity_species <- ddply(data_species, .(Sample), summarise, n.rare = sum(status.factor), a.rare = sum(status.abund))
# Diversity_table <- merge(Diversity_table, Rarity_species, by = "Sample", all = TRUE)

Taxonomic_indicies <- base::rbind(Taxonomic_indicies, Diversity_table)


print("1_BioQ - Done")
