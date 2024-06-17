#========================================================================
#
# INTERACTIONS INDICIE CALCULATIONS
#
# by Ole B. Brodnicke, DHI Aug 2023
#
#========================================================================


#===========================================================
# CO-OCCURENCE USING IGRAPH ACROSS STATIONS
#===========================================================

print(Species_Summary)
station_matrix <- Species_Summary

co_occurence <- data.frame()
result_df<-data.frame()
for (R in unique(Species_Summary$Sample)) {
  
  station_year_data <- station_matrix[station_matrix$Sample == R & station_matrix$Abundance !=0, ] #Removing 0's for the further analysis
  
  # Calculate dissimilarity matrix
  demo_matrix <- vegdist(station_year_data[, 3], method = "bray") #Maybe we can look at other distances
  
  # Convert dissimilarity matrix to a distance matrix
  distance_matrix <- as.dist(demo_matrix)  #The values stay the same here ? Not sure this is needed.
  
  # Create a network object
  network <- graph.adjacency(as.matrix(distance_matrix), weighted = TRUE, mode = "undirected", diag = FALSE)
  
  # Set node names
  V(network)$name <- station_year_data[, 2]
  
  # Calculate network density
  E.density <- edge_density(network, loops = T) 
  
  # Calculate diameter
  N.diam <- 1/diameter(network, directed = FALSE) #Because a smaller diameter is better we change the range so that a larger number is better (inverse)
  
  # Calculate average path length
  A.path <- 1/mean_distance(network, directed = FALSE) #Because a shorter path is better we change the range so that a larger number is better (inverse)
  
  
  
  # Append results to the dataframe
  result_df <- cbind(Sample=R, E.density, N.diam, A.path) # Adding the Sample back to the values
  
  co_occurence <- rbind(co_occurence, result_df)
  
}

# transform the avlues back to numerical variables
co_occurence[,2:4] <- lapply(co_occurence[, 2:4], as.numeric)

#===========================================================
# Trophic-guild interaction
#===========================================================

#Trophic traits from CEFAS


station_matrix$Genus <- unlist(lapply(strsplit(station_matrix$Species, " "), "[", 1))

station_matrix_genus <- ddply(station_matrix, .( Sample, Genus), summarise, 
                              Abundance = sum(Abundance, na.rm = TRUE), Abundance = sum(Abundance, na.rm=TRUE))

Platform_year_traits <- merge(station_matrix_genus, Trait_list, by = "Genus", all.x = TRUE)

station_matrix_genus_2<-Platform_year_traits[,c("Genus", "Sample", "Abundance","f_Suspension","f_Surface_deposit","f_Subsurface_deposit","f_Scavenger","f_Predator", "f_Parasite")]


# Create an entry per genus per trophic trait (f_)
df_long <- station_matrix_genus_2 %>%
  pivot_longer(
    cols = starts_with("f_"),  # Specify the columns to pivot
    names_to = "Trophic",      # Name of the new categorical variable
    values_to = "Trait_Value"     # Name of the value column
  ) %>%
  filter(Trait_Value > 0)

#Remove a column we dont need
df_long$Trait_Value<-NULL


combined_data <- data.frame()
# Loop through each station and year and make unique trait combinations
for (station in unique(df_long$Sample)) {
  station_year_data <- df_long[df_long$Sample == station, ]
  
  if (length(unique(station_year_data$Trophic)) > 1) {
    trait_combinations <- unique(t(combn(unique(station_year_data$Trophic), 2)))
    
    trait_combinations_df <- data.frame(Trait1 = trait_combinations[, 1],
                                        Trait2 = trait_combinations[, 2])
  } else { # in case there is only 1 trophic level per station:
    trait_combinations <- 0
    
    trait_combinations_df <- data.frame(Trait1 = 0,
                                        Trait2 = 0)
  }
  
  # trait_combinations <- unique(t(combn(unique(station_year_data$Trophic), 2)))}
  
  # trait_combinations_df <- data.frame(Trait1 = trait_combinations[, 1],
  #                                     Trait2 = trait_combinations[, 2])
  
  
  trait_combinations_df$ratio <- sapply(1:nrow(trait_combinations_df), function(x) {
    
    trait1_abundance <- sum(station_year_data$Abundance[station_year_data$Trophic == trait_combinations_df$Trait1[x]])
    
    trait2_abundance <- sum(station_year_data$Abundance[station_year_data$Trophic == trait_combinations_df$Trait2[x]])
    
    ratio <- trait1_abundance / trait2_abundance
    return(ratio)
  })
  
  Results <- cbind(Sample = station, trait_combinations_df)
  
  combined_data <- rbind(combined_data, Results)
}



# Calculate the the geometric mean of ratios
agg_ratio <- combined_data %>%
  group_by(Sample) %>%
  summarize(mean_ratio = exp(mean(log(ratio))))

#Merging with interactions
Interactiv_indicies <- merge(co_occurence, agg_ratio, by = "Sample", all = TRUE)


print("3_BioQ - Done")


