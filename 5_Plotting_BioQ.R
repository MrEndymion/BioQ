#========================================================================
#
# STANDARDIZATION OF INDICIE
#
# by Lars O. Mortensen, DHI Okt 2021
#
#========================================================================

#=============================================
# Adding variables
#=============================================

BioQ_scores$Distance <- unlist(lapply(strsplit(BioQ_scores$Sample, " "), "[", 2))
BioQ_scores$Distance <- ifelse(is.na(BioQ_scores$Distance), "Ref", BioQ_scores$Distance)
BioQ_scores$Distance <- factor(BioQ_scores$Distance, levels = c("100", "250", "750", "1500", "3000", "5000", "Ref"))


#=============================================
# Summarising plots
#=============================================

#Tax summary
g <- ggplot(BioQ_scores, aes(Distance, Tax_mean)) +
  geom_boxplot() +
  xlab("Distance") +
  ylab("Std alpha diveristy") +
  theme_bw()
g
ppi <- 200
ggsave(file.path(results.path, "Plots", paste0("Species_per_year_", Platform_pos, ".png")), g, width=10*ppi, height=8*ppi, units = "px")

#Fun summary
g <- ggplot(BioQ_scores, aes(Distance, Fun_mean)) +
  geom_boxplot() +
  xlab("Distance") +
  ylab("Std beta diveristy") +
  theme_bw()
g
ppi <- 200
ggsave(file.path(results.path, "Plots", paste0("Individuals_per_year_", Platform_pos, ".png")), g, width=10*ppi, height=8*ppi, units = "px")

#Inc summary
g <- ggplot(BioQ_scores, aes(Distance, Inc_mean)) +
  geom_boxplot() +
  xlab("Distance") +
  ylab("Std beta diveristy") +
  theme_bw()
g
ppi <- 200
ggsave(file.path(results.path, "Plots", paste0("Individuals_per_year_", Platform_pos, ".png")), g, width=10*ppi, height=8*ppi, units = "px")


#=============================================
#Plotting BioQ per year
#=============================================

#Impact plots
g <- ggplot(BioQ_scores, aes(Distance, BioQ)) +
  geom_boxplot() +
  xlab("Distance") +
  ylab("BioQ score") +
  theme_bw()
g
ppi <- 200
ggsave(file.path(results.path, "Plots", paste0("Individuals_per_year_", Platform_pos, ".png")), g, width=10*ppi, height=8*ppi, units = "px")

#=============================================
#Individual 
#=============================================


#Map
g <- ggplot() +
  geom_point(data = EBM_scores_pos, aes(x = x, y = y, col = BioQ), size = 4) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "lightblue"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "lightblue")) +
  xlab("Longitude (degrees)") +
  ylab("Latitude (degrees)") +
  ggtitle("EBM BioQ score per station") +
  scale_color_gradient(low = "red", high = "green", name = "EBM", limits = c(c(min(EBM_scores_pos$BioQ, na.rm=TRUE), 
                                                                               c(max(EBM_scores_pos$BioQ, na.rm=TRUE))))) +
  facet_wrap(~Year)
g
ggsave(file.path(results.path, "Plots", paste0("BioQ_per_year_", Platform_pos, ".png"), g))


#=============================================
#Plotting for distance
#=============================================
#Box_plot
g <- ggplot(EBM_scores_pos, aes(Distance, BioQ)) +
  geom_boxplot() +
  xlab("Year") +
  ylab("Standardized interactions diveristy Score") +
  theme_bw() +
  ggtitle("Ecosystem Scores") +
  facet_wrap(~Year)
g

ggsave(file.path(results.path, "Plots", "Ecosystem_scores_distance.png"), g)

#=============================================
# Transforming to raster
#=============================================
#Getting HD as well

Spat_EBM <- vect(EBM_scores_pos, geom=c("x", "y"), crs = "EPSG:4326")

Current_speed <- raster("L:\\Data\\HD_data_North_Sea\\summer\\hdukns_jncc_mean_cs_apr-sep18.asc")

Spat_EBM$Distance <- distance(Spat_EBM, Dan_f_center)
EBM_scores_model <- as.data.frame(Spat_EBM, geom = "xy")

EBM_scores_model <- subset(EBM_scores_model, Year == 2018)

test_spat <- gam(Mean ~ s(Distance, k = 5) + 
                   s(x, y, by = Year, k = 10), data = EBM_scores_model)
summary(test_spat)
plot(test_spat,scale=0,all.terms=T,pers=T,shade=T, page=1, rug = TRUE)


r <- rast()
crs(r) <- "EPSG:4326"
ext(r) <- ext(Spat_EBM)
res(r) <- 0.001

r_point <- as.points(r)

r_point$Distance <- distance(r_point, Dan_f_center)

r_point <- as.data.frame(r_point, geom = "xy")
r_point$Year <- 2015
r_point$pred <- predict(test_spat, r_point)

r_pred <- vect(r_point, geom=c("x", "y"), crs = "EPSG:4326")
r_pred$pred <- as.numeric(r_pred$pred)

r_pred <- rasterize(r_pred, r, "pred")

plot(r_pred)
wri

#Testing IDW moving window
r <- rast()
crs(r) <- "EPSG:4326"
ext(r) <- ext(Spat_EBM)
res(r) <- 0.001

EBM_IDW <- EBM_scores_pos_2018[,c("x", "y", "Mean")]
EBM_IDW <- cbind(EBM_IDW[,1], EBM_IDW[,2], EBM_IDW[,3])
test <- interpIDW(r, EBM_IDW, radius=0.1, power=3, smooth=3, maxPoints=2)


#Testing IDW
suppressMessages(library(spatstat))

obs_window <- owin(ext(Spat_EBM)[1:2], ext(Spat_EBM)[3:4])
ppp_EBM <- ppp(EBM_IDW$x, EBM_IDW$y, marks = EBM_IDW$Mean, window = obs_window)

idw_EBM <- idw(ppp_EBM, power=1, at="pixels")
idw_EBM_raster <- raster(idw_malaria)

plot(idw_malaria,
     col=heat.colors(20), 
     main="Interpolated EBM on IDW method \n (Power = 1)") 






