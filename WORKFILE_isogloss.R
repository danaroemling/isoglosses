# libraries
library(sf)
library(sp)
library(tidyverse)
library(gstat)
library(stringr) 
library(dplyr) 
library(scales)
library(classInt)
library(viridis)
library(viridisLite)
library(rnaturalearth)
library(spThin)
library(sfheaders)

# data
corpus <- read.csv(file = './data_ling/full_matrix_for_filtering.csv') 
token_at_location <- read.csv(file = './data_ling/tokens_at_location.csv')
colnames(token_at_location) <- c("City", "Tokencount") 
longlat <- read.csv(file = './data_maps/gsa_geo_filtered.csv')
#longlat_small <- read.csv(file = './data_maps/gsa_geo_filtered_nosmall22.csv')
#longlat_small <- read.csv(file = './data_maps/gsa_geo_filtered_nosmall85.csv')

# cities
cities <- data.frame(
  City = c("Köln", "München", "Wien", "Zürich", "Berlin", "Hamburg"),
  Long = c(6.9578, 11.5755, 16.3731, 8.5417, 13.3833, 10),
  Lat = c(50.9422, 48.1372, 48.2083, 47.3769, 52.5167, 53.55))
crs2 <- CRS("+init=epsg:4326")
cities_sf <- st_as_sf(cities, coords = c("Long", "Lat"), crs = crs2)

# with usual code / alternative
gsa_outline <- ne_countries(country = c("Austria", "Germany", "Switzerland"), returnclass="sf", scale = "large")
gsa_plot <- gsa_outline %>% dplyr::select(geometry) # just the outlines of the three
gsa_spatial <- as_Spatial(gsa_plot) # new object as spatial instead of sf

# large grid
sp_grid <- as.data.frame(spsample(gsa_spatial,
                                  n = 100000,
                                  type = "regular", # systematically aligned sampling
                                  offset = c(0.5,0.5))) # makes sure grid is the same every time
# This grid needs A CRS and to be turned into a spatial object, which this line does:
sp_grid_sf <- st_as_sf(sp_grid, coords=c("x1","x2"), crs = st_crs(4326))

# small grid
sp_grid_small <- as.data.frame(spsample(gsa_spatial,
                                        n = 500,
                                        type = "regular", # systematically aligned sampling
                                        offset = c(0.5,0.5))) # makes sure grid is the same every time
# This grid needs A CRS and to be turned into a spatial object, which this line does:
sp_grid_sf_small <- st_as_sf(sp_grid_small, coords=c("x1","x2"), crs = st_crs(4326))


# corpus search
# alternations
one_word <- corpus %>% filter(word == "schau")
colnames(one_word) <- c("Token", "City", "Frequency")
second_word <- corpus %>% filter(word == "guck")
colnames(second_word) <- c("Token", "City", "Frequency")

# combine into one and calculate proportion
merged <- merge(one_word, second_word, by ="City", all=T)
# introduce zeros for NAs
merged$Frequency.x[is.na(merged$Frequency.x)] <- 0
merged$Frequency.y[is.na(merged$Frequency.y)] <- 0
merged$Proportion <- (merged$Frequency.x / (merged$Frequency.x + merged$Frequency.y) )
merged$Percentage <- 100* (merged$Frequency.x / (merged$Frequency.x + merged$Frequency.y) )

# merge with geo data
merger <- merge(merged, longlat, by.x ="City", by.y = "City")
#write.csv(merger, "isogloss_data.csv", row.names=FALSE)


# variogram & kriging
trial <- merger %>% dplyr::select("lon", "lat", "Proportion") 
data <- trial[3]
coords <- trial[1:2]
crs <- CRS("+init=epsg:4326") 
spdf <- SpatialPointsDataFrame(coords = coords, data = data, proj4string = crs)
spdf = spdf[which(!duplicated(spdf@coords)), ]
vg <- variogram(Proportion~1, data = spdf, width = 1, cutoff = 300)
vg_fit <- fit.variogram(vg, vgm("Exp")) 
krig_comp_1 <- krige(Proportion~1, locations = spdf, newdata = sp_grid_sf_small, model = vg_fit)
krig_comp_2 <- krige(var1.pred~1, locations = krig_comp_1, newdata = sp_grid_sf, model = vg_fit)

# mapping
# set legends
legend_max_c <- round(max(krig_comp_2$var1.pred), digits = 2)
legend_min_c <- round(min(krig_comp_2$var1.pred), digits = 2)
summary(krig_comp_2$var1.pred)


# contourline with one breakpoint
contour <- krig_comp_2 %>% dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                                         lat = sf::st_coordinates(.)[,2])



# light mode map
tiff("schau_guck_comparison_sixth_isogloss.tiff", units="in", width=3.78, height=4.8, res=300)

ggplot() +
  geom_sf(data = krig_comp_2, aes(fill=var1.pred), shape = 21, size = 0.75, stroke = 0, lwd = 0) +
  geom_sf(data = gsa_plot, aes(geometry = geometry), color="snow2", fill=NA, size = 0.5) +
  geom_contour(data = contour, aes(x = lon, y = lat, z = var1.pred, colour = "red"), stat = "contour", breaks = c(0.5), position = "identity", show.legend = FALSE) +
  geom_sf_text(data = cities_sf, aes(label = City), size = 2.5, nudge_x = 0, nudge_y = -0.15, color = "snow", family = "Optima") +
  geom_sf(data = cities_sf, aes(geometry = geometry), shape = 4, color = "snow") +
  theme_minimal() +
  scale_fill_steps2(low = "#ECE74F", mid = "#111111", high = "#73ffc2", midpoint = 0.5, n.breaks = 15, labels = NULL) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank()) +
  theme(legend.position = c(0.90, 0.65), 
        legend.title = element_text(size = 6, color = "snow", family = "Optima"),
        legend.text = element_text(size = 6, color = "snow", family = "Optima", hjust = 1),
        legend.title.align = 0,
        panel.background = element_rect(fill = "gray10", color = "gray10"),
        plot.background = element_rect(fill = "gray10", color = "gray10"),
        legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.4,"cm")) +
  annotate(geom="text", x = 16.6, y = 52.4, label=legend_max_c, size = 1.5, family = "Optima", color = "snow") +
  annotate(geom="text", x = 16.6, y = 51.3, label=legend_min_c, size = 1.5, family = "Optima", color = "snow") +
  labs(fill = "brötchen / weggli") 


dev.off()



# map in dark mode blue/red
tiff("schau_smallgrid_isogloss.tiff", units="in", width=3.78, height=4.8, res=300)

ggplot() +
  geom_sf(data = krig_comp_1, aes(fill=var1.pred), shape = 21, size = 1.5, stroke = 0, lwd = 0) +
  geom_sf(data = gsa_plot, aes(geometry = geometry), color="snow2", fill=NA, size = 0.5) +
  geom_contour(data = contour, aes(x = lon, y = lat, z = var1.pred), stat = "contour", colour = "yellow", breaks = c(0.5), position = "identity", show.legend = FALSE) +
  geom_sf_text(data = cities_sf, aes(label = City), size = 2.5, nudge_x = 0, nudge_y = -0.15, color = "snow", family = "Optima") +
  geom_sf(data = cities_sf, aes(geometry = geometry), shape = 4, color = "snow") +
  theme_minimal() +
  scale_fill_steps2(low = "blue", mid = "#303030", high = "firebrick1", midpoint = 0.5, n.breaks = 15, labels = NULL) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank()) +
  theme(legend.position = c(0.90, 0.65), 
        legend.title = element_text(size = 6, color = "snow", family = "Optima"),
        legend.text = element_text(size = 6, color = "snow", family = "Optima", hjust = 1),
        legend.title.align = 0,
        panel.background = element_rect(fill = "#303030", color = "#303030"),
        plot.background = element_rect(fill = "#303030", color = "#303030"),
        legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.4,"cm")) +
  annotate(geom="text", x = 16.6, y = 52.3, label=legend_max_c, size = 1.5, family = "Optima", color = "snow") +
  annotate(geom="text", x = 16.6, y = 51.15, label=legend_min_c, size = 1.5, family = "Optima", color = "snow") +
  labs(fill = "schau / guck") 


dev.off()




# map in dark mode blue/red isogloss no fill
tiff("schau_isogloss_smallgrid_only.tiff", units="in", width=3.78, height=4.8, res=300)

ggplot() +
  geom_sf(data = gsa_plot, aes(geometry = geometry), color="snow2", fill=NA, size = 0.5) +
  geom_contour(data = contour, aes(x = lon, y = lat, z = var1.pred), stat = "contour", colour = "yellow", breaks = c(0.5), position = "identity", show.legend = FALSE) +
  geom_sf_text(data = cities_sf, aes(label = City), size = 2.5, nudge_x = 0, nudge_y = -0.15, color = "snow", family = "Optima") +
  geom_sf(data = cities_sf, aes(geometry = geometry), shape = 4, color = "snow") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank()) +
  theme(panel.background = element_rect(fill = "#303030", color = "#303030"),
        plot.background = element_rect(fill = "#303030", color = "#303030")) +
  #annotate(geom="text", x = 10.7, y = 48.7, label= "semmel", size = 2, family = "Optima", color = "snow", fontface = "italic") +
  #annotate(geom="text", x = 15.5, y = 47.6, label= "mistkübel", size = 2, family = "Optima", color = "snow", fontface = "italic") 
  annotate(geom="text", x = 10.5, y = 50., label= "schau", size = 2, family = "Optima", color = "snow", fontface = "italic") +
  annotate(geom="text", x = 8.4, y = 52.3, label= "guck", size = 2, family = "Optima", color = "snow", fontface = "italic") 

dev.off()




# map in light mode
jpeg("schau_guck_noline.jpg", units="in", width=3.78, height=4.8, res=300)

ggplot() +
  geom_sf(data = krig_comp_2, aes(fill=var1.pred), shape = 21, size = 0.75, stroke = 0, lwd = 0) +
  geom_sf(data = gsa_plot, aes(geometry = geometry), color="black", fill=NA, size = 0.5) +
  #geom_contour(data = contour, aes(x = lon, y = lat, z = var1.pred, colour = "red"), stat = "contour", breaks = c(0.5), position = "identity", show.legend = FALSE) +
  geom_sf_text(data = cities_sf, aes(label = City), size = 2.5, nudge_x = 0, nudge_y = -0.15, family = "Optima") +
  geom_sf(data = cities_sf, aes(geometry = geometry), shape = 4) +
  theme_minimal() +
  scale_fill_steps2(low = "#E0961A", mid = "#ffffff", high = "#310d59", midpoint = 0.5, n.breaks = 15, labels = NULL) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank()) +
  theme(legend.position = c(0.90, 0.65), 
        legend.title = element_text(size = 6, family = "Optima"), 
        legend.text = element_text(size = 6, family = "Optima", hjust = 1),
        legend.title.align = 0,
        legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.4,"cm")) +
  annotate(geom="text", x = 16.7, y = 52.35, label=legend_max_c, size = 1.5, family = "Optima") +
  annotate(geom="text", x = 16.7, y = 51.25, label=legend_min_c, size = 1.5, family = "Optima") +
  labs(fill = "schau/guck")
  #labs(fill = "mülleimer/\nmistkübel")


dev.off()









