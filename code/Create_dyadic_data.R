################################################################################
# Mantel correlation tests
################################################################################

  # This code describes the process to combine the extracted and adapted data
  # from Kalan et al. 2020 and Fontsere et al. 2022 into a dyadic format

  library(dplyr)


## LOAD BEHAVIOURAL AND ECOLOGICAL DATA

  # Adapted data file from Kalan et al. 2020 for behavioural and ecological data, 
  # subset of 35 populations with genetic data, and adapted behaviours to 15 foraging behaviours (see table 1),
  # included variable "behavioral class" according to classification (see Materials and Methods)

  d.behaviours <- read.delim("Data/d.behaviours.txt")
  str(d.behaviours)
  
## CONVERT TO DYADIC DATA
  
  # Create dyads
  m.dyads <- data.frame(t(combn(unique(d.behaviours$site), 2))) # 595 dyads
  colnames(m.dyads) <- c("site.x", "site.y")
  
  # Create variable dyad ID
  m.dyads$dyad_ID <- paste0(pmin(m.dyads$site.x, m.dyads$site.y), "_",
                            pmax(m.dyads$site.x, m.dyads$site.y))
  
  # repeat each row for each behaviour
  m.dyads <- m.dyads[rep(1:nrow(m.dyads),each=15),]
  m.dyads$behavior <- rep(unique(d.behaviours$behavior), times = 595)
  
  # add behavioral class
  m.dyads <- merge(m.dyads, 
                   distinct(d.behaviours, behavior, .keep_all = T)[,c("behavior", "behavioral_class")],
                   by = "behavior")

  # add population-level information
  m.dyads <- merge(m.dyads, distinct(d.behaviours, site, .keep_all = T)[,-c(6:8)],
                   by.x = "site.x", by.y = "site", all.x = T)
  m.dyads <- merge(m.dyads, distinct(d.behaviours, site, .keep_all = T)[,-c(6:8)],
                   by.x = "site.y", by.y = "site", all.x = T)
  
  # add presence of behaviours
  m.dyads <- merge(m.dyads, d.behaviours[,c("site", "behavior", "presence")],
                   by.x = c("site.x", "behavior"), by.y = c("site", "behavior"), all.x = T)
  m.dyads <- merge(m.dyads, d.behaviours[,c("site", "behavior", "presence")],
                   by.x = c("site.y", "behavior"), by.y = c("site", "behavior"), all.x = T)

  # create dyadic-level variables
  # same habitat
  m.dyads$same_habitat <- ifelse(m.dyads$habitat.x == m.dyads$habitat.y, 1, 0)
  
  # mean distance to refuagia
  m.dyads$mean_refugia_distance <- rowMeans(m.dyads[, c("refugia_distance.x", "refugia_distance.y")])
  
  # mean annual precipitation
  m.dyads$mean_annual_precipitation <- rowMeans(m.dyads[, c("annual_precipitation.x", "annual_precipitation.y")])
  
  # minimum observation time
  m.dyads$min_observation_time <- pmin(m.dyads$observation_time.x, m.dyads$observation_time.y)
  
  # shared behaviour
  m.dyads$shared <- ifelse(m.dyads$presence.x == 1 & m.dyads$presence.y == 1, 1,
                           ifelse(is.na(m.dyads$presence.x) | is.na(m.dyads$presence.y), NA, 0))
  
  
## LOAD GENETIC AND GEOGRAPHIC DATA
  
  # genetic data on IBD compiled from Fontsere et al. 2022 Fig. 3, Fig. S93 and Table S9 (column IBD source)
  # genetic data on NePRA extracted and processed from the GitHub repository from
  # Fontsere et al. 2022: https://github.com/kuhlwilm/rareCAGA
  # Geographic coordinates and subspecies extracted from Fontsere et al. 2022 Table S1
  # and subset to the populations in our study
  
  d.ibd <- read.delim("Data/d.ibd.txt")
  d.nepra <- read.delim("Data/d.nepra.txt")
  d.geographic <- read.delim("Data/d.geographic_data.txt")
  
  # merge with dyadic data
  m.dyads <- merge(m.dyads, d.ibd[,c("dyad_ID", "ibd_binary")],
                   by = "dyad_ID", all.x = T)
  m.dyads <- merge(m.dyads, d.nepra[,c("dyad_ID", "nepra_proportion")],
                   by = "dyad_ID", all.x = T)
  
  m.dyads <- merge(m.dyads, d.geographic,
                   by.x = "site.x", by.y = "site", all.x = T)
  m.dyads <- merge(m.dyads, d.geographic,
                   by.x = "site.y", by.y = "site", all.x = T)
  
  
  
  # Geographic distance in km
  m.dyads$distance <- NA
  
  for (i in (1:nrow(m.dyads))) {
    d <- as.numeric(geosphere::distm(c(m.dyads$lon.x[i], m.dyads$lat.x[i]), c(m.dyads$lon.y[i], m.dyads$lat.y[i])))
    m.dyads$distance[i] <- d
  }
  remove(i, d)
  
  # same subspecies
  m.dyads$same_subspecies <- ifelse(m.dyads$subspecies.x == m.dyads$subspecies.y, 1, 0)

  # nepra_binary
  hist(m.dyads$nepra_proportion, breaks = 100) # threshold at 0.00006
  m.dyads$nepra_binary <- ifelse(m.dyads$nepra_proportion > 0.0006, 1, 0)

  
## ORGANISE DATA AND SAVE
  
  m.dyads <- m.dyads[,c("site.x", "site.y", "dyad_ID", "behavior", "behavioral_class",
                        "shared", "nepra_binary", "nepra_proportion", "ibd_binary",
                        "same_subspecies", "distance", "same_habitat", "mean_refugia_distance",
                        "mean_annual_precipitation", "min_observation_time")]
  colnames(m.dyads)[1:2] <- c("site1", "site2")

  write.table(m.dyads, "Data/m.dyads.txt", sep = "\t")      
  
  
  
