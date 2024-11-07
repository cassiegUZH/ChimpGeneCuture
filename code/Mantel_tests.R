################################################################################
# Mantel correlation tests
################################################################################

  # Code to run mantel tests from distance matrices for each variable
  # using "vegan" and "ecodist" packages

rm(list=ls())
#library(dplyr)
'%!in%' <- function(x,y)!('%in%'(x,y))

################################################################################
## LOAD MATRICES -----

## GENETIC
  m.nepra <- as.dist(read.delim("Matrix_NePRA.txt"), diag=F, upper=F)
  m.ibd <- as.dist(read.delim("Matrix_IBD.txt", row.names = 1), diag=F, upper=F)

## GEOGRAPHIC
  m.region <- as.dist(read.delim("Matrix_same_subspecies.txt", row.names = 1), diag=F, upper=F)
  m.distance <- as.dist(read.delim("Matrix_distances.txt", row.names = 1), diag=F, upper=F)

## ECOLOGICAL
  m.habitat <- as.dist(read.delim("Matrix_same_habitat.txt", row.names = 1), diag=F, upper=F)
  m.refugia <- as.dist(read.delim("Matrix_distance_refugia.txt", row.names = 1), diag=F, upper=F)
  m.precipitation <- as.dist(read.delim("Matrix_annual_precipitation.txt", row.names = 1), diag=F, upper=F)

## OBSERVATION TIME
  m.time <- as.dist(read.delim("Matrix_observation_time.txt", row.names = 1), diag=F, upper=F)


## BEHAVIOURS
  # Load affiliation matrices
  m.nontool <- read.delim("Matrix_affiliation_nontool.txt", row.names = 1)
  m.simple <- read.delim("Matrix_affiliation_simple.txt", row.names = 1)
  m.complex <- read.delim("Matrix_affiliation_complex.txt", row.names = 1)

  # to calculate behavioural distance between populations, apply Jaccard dissimilarity 
  # from the vegan package
  library(vegan)
  m.nontool <- vegdist(m.nontool, method = "jaccard", diag = F, upper = F)
  m.nontool[is.na(m.nontool)] <- 1 # pairs that don't have any behaviours coded as NA, we want maximum dissimilarity (1)
  # create "similarity matrix" for more intuitive results
  m.nontool <- 1 - m.nontool

  # repeat for simple and complex
  m.simple <- vegdist(m.simple, method = "jaccard", diag = F, upper = F)
  m.simple[is.na(m.simple)] <- 1
  m.simple <- 1 - m.simple
  
  m.complex.original <- m.complex
  m.complex <- vegdist(m.complex, method = "jaccard", diag = F, upper = F)
  m.complex[is.na(m.complex)] <- 1
  m.complex <- 1 - m.complex
  
  # get matrices for each complex behaviour (can be repeated for simple/non-tool)

  m.honey <- vegdist(m.complex.original[,1], method = "jaccard", diag = F, upper = F)
  m.honey[is.na(m.honey)] <- 1
  m.honey <- 1 - m.honey
  
  m.nut <- vegdist(m.complex.original[,2], method = "jaccard", diag = F, upper = F)
  m.nut[is.na(m.nut)] <- 1
  m.nut <- 1 - m.nut
  
  m.underground <- vegdist(m.complex.original[,3], method = "jaccard", diag = F, upper = F)
  m.underground[is.na(m.underground)] <- 1
  m.underground <- 1 - m.underground


################################################################################
## RUN MANTEL TESTS -----
  
  # detach vegan, and get ecodist for mantel tests
  detach("package:vegan", unload=TRUE)
  library(ecodist)


## SIMPLE MANTEL TESTS

  # prepare list of all matrices that are tested pair-wise
  matrix.list <- list(m.habitat, m.refugia, m.precipitation, 
                      m.region, m.distance,
                      m.nepra, m.ibd, m.nontool, m.simple, m.complex)
  
  # number of matrices to be tested
  n = length(matrix.list)
  
  # prepare matrices to store results on Mantel correlation coefficient r, 
  # CI and p-value
  matrix.r <-matrix(1,nrow=n,ncol=n)
  rownames(matrix.r)<- c('Habitat', "Distance to\n refugia", "Annual\n precipitation", 
                         'Region','Geographic\ndistance', #"Observation\ntime",
                         'NePRA', "IBD", 'Non-tool','Simple','Complex')
  
  colnames(matrix.r)<- c('Habitat', "Distance to\n refugia", "Annual\n precipitation", 
                         'Region','Geographic\ndistance', #"Observation\ntime",
                         'NePRA', "IBD", 'Non-tool','Simple','Complex')
  
  matrix.ci2.5 <- matrix.r
  matrix.ci97.5 <- matrix.r
  matrix.pval <- matrix.r
  
  # Loop through all matrices in matrix.list for each pair-wise test
  
  for(i in 1:n){
    k = i
    matrix1 = matrix.list[[i]]
    while (k<n){
      k <- k+1
      matrix2 = matrix.list[[k]]
      mantel.test <- mantel(matrix1 ~ matrix2, nperm = 500, mrank = T)
      # store test results
      matrix.r[k,i] <- mantel.test[1]
      matrix.r[i,k] <- mantel.test[1] # get symmetric results
      matrix.ci2.5[k,i] <- mantel.test[5]
      matrix.ci2.5[i,k] <- mantel.test[5]
      matrix.ci97.5[k,i] <- mantel.test[6]
      matrix.ci97.5[i,k] <- mantel.test[6]
      matrix.pval[k,i] <- mantel.test[4]
      matrix.pval[i,k] <- mantel.test[4]
    }
    # store all results in a list
    l.results <- list(matrix.r, matrix.ci2.5, matrix.ci97.5, matrix.pval)
  }
  
  


## SIMPLE MANTEL TESTS FOR COMPLEX BEHAVIOURS

  # repeat above for simple Mantel test including individual complex behaviours
  matrix.list <- list(m.honey, m.nut, m.underground,
                      m.nepra, m.ibd)
  n = length(matrix.list)
  
  matrix.r <-matrix(1,nrow=n,ncol=n)
  rownames(matrix.r)<- c("Honey dipping\ntoolset", "Nut cracking", "Underground\ntoolset", 'NePRA','IBD')
  colnames(matrix.r)<- c("Honey dipping\ntoolset", "Nut cracking", "Underground\ntoolset", 'NePRA','IBD')
  
  matrix.ci2.5 <- matrix.r
  matrix.ci97.5 <- matrix.r
  matrix.pval <- matrix.r
  
  
  for(i in 1:n){
    k = i
    matrix1 = matrix.list[[i]]
    while (k<n){
      k <- k+1
      matrix2 = matrix.list[[k]]
      mantel.test <- mantel(matrix1 ~ matrix2, nperm = 500, mrank = T)
      matrix.r[k,i] <- mantel.test[1]
      matrix.r[i,k] <- mantel.test[1]
      matrix.ci2.5[k,i] <- mantel.test[5]
      matrix.ci2.5[i,k] <- mantel.test[5]
      matrix.ci97.5[k,i] <- mantel.test[6]
      matrix.ci97.5[i,k] <- mantel.test[6]
      matrix.pval[k,i] <- mantel.test[4]
      matrix.pval[i,k] <- mantel.test[4]
    }
    l.results.complex <- list(matrix.r, matrix.ci2.5, matrix.ci97.5, matrix.pval)
  }

# PARTIAL MANTEL CONTROL FOR OBSERVATION TIME

  matrix.list <- list(m.complex, m.simple, m.nontool, m.nepra, m.ibd)
  n = length(matrix.list)
  
  matrix.r <-matrix(1,nrow=n,ncol=n)
  rownames(matrix.r)<- c("Complex", "Simple", "Non-tool", 'NePRA','IBD')
  
  colnames(matrix.r)<- c("Complex", "Simple", "Non-tool", 'NePRA','IBD')
  
  # matrix.ci2.5 <- matrix.r
  # matrix.ci97.5 <- matrix.r
  matrix.pval <- matrix.r
  
  
  for(i in 1:n){
    k = i
    matrix1 = matrix.list[[i]]
    while (k<n){
      k <- k+1
      matrix2 = matrix.list[[k]]
      mantel.test <- mantel(matrix1 ~ matrix2 + m.precipitation, nperm = 500, mrank = T)
      # STORE TEST RESULTS
      matrix.r[k,i] <- mantel.test[1]
      matrix.r[i,k] <- mantel.test[1]
      # matrix.ci2.5[k,i] <- mantel.test[5]
      # matrix.ci2.5[i,k] <- mantel.test[5]
      # matrix.ci97.5[k,i] <- mantel.test[6]
      # matrix.ci97.5[i,k] <- mantel.test[6]
      matrix.pval[k,i] <- mantel.test[4]
      matrix.pval[i,k] <- mantel.test[4]
    }
    l.results.obstime <- list(matrix.r, matrix.ci2.5, matrix.ci97.5, matrix.pval)
  }

# PARTIAL MANTEL CONTROL FOR HABITAT

  matrix.list <- list(m.complex, m.simple, m.nontool,
                      m.nepra, m.ibd)
  n = length(matrix.list)
  
  matrix.r <-matrix(1,nrow=n,ncol=n)
  rownames(matrix.r)<- c("Complex", "Simple", "Non-tool", 'NePRA','IBD')
  
  colnames(matrix.r)<- c("Complex", "Simple", "Non-tool", 'NePRA','IBD')
  
  # matrix.ci2.5 <- matrix.r
  # matrix.ci97.5 <- matrix.r
  matrix.pval <- matrix.r
  
  
  for(i in 1:n){
    k = i
    matrix1 = matrix.list[[i]]
    while (k<n){
      k <- k+1
      matrix2 = matrix.list[[k]]
      mantel.test <- mantel(matrix1 ~ matrix2 + m.habitat, nperm = 500, mrank = T)
      # STORE TEST RESULTS
      matrix.r[k,i] <- mantel.test[1]
      matrix.r[i,k] <- mantel.test[1]
      # matrix.ci2.5[k,i] <- mantel.test[5]
      # matrix.ci2.5[i,k] <- mantel.test[5]
      # matrix.ci97.5[k,i] <- mantel.test[6]
      # matrix.ci97.5[i,k] <- mantel.test[6]
      matrix.pval[k,i] <- mantel.test[4]
      matrix.pval[i,k] <- mantel.test[4]
    }
    l.results.habitat <- list(matrix.r, matrix.ci2.5, matrix.ci97.5, matrix.pval)
  }


# PARTIAL MANTEL CONTROL FOR REFUGIA

  matrix.list <- list(m.complex, m.simple, m.nontool,
                      m.nepra, m.ibd)
  n = length(matrix.list)
  
  matrix.r <-matrix(1,nrow=n,ncol=n)
  rownames(matrix.r)<- c("Complex", "Simple", "Non-tool", 'NePRA','IBD')
  
  colnames(matrix.r)<- c("Complex", "Simple", "Non-tool", 'NePRA','IBD')
  
  # matrix.ci2.5 <- matrix.r
  # matrix.ci97.5 <- matrix.r
  matrix.pval <- matrix.r
  mantel(m.complex ~ m.ibd, mrank = T, nperm = 500)
  
  for(i in 1:n){
    k = i
    matrix1 = matrix.list[[i]]
    while (k<n){
      k <- k+1
      matrix2 = matrix.list[[k]]
      mantel.test <- mantel(matrix1 ~ matrix2 + m.refugia, nperm = 500, mrank = T)
      # STORE TEST RESULTS
      matrix.r[k,i] <- mantel.test[1]
      matrix.r[i,k] <- mantel.test[1]
      # matrix.ci2.5[k,i] <- mantel.test[5]
      # matrix.ci2.5[i,k] <- mantel.test[5]
      # matrix.ci97.5[k,i] <- mantel.test[6]
      # matrix.ci97.5[i,k] <- mantel.test[6]
      matrix.pval[k,i] <- mantel.test[4]
      matrix.pval[i,k] <- mantel.test[4]
    }
    l.results.refugia <- list(matrix.r, matrix.ci2.5, matrix.ci97.5, matrix.pval)
  }

# PARTIAL MANTEL CONTROL FOR REFUGIA

  matrix.list <- list(m.complex, m.simple, m.nontool,
                      m.nepra, m.ibd)
  n = length(matrix.list)
  
  matrix.r <-matrix(1,nrow=n,ncol=n)
  rownames(matrix.r)<- c("Complex", "Simple", "Non-tool", 'NePRA','IBD')
  
  colnames(matrix.r)<- c("Complex", "Simple", "Non-tool", 'NePRA','IBD')
  
  # matrix.ci2.5 <- matrix.r
  # matrix.ci97.5 <- matrix.r
  matrix.pval <- matrix.r
  
  
  for(i in 1:n){
    k = i
    matrix1 = matrix.list[[i]]
    while (k<n){
      k <- k+1
      matrix2 = matrix.list[[k]]
      mantel.test <- mantel(matrix1 ~ matrix2 + m.precipitation, nperm = 500, mrank = T)
      # STORE TEST RESULTS
      matrix.r[k,i] <- mantel.test[1]
      matrix.r[i,k] <- mantel.test[1]
      # matrix.ci2.5[k,i] <- mantel.test[5]
      # matrix.ci2.5[i,k] <- mantel.test[5]
      # matrix.ci97.5[k,i] <- mantel.test[6]
      # matrix.ci97.5[i,k] <- mantel.test[6]
      matrix.pval[k,i] <- mantel.test[4]
      matrix.pval[i,k] <- mantel.test[4]
    }
    l.results.precipitation <- list(matrix.r, matrix.ci2.5, matrix.ci97.5, matrix.pval)
  }


################################################################################
## VISUALISATION -----

## PREPARATIONS
  library(corrplot)
  
  # colour palette
  palet <- colorRampPalette(c(shades::saturation("tomato2",1), shades::saturation("antiquewhite",0.05) ,"deepskyblue"))

  # dev.off()
  
## SIMPLE TEST
  jpeg(filename = "Figure S7.jpeg", width = 20, height = 21, 
       units = "cm", 
       res = 300)
  layout(matrix(c(1, 2, 3, 4,5,6), nrow = 3, ncol = 2), widths = c(12, 8), heights = c(10,5,5,5,5,5))
  
  corrplot(l.results[[1]], p.mat = l.results[[4]], 
           diag=FALSE, method = 'pie', type = 'lower', col = palet(200), tl.col = 1, tl.cex = 0.9, 
           insig='label_sig',
           sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
           #cl.cex = 0.7, cl.offset = 0.5, 
           cl.pos = 'r',
           mar = c(0.2,0,2,0.2)) #mar = c(b,l,t,r)
  title("A) Simple Mantel tests", font.main = 2, cex.main = 1, adj = 0)


# plot 3: partial observation time
corrplot(l.results.obstime[[1]], p.mat = l.results.obstime[[4]], 
         diag=FALSE, method = 'pie', type = 'upper', col = palet(200), tl.col = 1, tl.cex = 0.9, 
         insig='label_sig',sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         cl.pos = 'n', mar = c(0.2,0,3.5,0.2))
title("C) Partial controlling observation time", font.main = 2, cex.main = 1, adj = 0)


# plot 5: partial proximity
corrplot(l.results.refugia[[1]], p.mat = l.results.refugia[[4]], 
         diag=FALSE, method = 'pie', type = 'upper', col = palet(200), tl.col = 1, tl.cex = 0.9, 
         insig='label_sig',sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         cl.pos = 'n', mar = c(0.2,0,3.5,0))
title("E) Partial controlling distance to refugia", font.main = 2, cex.main = 1, adj = 0)


# plot 2: complex
corrplot(l.results.complex[[1]], p.mat = l.results.complex[[4]], 
         diag=FALSE, method = 'pie', type = 'upper', col = palet(200), tl.col = 1, tl.cex = 0.9, 
         insig='label_sig',sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         cl.pos = 'n', mar = c(4.5,3.5,2,3.5))


title("B) Simple with complex tools", font.main = 2, cex.main = 1, adj = 0)

# plot 4: habitat type
corrplot(l.results.habitat[[1]], p.mat = l.results.habitat[[4]], 
         diag=FALSE, method = 'pie', type = 'upper', col = palet(200), tl.col = 1, tl.cex = 0.9, 
         insig='label_sig',sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         cl.pos = 'n', mar = c(0.2,0.2,3,0))
title("D) Partial controlling same habitat", font.main = 2, cex.main = 1, adj = 0)


# plot 6: region
corrplot(l.results.precipitation[[1]], p.mat = l.results.precipitation[[4]], 
         diag=FALSE, method = 'pie', type = 'upper', col = palet(200), tl.col = 1, tl.cex = 0.9, 
         insig='label_sig',sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         cl.pos = 'n', mar = c(0.2,0.2,3,0))
title("F) Partial controlling annual precipitation", font.main = 2, cex.main = 1, adj = 0)


layout(1)

dev.off()


