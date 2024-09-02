################################################################################
# Process nearly-private rare allele data
################################################################################

# To extract nearly-private rare allele similarity, download raw files from 
# Fontsere et al. 2022: GitHub repository from Martin Kuhlwilm in data
# https://github.com/kuhlwilm/rareCAGA
# Load files from Data named "aprivatT9" and "coln"

# PROCESS SNP DATA

  # get all loci and nearly-private rare alleles (NePRA): 1 if NePRA, 0 if no NePRA, NA if not sequenced
  m.snps <-do.call(rbind,apri)
  table(m.snps, useNA = "always")
  
  # add population names to columns
  colnames(m.snps)<-coln

  # remove populations we are not interested in: Campo Ma'an does not have sufficient information
  # and Chinko and Tayna are not in behavioural data
  m.snps <- m.snps[,-which(coln %in% c("Campo_Ma'an", "Chinko", "Tayna"))] 

  # convert to data frame
  d.snps <- as.data.frame(m.snps)
  colnames(d.snps)

  # rename to match other data files
  colnames(d.snps)[c(3,11,16,20,24:26,29,30,34)] <- 
  c("Bili-Uele","East Nimba","Issa-Ugalla", "La Belgique", "Mt. Cameroon",
    "Mt. Sangbe","Mt. Cristal","Outamba-Kilimi","Rubi-Tele", "Bakoun-Sobory")


# GET DYADIC NEPRA SIMILARITY
  
  # Function to extract number of shared NePRA (ones) and total number of sequenced SNPs (without NA)
  get_nepra <- function(data){
    
    x <- data.frame(site1 = rep(colnames(data), each = ncol(data)))
    x$site2 <- rep(colnames(data), times = ncol(data))
    x$dyad_ID <- paste0(pmin(x$site1, x$site2), "_",
                       pmax(x$site1, x$site2))
    x$ones <- NA
    x$total <- NA
    c = 0
    for (i in 1:ncol(data)) {
      for (j in 1:ncol(data)) {
        c = c+1
        x[c,4] <- length(which(data[,i] == 1 & data[,j] ==1))
        x[c,5] <- length(which(data[,i] %in% c(0,1) & data[,j] %in% c(0,1)))
        print(x[c,])
        
      }
    }
    # subset from each snp to dyads
    x <- distinct(x, dyad_ID, .keep_all = T)
    return(x)
  }

  # Run function
  d.nepra <- get_nepra(d.snps)
  
  # calcluate shared nr of nepra by total number of sequenced snps in both sites
  
  d.nepra$nepra_proportion <- d.nepra$ones / d.nepra$total
  
  # since we are interested in between-population properties, exclude dyads from the same site
  d.nepra <- filter(d.nepra, site1 != site2)
  
  # save file
  d.nepra <- select(d.nepra, c("site1", "site2", "dyad_ID", "nepra_proportion"))
  write.table(d.nepra, "Data/d.nepra.txt", sep = "\t")
  