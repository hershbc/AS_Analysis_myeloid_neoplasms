

library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(matrixStats)


#This script will
#1) Join all rMATS data by all chromosome coordinates provided
#2) What is the frequency samples with no reads to support a splice junction?
#3) Pick a cutoff and filter table

###Read in rMATS data from all rMATS runs of a specific splicing type

rMATS <- list()
rMATS_gtf <- list()

setwd("/Volumes/padgetrlab/courtney/MLL/all/rmats_nostatt/")

rmatsfiles <- list.files("/Volumes/padgetrlab/courtney/MLL/all/rmats_nostatt/", pattern = "A3SS.MATS.JC.txt", recursive=TRUE)
gtffiles <- list.files("/Volumes/padgetrlab/courtney/MLL/all/rmats_nostatt/", pattern = "fromGTF.A3SS.txt", recursive=TRUE)
files1 <- list.files("/Volumes/padgetrlab/courtney/MLL/all/rmats_nostatt/", pattern = "file1.txt", recursive=TRUE)
files2 <- list.files("/Volumes/padgetrlab/courtney/MLL/all/rmats_nostatt/", pattern = "file2.txt", recursive=TRUE)


for(i in 1:length(rmatsfiles)){ 
  rMATS[[i]] <- read.csv(rmatsfiles[[i]], header=TRUE, sep="\t")
  
  rMATS_gtf[[i]] <- read.csv(gtffiles[[i]], header=TRUE, sep="\t")
  
  rMATS[[i]] <- merge(rMATS_gtf[[i]], rMATS[[i]], by = "ID")
  
}


###Unique Identifier for each splicing type


rMATS_2 <- rMATS
rMATS_2ID <- list()

if(grepl("A3SS|A5SS", rmatsfiles[[1]])){ 
  rMATS_2ID<- lapply(rMATS_2, function(x) paste(x$chr, x$strand, x$longExonStart_0base, x$longExonEnd, x$shortES, x$shortEE, x$flankingES, x$flankingEE, sep=","))
}

if(grepl("MXE", rmatsfiles[[1]])){ 
  rMATS_2ID<- lapply(rMATS_2, function(x) paste(x$chr, x$strand, x$X1stExonStart_0base, x$X1stExonEnd, x$X2ndExonStart_0base, x$X2ndExonEnd, x$upstreamES, x$upstreamEE, x$downstreamES, x$downstreamEE, sep=","))
}

if(grepl("RI", rmatsfiles[[1]])){ 
  rMATS_2ID<- lapply(rMATS_2, function(x) paste(x$chr, x$strand, x$riExonStart_0base, x$riExonEnd, x$upstreamES, x$upstreamEE, x$downstreamES, x$downstreamEE, sep=","))
}

if(grepl("SE", rmatsfiles[[1]])){ 
  rMATS_2ID<- lapply(rMATS_2, function(x) paste(x$chr, x$strand, x$exonStart_0base, x$exonEnd, x$upstreamES, x$upstreamEE, x$downstreamES, x$downstreamEE, sep=","))
}


#Create coord column
for (i in 1:length(rMATS_2)){ 
  rMATS_2[[i]]$coord <- rMATS_2ID[[i]]
}




###Split comma delimited columns so each value is in a column
#This needs done for IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2 and SJC_SAMPLE_2

#Parameters: columns, split into new named columns, separate on commas, remove column from dataframe, convert into integer, warn if there are more commas than new column names, warn if there are too few commas for new column names.

#Join Condition column from this dataset with control column from previous dataset

condition_IJC <- list()
condition_SJC <- list()
controls_IJC <- list()
controls_SJC <- list()
IncLevel1 <- list()
IncLevel2 <- list()


for(df in 1:length(rMATS_2)){ 
  setwd("/Volumes/padgetrlab/courtney/MLL/all/rmats_nostatt/")
  
  #Condition file names
  columnname1_file <- files1[[df]]
  columnname1 <- readChar(columnname1_file, file.info(columnname1_file)$size)
  
  #Control file names
  columnname2_file <- files2[[df]]
  columnname2 <- readChar(columnname2_file, file.info(columnname2_file)$size)
  
  list_of_split2 <- list()
  
  list_of_columns <- c("IJC_SAMPLE_1", "SJC_SAMPLE_1" , "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel1", "IncLevel2")
  
  
  for(i in 1:length(list_of_columns)){
    
    column <- list_of_columns[[i]]
    
    #If column name includes a 1, use list from file1, if not, use list from file2
    columnname <- ifelse(grepl("1",column), columnname1,columnname2)
    
    name <-strsplit(columnname, ",")
    
    name <- unlist(name)
    
    name <- gsub("(.+_.+?)(\\_.*)", "\\1", name)
    
    list_of_split2[[i]] <- rMATS_2[[df]][,c("coord",column)]
    
    list_of_split2[[i]] <- separate(list_of_split2[[i]], column, name, sep = ",", remove=TRUE, convert=TRUE, extra = "warn", fill = "warn")
    
    
  }
  
  condition_IJC[[df]] <- list_of_split2[[1]]
  condition_SJC[[df]] <- list_of_split2[[2]]
  IncLevel1[[df]] <- list_of_split2[[5]]
  
  controls_IJC[[df]] <- list_of_split2[[3]]
  controls_SJC[[df]] <- list_of_split2[[4]]
  IncLevel2[[df]] <- list_of_split2[[6]]
  
}





###Join the dataframes by coordinate to get counts table for all splicing changes


all_IJC1 <- Reduce(function(x, y) merge(x, y), condition_IJC)
rownames(all_IJC1) <- all_IJC1[,1]
all_IJC1 <- all_IJC1[,-1]

all_IJC2 <- Reduce(function(x, y) merge(x, y), controls_IJC)
rownames(all_IJC2) <- all_IJC2[,1]
all_IJC2 <- all_IJC2[,-1]

all_IJC <- merge(all_IJC1, all_IJC2, by="row.names")

all_SJC1 <- Reduce(function(x, y) merge(x, y), condition_SJC)
rownames(all_SJC1) <- all_SJC1[,1]
all_SJC1 <- all_SJC1[,-1]

all_SJC2 <- Reduce(function(x, y) merge(x, y), controls_SJC)
rownames(all_SJC2) <- all_SJC2[,1]
all_SJC2 <- all_SJC2[,-1]

all_SJC <- merge(all_SJC1, all_SJC2, by="row.names")

dim(all_IJC)
dim(all_SJC)


###Join dataframes by coordinates to get inclusion levels for all samples

#Take coordinate information and add Inc Level columns
all_IncLevel1 <- Reduce(function(x, y) merge(x, y), IncLevel1)
rownames(all_IncLevel1) <- all_IncLevel1[,1]
all_IncLevel1 <- all_IncLevel1[,-1]

all_IncLevel2 <- Reduce(function(x,y) merge(x, y), IncLevel2)
rownames(all_IncLevel2) <- all_IncLevel2[,1]
all_IncLevel2 <- all_IncLevel2[,-1]

all_IncLevel <- merge(all_IncLevel1, all_IncLevel2, by="row.names")
dim(all_IncLevel)



###Combine IJC and SJC to get a read coverage for a splicing event


rownames(all_IJC) <- all_IJC[,1]
all_IJC <- all_IJC[,-1]
all_IJC <- as.matrix(all_IJC)

rownames(all_SJC) <- all_SJC[,1]
all_SJC <- all_SJC[,-1]
all_SJC <- as.matrix(all_SJC)

dim(all_IJC)
dim(all_SJC)

sum <- all_IJC + all_SJC


###How many splicing events have less than 10 reads at a given splicing event?
####AKA occurences of a read count between 0-9

sum$count <- apply(sum, 1, function(x) sum(x %in% c(0:9)))

hist(sum$count)
sum$count



###Remove reads that have greater than ten percent of samples with read counts 10 or lower

ninetypercent <- 140

rownames(all_IncLevel) <- all_IncLevel[,1]
all_IncLevel <- all_IncLevel[,-1]
all_IncLevel <- as.matrix(all_IncLevel)

all_IncLevel <- cbind(sum$count, all_IncLevel)
all_IncLevel <- data.frame(all_IncLevel)
all_IncLevel <- all_IncLevel[which(all_IncLevel$V1<=ninetypercent),]
all_IncLevel <- all_IncLevel[,-1]


###Remove reads that have greater than 10 percent of samples samples with read counts 10 or lower for SJC or IJC


all_SJC <- cbind(sum$count, all_SJC)
all_SJC <- data.frame(all_SJC)
all_SJC <- all_SJC[which(all_SJC$V1<=ninetypercent),]
all_SJC <- all_SJC[,-1]

all_IJC <- cbind(sum$count, all_IJC)
all_IJC <- data.frame(all_IJC)
all_IJC <- all_IJC[which(all_IJC$V1<=ninetypercent),]
all_IJC <- all_IJC[,-1]



###Create dataframe of coord ID, skip length and inclusion length, filter for low counts

rMATS_lengths <- lapply(rMATS_2, function(x) subset(x, select=c("coord","IncFormLen","SkipFormLen")))

incl_lengths <- Reduce(function(x, y) merge(x, y), rMATS_lengths)

rownames(incl_lengths) <- incl_lengths[,1]
incl_lengths <- incl_lengths[,-1]
incl_lengths <- as.matrix(incl_lengths)

incl_lengths <- cbind(sum$count, incl_lengths)
incl_lengths <- data.frame(incl_lengths)
incl_lengths <- incl_lengths[which(incl_lengths$V1<=ninetypercent),]
incl_lengths <- incl_lengths[,-1]


###Final table information

dim(all_IncLevel)
dim(all_SJC)
dim(all_IJC)
dim(incl_lengths)

write.table(all_IncLevel, "/Volumes/padgetrlab/courtney/MLL/all/combinedrmatstable/A3SS_filteredlowcount.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)

write.table(all_SJC, "/Volumes/padgetrlab/courtney/MLL/all/combinedrmatstable/A3SS_filteredlowcount_SJC.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)

write.table(all_IJC, "/Volumes/padgetrlab/courtney/MLL/all/combinedrmatstable/A3SS_filteredlowcount_IJC.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)

write.table(incl_lengths, "/Volumes/padgetrlab/courtney/MLL/all/combinedrmatstable/A3SS_inclengths.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)
