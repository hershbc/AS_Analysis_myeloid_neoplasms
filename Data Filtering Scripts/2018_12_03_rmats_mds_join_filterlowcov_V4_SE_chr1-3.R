

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

setwd("/home/hershbc/lustre/testingrscripts/rmats_nostatt")

rmatsfiles <- list.files("/home/hershbc/lustre/testingrscripts/rmats_nostatt/", pattern = "SE.MATS.JC.txt", recursive=TRUE)
gtffiles <- list.files("/home/hershbc/lustre/testingrscripts/rmats_nostatt/", pattern = "fromGTF.SE.txt", recursive=TRUE)
files1 <- list.files("/home/hershbc/lustre/testingrscripts/rmats_nostatt/", pattern = "file1.txt", recursive=TRUE)
files2 <- list.files("/home/hershbc/lustre/testingrscripts/rmats_nostatt/", pattern = "file2.txt", recursive=TRUE)

print("read in files")

for(i in 1:length(rmatsfiles)){ 
  rMATS[[i]] <- read.csv(rmatsfiles[[i]], header=TRUE, sep="\t")
   
  rMATS_gtf[[i]] <- read.csv(gtffiles[[i]], header=TRUE, sep="\t")
  
  rMATS[[i]] <- merge(rMATS_gtf[[i]], rMATS[[i]], by = "ID")
  rMATS[[i]] <- rMATS[[i]][which(grepl("chr1$|chr2$|chr3", rMATS[[i]]$chr)),]  
}
#At this spot, the rMATS list contains all the files that need to be divided. 



###Unique Identifier for each splicing type


rMATS_2 <- rMATS
rm(rMATS)
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

namecolumn <- function(column, columnname1, columnname2, df){
  #If column name includes a 1, use list from file1, if not, use list from file2
  columnname <- ifelse(grepl("1",column), columnname1,columnname2)
  name <-strsplit(columnname, ",")
  name <- unlist(name)
  name <- gsub("(.+_.+?)(\\_.*)", "\\1", name)
  thing  <- df[,c("coord",column)]
  output <- separate(thing, column, name, sep = ",", remove=TRUE, convert=TRUE, extra = "warn", fill = "warn")
  return(output)
}

split_columns <- function(files) {
  setwd("/home/hershbc/lustre/testingrscripts/rmats_nostatt")
  #Condition file names
  columnname1 <- readChar(files$files1[[1]], file.info(files$files1[[1]])$size)
  #Control file names
  columnname2 <- readChar(files$files2[[1]], file.info(files$files2[[1]])$size)

  list_of_columns <- c("IJC_SAMPLE_1", "SJC_SAMPLE_1" , "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel1", "IncLevel2")
  split_list <- lapply(list_of_columns, namecolumn, columnname1 = columnname1, columnname2 = columnname2, df = files$rMATS_2)
  return(split_list)
}

files_df <- cbind(files1, files2, rMATS_2)
list_of_split_lists <- apply(files_df, 1, split_columns)

condition_IJC <- lapply(list_of_split_lists, "[[", 1)
condition_SJC <- lapply(list_of_split_lists, "[[", 2)
controls_IJC <- lapply(list_of_split_lists, "[[", 3)
controls_SJC <- lapply(list_of_split_lists, "[[", 4)
IncLevel1 <- lapply(list_of_split_lists, "[[", 5)
IncLevel2 <- lapply(list_of_split_lists, "[[", 6)

###Create dataframe of coord ID, skip length and inclusion length

rMATS_lengths <- lapply(rMATS_2, function(x) subset(x, select=c("coord","IncFormLen","SkipFormLen")))
rm(rMATS_2)
incl_lengths <- Reduce(function(x, y) merge(x, y), rMATS_lengths)
rm(rMATS_lengths)





###Join the dataframes by coordinate to get counts table for all splicing changes


all_IJC1 <- Reduce(function(x, y) merge(x, y), condition_IJC)
rownames(all_IJC1) <- all_IJC1[,1]
all_IJC1 <- all_IJC1[,-1]

all_IJC2 <- Reduce(function(x, y) merge(x, y), controls_IJC)
rownames(all_IJC2) <- all_IJC2[,1]
all_IJC2 <- all_IJC2[,-1]

all_IJC <- merge(all_IJC1, all_IJC2, by="row.names")
rm(all_IJC1)
rm(all_IJC2)

all_SJC1 <- Reduce(function(x, y) merge(x, y), condition_SJC)
rownames(all_SJC1) <- all_SJC1[,1]
all_SJC1 <- all_SJC1[,-1]

all_SJC2 <- Reduce(function(x, y) merge(x, y), controls_SJC)
rownames(all_SJC2) <- all_SJC2[,1]
all_SJC2 <- all_SJC2[,-1]

all_SJC <- merge(all_SJC1, all_SJC2, by="row.names")
rm(all_SJC1)
rm(all_SJC2)

dim(all_IJC)
dim(all_SJC)

print("dataframes joined, counts tables generated")

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

print("dataframes joined, inclusion levels generated")


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

print("counts joined and sum table created")
dim(sum)

###How many splicing events have less than 10 reads at a given splicing event?
####AKA occurences of a read count between 0-9

sum$count <- apply(sum, 1, function(x) sum(x %in% c(0:9)))

#hist(sum$count)




###Remove reads that have greater than ten percent of samples with read counts 10 or lower

ninetypercent <- 140

rownames(all_IncLevel) <- all_IncLevel[,1]
all_IncLevel <- all_IncLevel[,-1]
all_IncLevel <- as.matrix(all_IncLevel)

all_IncLevel <- cbind(sum$count, all_IncLevel)
all_IncLevel <- data.frame(all_IncLevel)
all_IncLevel <- all_IncLevel[which(all_IncLevel$V1<=ninetypercent),]
all_IncLevel <- all_IncLevel[,-1]
print("Inclusion levels filtered")

###Remove reads that have greater than 10 percent of samples samples with read counts 10 or lower for SJC or IJC

all_SJC <- cbind(sum$count, all_SJC)
all_SJC <- data.frame(all_SJC)
all_SJC <- all_SJC[which(all_SJC$V1<=ninetypercent),]
all_SJC <- all_SJC[,-1]

all_IJC <- cbind(sum$count, all_IJC)
all_IJC <- data.frame(all_IJC)
all_IJC <- all_IJC[which(all_IJC$V1<=ninetypercent),]
all_IJC <- all_IJC[,-1]

print("files filtered for low coverage")

###Filter skip length and inclusion length, filter for low counts

rownames(incl_lengths) <- incl_lengths[,1]
incl_lengths <- incl_lengths[,-1]
incl_lengths <- as.matrix(incl_lengths)

incl_lengths <- cbind(sum$count, incl_lengths)
incl_lengths <- data.frame(incl_lengths)
incl_lengths <- incl_lengths[which(incl_lengths$V1<=ninetypercent),]
incl_lengths <- incl_lengths[,-1]

print("lengths selected and filtered")

###Final table information

print("final dimensions")
dim(all_IncLevel)
dim(all_SJC)
dim(all_IJC)
dim(incl_lengths)

write.table(all_IncLevel, "/home/hershbc/lustre/testingrscripts/combinedrmatstable/SE_filteredlowcount_chr1-3.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)

write.table(all_SJC, "/home/hershbc/lustre/testingrscripts/combinedrmatstable/SE_filteredlowcount_SJC_chr1-3.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)

write.table(all_IJC, "/home/hershbc/lustre/testingrscripts/combinedrmatstable/SE_filteredlowcount_IJC_chr1-3.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)

write.table(incl_lengths, "/home/hershbc/lustre/testingrscripts/combinedrmatstable/SE_inclengths_chr1-3.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)

print("tables written")
