#R script for clustering on the cluster!

setwd("/home/hershbc/lustre/clustering_03_20_19/GitHub/scripts/")

#Get date
date <- Sys.Date()
format(date, format="%Y_%m_%d")
date
getwd()

#libraries
library(data.table)
library(stringr)
library(dplyr)
library(pheatmap)
library(matrixStats)
library(parallel)

###Read in final mutations file
#somatic mutations
binary <- read.csv("../inputfiles/filtered_tidied/2019_03_01_combined_mutations_deletions_controls_binary.txt", header = TRUE, sep="\t", stringsAsFactors = FALSE)

mutations <- read.csv("../inputfiles/filtered_tidied/2019_03_01_combined_mutations_deletions_controls.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)


#select for splicing factor mutations
binary_sf <- binary[,grepl("SF3B1|SRSF2|U2AF1",colnames(binary))]
#clean up splicing factors
binary_sf$SF3B1 <- ifelse(binary_sf$SF3B1_Hotspot==1 | binary_sf$SF3B1_missense==1, "1","0")
binary_sf$SF3B1_Hotspot <- NULL
binary_sf$SF3B1_missense <- NULL
binary_sf$SRSF2_missense <- NULL
binary_sf$U2AF1_missense <- NULL 

#select for non-splicing factor mutations
#binary_nosf <- binary[,!grepl("SF3B1|SRSF2|U2AF1|LUC7L2|ZRSR2|PRPF8|DDX41",colnames(binary))]
#binary_nosf <- binary_nosf[,which(colSums(binary_nosf)>=20)]

##Add expression information to binary file
###Create binary file for Expression of LUC7L2, DDX41, ZRSR2, PRPF8
expression <- fread("../inputfiles/filtered_tidied/2019_03_01_expression.txt",  sep = "\t", stringsAsFactors = FALSE, data.table=FALSE)

rownames(expression) <- expression[,1]
expression$V1 <- NULL
expression <- data.frame(t(expression))
expression <- expression[which(!rownames(expression) %in% mutations[which(mutations$mkh=="healthy_bm"),"array_id"]),]


expression_short <- expression[,c("LUC7L2","PRPF8","ZRSR2","DDX41"),]

#Getting lists for final table
expression_short <- expression_short[order(expression_short$LUC7L2),]
LUC7L2low <- rownames(expression_short)[1:140]

expression_short <- expression_short[order(expression_short$PRPF8),]
PRPF8low <- rownames(expression_short)[1:140]

expression_short <- expression_short[order(expression_short$ZRSR2),]
ZRSR2low <- rownames(expression_short)[1:140]

expression_short <- expression_short[order(expression_short$DDX41),]
DDX41low <- rownames(expression_short)[1:140]

binary_sf$LUC7L2_low <- ifelse(rownames(binary_sf)%in%LUC7L2low, "1","0")
binary_sf$DDX41_low <- ifelse(rownames(binary_sf)%in%DDX41low,"1","0")
binary_sf$PRPF8_low <- ifelse(rownames(binary_sf)%in%PRPF8low,"1","0")
binary_sf$ZRSR2_low <- ifelse(rownames(binary_sf)%in%ZRSR2low,"1","0")

#binary <- cbind(binary_sf,binary_nosf)
binary <- binary_sf

###Preparing binary file of Mutations
k <- which(binary>=1, arr.ind=TRUE)
binary[k] <- "T"

binary[binary == 0] <- "F"

###Read in rMATS output
#Get inclusion levels
inclusionlevels <- fread("../inputfiles/filtered_tidied/2019_03_01_filtered_inclusionlevels.txt", sep="\t", stringsAsFactors = FALSE, data.table=FALSE)
rownames(inclusionlevels) <- inclusionlevels$V1
inclusionlevels$V1 <- NULL

#get AS Event type
inclengths <- fread("../inputfiles/filtered_tidied/2019_03_01_filtered_inclengths.txt", sep="\t", stringsAsFactors = FALSE, data.table=FALSE)
rownames(inclengths) <- inclengths$V1
inclengths$V1 <- NULL
AS_type <- inclengths[,"type",drop=FALSE]
inclengths <- NULL

###List for disease subtype
AML <- mutations[which(mutations$mkh=="AML"),"array_id"]
MDS <- mutations[which(mutations$mkh=="MDS"),"array_id"]
CMML <- mutations[which(mutations$mkh=="CMML"),"array_id"]
MDS_MPN_RS_T <- mutations[which(mutations$mkh=="MDS/MPN-RS-T"),"array_id"]
MDS_MPN_U <- mutations[which(mutations$mkh=="MDS/MPN-U"),"array_id"]
healthybm <- mutations[which(mutations$mkh=="healthy_bm"),"array_id"]

###Make sample annotations
sample_annotations <- mutations[,c("array_id","mkh")]
sample_annotations <- unique(sample_annotations)
rownames(sample_annotations) <- sample_annotations$array_id
sample_annotations$array_id <- NULL

dat <- sample_annotations
dat$mkh <- as.factor(dat$mkh)

yay <- model.matrix(~.+0,dat)

sample_annotations <- as.data.frame(yay)
colnames(sample_annotations) <- sub("mkh", "", colnames(sample_annotations))

sample_annotations[sample_annotations==0] <- "F"
sample_annotations[sample_annotations==1] <- "T"

###Combine diangosis and mutations information and re-order
healthybm <- sample_annotations[,"healthy_bm",drop=FALSE]
sample_annotations$healthy_bm <- NULL

sample_annotations <- cbind(healthybm, sample_annotations, binary)

#Make annotation colors
make_annotcolors <- function(column){
 
  if(grepl("SF3B1|SRSF2|U2AF1|LUC7L2|ZRSR2|PRPF8|DDX41",column)){
    x <- list(column=c(T="black",F="white"))
  }
  else if(grepl("healthy_bm",column)){
    x <- list(column=c(T="#117733",F="white"))
  } 
  
  
  else if(grepl("MDS|AML|CMML|MDS/MPN-RS-T|MDS/MPN-U",column)){
    x <- list(column=c(T="#882255",F="white"))
  }
  
  else{
   x <- list(column=c(T="dark gray",F="white"))
  }

  return(x)
}

x <- sapply(colnames(sample_annotations), make_annotcolors)

names(x) <- colnames(sample_annotations)
annotation_colors <- x

#Rank by range
group <- cbind(AS_type, inclusionlevels)

group$values <- rowMaxs(as.matrix(group[,-1]),na.rm = TRUE) - rowMins(as.matrix(group[,-1]),na.rm=TRUE)

group<- group[order(-group$values),] 
group$values <- NULL

#separate rMATS output by AS type
group$type <- as.character(group$type)

A3SS <- group[which(group$type=="A3SS"),]
A3SS$type <- NULL
A5SS <- group[which(group$type=="A5SS"),]
A5SS$type <- NULL
SE <- group[which(group$type=="SE"),]
SE$type <- NULL
RI <- group[which(group$type=="RI"),]
RI$type <- NULL
MXE <- group[which(group$type=="MXE"),]
MXE$type <- NULL

rmats_subsets <- list(A3SS, A5SS, SE, RI, MXE)
rmats_labels <- list("A3SS", "A5SS", "SE", "RI", "MXE")

###make heatmap
makeheatmap <- function(rmats1, name){
   if(nrow(rmats1)>20000){  
	rmats1 <- rmats1[1:20000,]
	}
  file <- paste(date,"_heatmap_", name, "_20.jpg", sep="")
  print(file)
  jpeg(file, width = 12, height = 12, units = 'in', res = 300)
  pheatmap(rmats1, annotation_col = sample_annotations,annotation_colors = annotation_colors,show_rownames = FALSE, show_colnames = FALSE, main = name)
  dev.off()
  }

mcmapply(makeheatmap, rmats_subsets, rmats_labels)