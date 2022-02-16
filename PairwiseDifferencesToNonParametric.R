#!/usr/bin/env Rscript

####################################################################################
################################## Before to start #################################
####################################################################################

# skip lines related to installation of libraries because there are supposed to be already installed
skip_instalation <- scan(what="character", quiet = TRUE)
# install libraries
install.packages("benchmarkme")
install.packages("data.table")
install.packages("spaa")
install.packages("scales")

rm(skip_instalation)

# clean environment
rm(list=ls())
# reset graph devices
graphics.off()

# load packages avoiding warning messages
suppressPackageStartupMessages(library(benchmarkme))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(spaa))
suppressPackageStartupMessages(library(scales))

# get arguments
args = commandArgs(trailingOnly=TRUE)

# test if there are two arguments: if not, return an error
if (length(args)!=2) {
  stop("Please, provide a first cvs file of profiles and a second csv file of types\n
       USAGE: Rscript PairwiseDifferencesToNonParametric.R Profiles.csv Types.csv", call.=FALSE)
}

####################################################################################
############################### Profiles to dataframe ##############################
####################################################################################

# prepare profiles
## S stands for sample (i.e. rows): n = 20
## L stands for locus (i.e. columns): n = 30
## G stands for genotype (i.e. data): n= 600
## read dataframe of profiles (i.e. Profiles.csv) with missing data (encoded 0)
dfpm <- read.table(args[1], dec = ".", header=TRUE, sep = ",", quote = "")
## make sure that each variable of the dataframe is a character
dfpm <- data.frame(lapply(dfpm, as.character))
## replace missing data encoded 0 with NA
dfpm[ dfpm == "0" ] <- NA
## remove columns harboring NAs (i.e. locus with missing data)
dfp <- dfpm[ , colSums(is.na(dfpm)) == 0]

# transform profiles into matrix of pairwise differences
## identify logical CPUs available of your machine
CPU <- get_cpu()$no_of_cores
## set to use logical CPUs available minus 2 CPU for OS
setDTthreads(threads = CPU-2, restore_after_fork = TRUE, throttle = 1024)
## transpose dataframe
tdfp <- transpose(dfp, keep.names = "locus", make.names = "sample")

## calculate pairwise differences between all vectors (i.e. independently of the sample amount)
### multiple for-loops adding results in an empty vector
v <- integer() # create empty vector v
for (i in tdfp[, 2:ncol(tdfp)]) { # first for loop from the column 2 of the dataframe
  for (j in tdfp[, 2:ncol(tdfp)]){ # second for loop from the column 2 of the dataframe
    #print(i) # for i checking
    #print(j) # for j checking
    output <- sum(!i == j) # count the number of FALSE (i.e. reverse (!) of paired vectors i versus j)
    v <- c(v, output) # over right the output vector into the empty vector v
  }
}
### create a matrix adding the vector v by row independently of the sample amount
mat <- matrix(v, nrow = (ncol(tdfp)-1) , ncol = (ncol(tdfp)-1), byrow = TRUE)
### add row and column names of the matrix
#### keep the names of sample identifiers into a vector
samples = dfp$sample
#### add sample names to matrix rows
rownames(mat) <- samples
#### add sample names to matrix columns
colnames(mat) <- samples

### export the matrix into a csv file
write.csv(mat,file="PairwiseDifferencesMatrix.csv", row.names=TRUE, quote=FALSE)

# transform the matrix of pairwise differences into distobject and dataframe
## transform the matrix into a dist object
distobj = as.dist(mat)
## transform the dist object into a long dataframe of pairwise differences
dfl = dist2list(distobj)
## rename all variable
colnames(dfl) <- c("FirstSample","SecondSample","Differences")

# remove pairwise differences from identical samples (i.e. matrix diagonal = 0)
dfl$Diagonal <- ifelse((dfl$FirstSample == dfl$SecondSample), "YES", "NO")
dfl <- subset(dfl,dfl$Diagonal %in% c("NO"))
dfl$Diagonal <- NULL

# remove duplicates of pairwise differences
dfl[1:2] <- t( apply(dfl[1:2], 1, sort) )
dfl <- dfl[!duplicated(dfl[1:2]),]

# prepare dataframes of sample type
## PC stands for positive controls
## TS stands for tested samples
## read the dataframe of samples types
dft <- read.table(args[2], dec = ".", header=TRUE, sep = ",", quote = "")
## subset PC
dfPC <- subset(dft,dft$Type %in% c("PC"))
## subset TS
dfTS <- subset(dft,dft$Type %in% c("TS"))

# derive variables "FirstSample" and "SecondSample" into variables "FirstType" and "SecondType" flagging outbreak controls (OC), positive controls (PC), and negative controls (NC)
## derivation of the "FirstSample" column
dfl$FirstType <- ifelse(dfl$FirstSample %in% dfPC$Sample, "PC",
                        ifelse(dfl$FirstSample %in% dfTS$Sample, "TS",
                               "error"))
## derivation of the "SecondSample" column
dfl$SecondType <- ifelse(dfl$SecondSample %in% dfPC$Sample, "PC",
                         ifelse(dfl$SecondSample %in% dfTS$Sample, "TS",
                                "error"))

## flag the elements below
dfl$Status <- 
  ifelse((dfl$FirstType == "PC") & (dfl$SecondType == "PC"), "related",
         ifelse((dfl$FirstType == "PC") & (dfl$SecondType == "TS"), "tested",
                ifelse((dfl$FirstType == "TS") & (dfl$SecondType == "PC"), "tested",
                       ifelse((dfl$FirstType == "TS") & (dfl$SecondType == "TS"), "untested",
                              "error"))))

## export the dataframe into csv files
write.csv(dfl,file="PairwiseDifferencesDataframe.csv", row.names=FALSE, quote=FALSE)

####################################################################################
################################ Non-parametric tests ##############################
####################################################################################

# dataframe splitting
## positive controls (PC)
dflPC <- subset(dfl,dfl$Status %in% c("related"))
## tested samples (TS)
dflTS <- subset(dfl,dfl$Status %in% c("tested"))

# retrieve tested samples in a variable
dflTS$Tested <- 
  ifelse((dflTS$FirstType == "TS"), dflTS$FirstSample,
         ifelse((dflTS$SecondType == "TS"), dflTS$SecondSample,
                "error"))

# transform long dataframes into lists
## positive controls (PC)
lsTS <- split.data.frame(dflTS,dflTS$Tested)
lsTS <- lapply(lsTS,"[[","Differences")

# retrieve vector of differences of the outbreak controles (OC)
vPC <- dflPC$Differences

# calculate p-values
## Wilcoxon-Mann-Whitney tests (i.e. Wilcoxon rank-sum tests)
pWilcoxon <- integer() # create empty vector
for (i in lsTS) { # initiate loop
  output <- suppressWarnings((wilcox.test(vPC, i, paired = FALSE))$p.value) # perform test
  pWilcoxon <- c(pWilcoxon, output) # over right the output vector into the empty vector
}
## Kolmogorov-Smirnov tests
pKolmogorov <- integer() # create empty vector
for (i in lsTS) { # initiate loop
  output <- suppressWarnings((ks.test(vPC, i, paired = FALSE))$p.value) # perform test
  pKolmogorov <- c(pKolmogorov, output) # over right the output vector into the empty vector
}
## Kruskal-Wallis tests
pKruskal <- integer() # create empty vector
for (i in lsTS) { # initiate loop
  output <- suppressWarnings((kruskal.test(list(vPC, i)))$p.value) # perform test
  pKruskal <- c(pKruskal, output) # over right the output vector into the empty vector
}

# combine and export
## combine samples and p-values
TestedSample <- names(lsTS)
dfpvalues <- data.frame(TestedSample, pWilcoxon, pKolmogorov, pKruskal)
## reorder columns
dfpvalues <- dfpvalues[, c("TestedSample", "pKolmogorov", "pKruskal", "pWilcoxon")]
## control scientific format and digits
dfpvalues$pKolmogorov <- scientific(dfpvalues$pKolmogorov, digits = 10)
dfpvalues$pKruskal <- scientific(dfpvalues$pKruskal, digits = 10)
dfpvalues$pWilcoxon <- scientific(dfpvalues$pWilcoxon, digits = 10)
## force to keep scientific number
dfpvalues$pKolmogorov <- as.character(dfpvalues$pKolmogorov)
dfpvalues$pKruskal <- as.character(dfpvalues$pKruskal)
dfpvalues$pWilcoxon <- as.character(dfpvalues$pWilcoxon)

## export the dataframe into csv files
write.csv(dfpvalues,file="PairwiseDifferencesTests.csv", row.names=FALSE, quote=FALSE)

# add a message
cat("Developped by Nicolas Radomski on February 16 (2022) with the R version 4.1.2 (2021-11-01)","\n")
