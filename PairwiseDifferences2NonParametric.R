####################################################################################
################################## Before to start #################################
####################################################################################

# install packages
install.packages("benchmarkme")
install.packages("data.table")
install.packages("spaa")
install.packages("scales")

# clean environment
rm(list=ls())
# reset graph devices
graphics.off()

# set working directory for Linux and Mac
setwd("/home/IZSNT/n.radomski/Documents/RstudioWorkingDirectory/PairwiseDifferencesToNonParametric-20220216")

# call library
library(benchmarkme)
library(data.table)
library(spaa)
library(scales)

####################################################################################
############################### Profiles to dataframe ##############################
####################################################################################

# prepare profiles
## S stands for sample (i.e. rows): n = 20
## L stands for locus (i.e. columns): n = 30
## G stands for genotype (i.e. data): n= 600
## read dataframe of profiles (i.e. Profiles.csv) with missing data (encoded 0)
dfpm <- read.table("Profiles.csv", dec = ".", header=TRUE, sep = ",", quote = "")
## make sure that each variable of the dataframe is a character
dfpm <- data.frame(lapply(dfpm, as.character))
## check nature of variables (must be character for each variable)
str(dfpm)
## check dimension (i.e. [1] 20 31)
dim(dfpm)
## check dataframe
dfpm
comment <- scan(what="character", quiet = TRUE)
   sample  L1  L2  L3  L4  L5  L6  L7  L8  L9 L10 L11 L12 L13 L14 L15 L16 L17 L18 L19 L20 L21 L22 L23 L24 L25 L26 L27 L28 L29 L30
1      S1 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G5 G86 G54 G47  G1  G2  G3  G4  G5  G6   0 G91 G92 G93 G98  G3 G67 G12  G1
2      S2 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90   0 G92 G93 G98  G3 G67 G12  G1
3      S3 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G86 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91   0 G93 G98  G3 G67 G12  G1
4      S4 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92   0 G98  G3 G67 G12  G1
5      S5 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G98 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
6      S6 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G5 G86 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
7      S7 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
8      S8 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G86 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
9      S9 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
10    S10 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G98 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
11    S11 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
12    S12 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
13    S13 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G5 G86 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
14    S14 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G86 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
15    S15 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G98 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
16    S16 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G5 G86 G54 G47  G1  G2  G2  G5  G6  G7 G90 G91 G92 G93 G98  G3 G67 G12  G1
17    S17 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G86 G54 G47  G2  G1  G2  G5  G6  G7 G90 G91 G92 G93 G98  G3 G67 G12  G1
18    S18 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G98 G54 G47  G2  G1  G2  G5  G6  G7 G90 G91 G92 G93 G98  G3 G67 G12  G1
19    S19 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G87 G54 G47  G2  G1  G2  G5  G6  G7 G90 G91 G92 G93 G98  G3 G67 G12  G1
20    S20 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G87 G54 G47  G2  G1  G2  G5  G6  G7 G90 G91 G92 G93 G98  G3 G67 G12  G1

rm(comment)

## replace missing data encoded 0 with NA
dfpm[ dfpm == "0" ] <- NA
## remove columns harboring NAs (i.e. locus with missing data)
dfp <- dfpm[ , colSums(is.na(dfpm)) == 0]
## check dimension (i.e. [1] 20 27)
dim(dfp)
## check dataframe
dfp
comment <- scan(what="character", quiet = TRUE)
   sample  L1  L2  L3  L4  L5  L6  L7  L8  L9 L10 L11 L12 L13 L14 L15 L16 L17 L18 L19 L20 L21 L26 L27 L28 L29 L30
1      S1 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G5 G86 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
2      S2 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G87 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
3      S3 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G86 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
4      S4 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G87 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
5      S5 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G98 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
6      S6 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G5 G86 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
7      S7 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G87 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
8      S8 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G86 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
9      S9 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G87 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
10    S10 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G98 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
11    S11 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G87 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
12    S12 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G87 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
13    S13 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G5 G86 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
14    S14 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G86 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
15    S15 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G98 G54 G47  G1  G2  G3  G4  G5  G6 G98  G3 G67 G12  G1
16    S16 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G5 G86 G54 G47  G1  G2  G2  G5  G6  G7 G98  G3 G67 G12  G1
17    S17 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G86 G54 G47  G2  G1  G2  G5  G6  G7 G98  G3 G67 G12  G1
18    S18 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G98 G54 G47  G2  G1  G2  G5  G6  G7 G98  G3 G67 G12  G1
19    S19 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G87 G54 G47  G2  G1  G2  G5  G6  G7 G98  G3 G67 G12  G1
20    S20 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G87 G54 G47  G2  G1  G2  G5  G6  G7 G98  G3 G67 G12  G1

rm(comment)

# transform profiles into matrix of pairwise differences
## identify logical CPUs available of your machine
CPU <- get_cpu()$no_of_cores
## check threads used by the function data.table
getDTthreads(verbose=TRUE)
## set to use logical CPUs available minus 2 CPU for OS
setDTthreads(threads = CPU-2, restore_after_fork = TRUE, throttle = 1024)
## check threads used by the function data.table
getDTthreads(verbose=TRUE)
## transpose dataframe
tdfp <- transpose(dfp, keep.names = "locus", make.names = "sample")
## check dimension (i.e. [1] 26 21)
dim(tdfp)
## check transposed dataframe
tdfp
comment <- scan(what="character", quiet = TRUE)
   locus  S1  S2  S3  S4  S5  S6  S7  S8  S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20
1     L1 G20 G20 G20 G20 G20 G20 G20 G20 G20 G20 G20 G20 G20 G20 G20 G20 G20 G20 G20 G20
2     L2 G15 G15 G15 G31 G31 G15 G15 G15 G31 G31 G15 G31 G15 G15 G31 G15 G15 G31 G15 G31
3     L3 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55
4     L4 G12 G12 G12 G30 G30 G12 G12 G12 G30 G30 G12 G30 G12 G12 G30 G12 G12 G30 G12 G30
5     L5 G30 G30 G30 G30 G30 G30 G30 G30 G30 G30 G30 G30 G30 G30 G30 G30 G30 G30 G30 G30
6     L6 G11 G11 G11 G11 G11 G11 G11 G11 G11 G11 G11 G11 G11 G11 G11 G11 G11 G11 G11 G11
7     L7 G24 G24 G24 G55 G55 G24 G24 G24 G55 G55 G24 G55 G24 G24 G55 G24 G24 G55 G24 G55
8     L8 G66 G66 G66 G66 G66 G22 G22 G22 G22 G22 G66 G66 G22 G22 G22 G66 G66 G66 G22 G22
9     L9 G12 G12 G12 G55 G55 G12 G12 G12 G55 G55 G12 G55 G12 G12 G55 G12 G12 G55 G12 G55
10   L10 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55 G55
11   L11 G66 G66 G66 G66 G66 G66 G66 G66 G66 G66 G66 G66 G66 G66 G66 G66 G66 G66 G66 G66
12   L12  G5  G2  G2  G5  G5  G5  G2  G2  G5  G5  G2  G5  G5  G2  G5  G5  G2  G5  G2  G5
13   L13 G86 G87 G86 G87 G98 G86 G87 G86 G87 G98 G87 G87 G86 G86 G98 G86 G86 G98 G87 G87
14   L14 G54 G54 G54 G54 G54 G54 G54 G54 G54 G54 G54 G54 G54 G54 G54 G54 G54 G54 G54 G54
15   L15 G47 G47 G47 G47 G47 G47 G47 G47 G47 G47 G47 G47 G47 G47 G47 G47 G47 G47 G47 G47
16   L16  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G2  G2  G2  G2
17   L17  G2  G2  G2  G2  G2  G2  G2  G2  G2  G2  G2  G2  G2  G2  G2  G2  G1  G1  G1  G1
18   L18  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G2  G2  G2  G2  G2
19   L19  G4  G4  G4  G4  G4  G4  G4  G4  G4  G4  G4  G4  G4  G4  G4  G5  G5  G5  G5  G5
20   L20  G5  G5  G5  G5  G5  G5  G5  G5  G5  G5  G5  G5  G5  G5  G5  G6  G6  G6  G6  G6
21   L21  G6  G6  G6  G6  G6  G6  G6  G6  G6  G6  G6  G6  G6  G6  G6  G7  G7  G7  G7  G7
22   L26 G98 G98 G98 G98 G98 G98 G98 G98 G98 G98 G98 G98 G98 G98 G98 G98 G98 G98 G98 G98
23   L27  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3  G3
24   L28 G67 G67 G67 G67 G67 G67 G67 G67 G67 G67 G67 G67 G67 G67 G67 G67 G67 G67 G67 G67
25   L29 G12 G12 G12 G12 G12 G12 G12 G12 G12 G12 G12 G12 G12 G12 G12 G12 G12 G12 G12 G12
26   L30  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1  G1

rm(comment)

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
### check the output vector v
print(v)
### check size of the vector v (i.e. [1] 400)
length(v)
### create a matrix adding the vector v by row independently of the sample amount
mat <- matrix(v, nrow = (ncol(tdfp)-1) , ncol = (ncol(tdfp)-1), byrow = TRUE)
### add row and column names of the matrix
#### keep the names of sample identifiers into a vector
samples = dfp$sample
#### add sample names to matrix rows
rownames(mat) <- samples
#### add sample names to matrix columns
colnames(mat) <- samples
### check class (i.e. [1] "matrix" "array")
class(mat)
### check matrix
mat
comment <- scan(what="character", quiet = TRUE)
    S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20
S1   0  2  1  5  5  1  3  2  6   6   2   5   1   2   6   4   7  11   9  12
S2   2  0  1  5  6  3  1  2  6   7   0   5   3   2   7   6   7  12   7  12
S3   1  1  0  6  6  2  2  1  7   7   1   6   2   1   7   5   6  12   8  13
S4   5  5  6  0  1  6  6  7  1   2   5   0   6   7   2   9  12   7  12   7
S5   5  6  6  1  0  6  7  7  2   1   6   1   6   7   1   9  12   6  13   8
S6   1  3  2  6  6  0  2  1  5   5   3   6   0   1   5   5   8  12   8  11
S7   3  1  2  6  7  2  0  1  5   6   1   6   2   1   6   7   8  13   6  11
S8   2  2  1  7  7  1  1  0  6   6   2   7   1   0   6   6   7  13   7  12
S9   6  6  7  1  2  5  5  6  0   1   6   1   5   6   1  10  13   8  11   6
S10  6  7  7  2  1  5  6  6  1   0   7   2   5   6   0  10  13   7  12   7
S11  2  0  1  5  6  3  1  2  6   7   0   5   3   2   7   6   7  12   7  12
S12  5  5  6  0  1  6  6  7  1   2   5   0   6   7   2   9  12   7  12   7
S13  1  3  2  6  6  0  2  1  5   5   3   6   0   1   5   5   8  12   8  11
S14  2  2  1  7  7  1  1  0  6   6   2   7   1   0   6   6   7  13   7  12
S15  6  7  7  2  1  5  6  6  1   0   7   2   5   6   0  10  13   7  12   7
S16  4  6  5  9  9  5  7  6 10  10   6   9   5   6  10   0   3   7   5   8
S17  7  7  6 12 12  8  8  7 13  13   7  12   8   7  13   3   0   6   2   7
S18 11 12 12  7  6 12 13 13  8   7  12   7  12  13   7   7   6   0   7   2
S19  9  7  8 12 13  8  6  7 11  12   7  12   8   7  12   5   2   7   0   5
S20 12 12 13  7  8 11 11 12  6   7  12   7  11  12   7   8   7   2   5   0

rm(comment)

### export the matrix into a csv file
write.csv(mat,file="PairwiseDifferencesMatrix.csv", row.names=TRUE, quote=FALSE)

# transform the matrix of pairwise differences into distobject and dataframe
## transform the matrix into a dist object
distobj = as.dist(mat)
## check class (i.e. [1] "dist")
class(distobj)
## check the dist object
distobj
comment <- scan(what="character", quiet = TRUE)
    S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19
S2   2                                                                
S3   1  1                                                             
S4   5  5  6                                                          
S5   5  6  6  1                                                       
S6   1  3  2  6  6                                                    
S7   3  1  2  6  7  2                                                 
S8   2  2  1  7  7  1  1                                              
S9   6  6  7  1  2  5  5  6                                           
S10  6  7  7  2  1  5  6  6  1                                        
S11  2  0  1  5  6  3  1  2  6   7                                    
S12  5  5  6  0  1  6  6  7  1   2   5                                
S13  1  3  2  6  6  0  2  1  5   5   3   6                            
S14  2  2  1  7  7  1  1  0  6   6   2   7   1                        
S15  6  7  7  2  1  5  6  6  1   0   7   2   5   6                    
S16  4  6  5  9  9  5  7  6 10  10   6   9   5   6  10                
S17  7  7  6 12 12  8  8  7 13  13   7  12   8   7  13   3            
S18 11 12 12  7  6 12 13 13  8   7  12   7  12  13   7   7   6        
S19  9  7  8 12 13  8  6  7 11  12   7  12   8   7  12   5   2   7    
S20 12 12 13  7  8 11 11 12  6   7  12   7  11  12   7   8   7   2   5

## transform the dist object into a long dataframe of pairwise differences
dfl = dist2list(distobj)
## check dimension (i.e. [1] 400   3)
dim(dfl)
## rename all variable
colnames(dfl) <- c("FirstSample","SecondSample","Differences")

# remove pairwise differences from identical samples (i.e. matrix diagonal = 0)
dfl$Diagonal <- ifelse((dfl$FirstSample == dfl$SecondSample), "YES", "NO")
dfl <- subset(dfl,dfl$Diagonal %in% c("NO"))
dfl$Diagonal <- NULL
dim(dfl) # [1] 380   3

# remove duplicates of pairwise differences
dfl[1:2] <- t( apply(dfl[1:2], 1, sort) )
dfl <- dfl[!duplicated(dfl[1:2]),]
dim(dfl) # [1] 190   3

# prepare dataframes of sample type
## PC stands for positive controls
## TS stands for tested samples
## read the dataframe of samples types
dft <- read.table("Types.csv", dec = ".", header=TRUE, sep = ",", quote = "")

### check dataframe
dft
comment <- scan(what="character", quiet = TRUE)
   Sample Type
1      S1   PC
2      S2   PC
3      S3   PC
4      S4   PC
5      S5   PC
6      S6   PC
7      S7   PC
8      S8   PC
9      S9   PC
10    S10   PC
11    S11   TS
12    S12   TS
13    S13   TS
14    S14   TS
15    S15   TS
16    S16   TS
17    S17   TS
18    S18   TS
19    S19   TS
20    S20   TS

rm(comment)

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

## check long dataframe
dfl
comment <- scan(what="character", quiet = TRUE)
    FirstSample SecondSample Differences FirstType SecondType   Status
2            S1           S2           2        PC         PC  related
3            S1           S3           1        PC         PC  related
4            S1           S4           5        PC         PC  related
5            S1           S5           5        PC         PC  related
6            S1           S6           1        PC         PC  related
7            S1           S7           3        PC         PC  related
8            S1           S8           2        PC         PC  related
9            S1           S9           6        PC         PC  related
10           S1          S10           6        PC         PC  related
11           S1          S11           2        PC         TS   tested
12           S1          S12           5        PC         TS   tested
13           S1          S13           1        PC         TS   tested
14           S1          S14           2        PC         TS   tested
15           S1          S15           6        PC         TS   tested
16           S1          S16           4        PC         TS   tested
17           S1          S17           7        PC         TS   tested
18           S1          S18          11        PC         TS   tested
19           S1          S19           9        PC         TS   tested
20           S1          S20          12        PC         TS   tested
23           S2           S3           1        PC         PC  related
24           S2           S4           5        PC         PC  related
25           S2           S5           6        PC         PC  related
26           S2           S6           3        PC         PC  related
27           S2           S7           1        PC         PC  related
28           S2           S8           2        PC         PC  related
29           S2           S9           6        PC         PC  related
30          S10           S2           7        PC         PC  related
31          S11           S2           0        TS         PC   tested
32          S12           S2           5        TS         PC   tested
33          S13           S2           3        TS         PC   tested
34          S14           S2           2        TS         PC   tested
35          S15           S2           7        TS         PC   tested
36          S16           S2           6        TS         PC   tested
37          S17           S2           7        TS         PC   tested
38          S18           S2          12        TS         PC   tested
39          S19           S2           7        TS         PC   tested
40           S2          S20          12        PC         TS   tested
44           S3           S4           6        PC         PC  related
45           S3           S5           6        PC         PC  related
46           S3           S6           2        PC         PC  related
47           S3           S7           2        PC         PC  related
48           S3           S8           1        PC         PC  related
49           S3           S9           7        PC         PC  related
50          S10           S3           7        PC         PC  related
51          S11           S3           1        TS         PC   tested
52          S12           S3           6        TS         PC   tested
53          S13           S3           2        TS         PC   tested
54          S14           S3           1        TS         PC   tested
55          S15           S3           7        TS         PC   tested
56          S16           S3           5        TS         PC   tested
57          S17           S3           6        TS         PC   tested
58          S18           S3          12        TS         PC   tested
59          S19           S3           8        TS         PC   tested
60          S20           S3          13        TS         PC   tested
65           S4           S5           1        PC         PC  related
66           S4           S6           6        PC         PC  related
67           S4           S7           6        PC         PC  related
68           S4           S8           7        PC         PC  related
69           S4           S9           1        PC         PC  related
70          S10           S4           2        PC         PC  related
71          S11           S4           5        TS         PC   tested
72          S12           S4           0        TS         PC   tested
73          S13           S4           6        TS         PC   tested
74          S14           S4           7        TS         PC   tested
75          S15           S4           2        TS         PC   tested
76          S16           S4           9        TS         PC   tested
77          S17           S4          12        TS         PC   tested
78          S18           S4           7        TS         PC   tested
79          S19           S4          12        TS         PC   tested
80          S20           S4           7        TS         PC   tested
86           S5           S6           6        PC         PC  related
87           S5           S7           7        PC         PC  related
88           S5           S8           7        PC         PC  related
89           S5           S9           2        PC         PC  related
90          S10           S5           1        PC         PC  related
91          S11           S5           6        TS         PC   tested
92          S12           S5           1        TS         PC   tested
93          S13           S5           6        TS         PC   tested
94          S14           S5           7        TS         PC   tested
95          S15           S5           1        TS         PC   tested
96          S16           S5           9        TS         PC   tested
97          S17           S5          12        TS         PC   tested
98          S18           S5           6        TS         PC   tested
99          S19           S5          13        TS         PC   tested
100         S20           S5           8        TS         PC   tested
107          S6           S7           2        PC         PC  related
108          S6           S8           1        PC         PC  related
109          S6           S9           5        PC         PC  related
110         S10           S6           5        PC         PC  related
111         S11           S6           3        TS         PC   tested
112         S12           S6           6        TS         PC   tested
113         S13           S6           0        TS         PC   tested
114         S14           S6           1        TS         PC   tested
115         S15           S6           5        TS         PC   tested
116         S16           S6           5        TS         PC   tested
117         S17           S6           8        TS         PC   tested
118         S18           S6          12        TS         PC   tested
119         S19           S6           8        TS         PC   tested
120         S20           S6          11        TS         PC   tested
128          S7           S8           1        PC         PC  related
129          S7           S9           5        PC         PC  related
130         S10           S7           6        PC         PC  related
131         S11           S7           1        TS         PC   tested
132         S12           S7           6        TS         PC   tested
133         S13           S7           2        TS         PC   tested
134         S14           S7           1        TS         PC   tested
135         S15           S7           6        TS         PC   tested
136         S16           S7           7        TS         PC   tested
137         S17           S7           8        TS         PC   tested
138         S18           S7          13        TS         PC   tested
139         S19           S7           6        TS         PC   tested
140         S20           S7          11        TS         PC   tested
149          S8           S9           6        PC         PC  related
150         S10           S8           6        PC         PC  related
151         S11           S8           2        TS         PC   tested
152         S12           S8           7        TS         PC   tested
153         S13           S8           1        TS         PC   tested
154         S14           S8           0        TS         PC   tested
155         S15           S8           6        TS         PC   tested
156         S16           S8           6        TS         PC   tested
157         S17           S8           7        TS         PC   tested
158         S18           S8          13        TS         PC   tested
159         S19           S8           7        TS         PC   tested
160         S20           S8          12        TS         PC   tested
170         S10           S9           1        PC         PC  related
171         S11           S9           6        TS         PC   tested
172         S12           S9           1        TS         PC   tested
173         S13           S9           5        TS         PC   tested
174         S14           S9           6        TS         PC   tested
175         S15           S9           1        TS         PC   tested
176         S16           S9          10        TS         PC   tested
177         S17           S9          13        TS         PC   tested
178         S18           S9           8        TS         PC   tested
179         S19           S9          11        TS         PC   tested
180         S20           S9           6        TS         PC   tested
191         S10          S11           7        PC         TS   tested
192         S10          S12           2        PC         TS   tested
193         S10          S13           5        PC         TS   tested
194         S10          S14           6        PC         TS   tested
195         S10          S15           0        PC         TS   tested
196         S10          S16          10        PC         TS   tested
197         S10          S17          13        PC         TS   tested
198         S10          S18           7        PC         TS   tested
199         S10          S19          12        PC         TS   tested
200         S10          S20           7        PC         TS   tested
212         S11          S12           5        TS         TS untested
213         S11          S13           3        TS         TS untested
214         S11          S14           2        TS         TS untested
215         S11          S15           7        TS         TS untested
216         S11          S16           6        TS         TS untested
217         S11          S17           7        TS         TS untested
218         S11          S18          12        TS         TS untested
219         S11          S19           7        TS         TS untested
220         S11          S20          12        TS         TS untested
233         S12          S13           6        TS         TS untested
234         S12          S14           7        TS         TS untested
235         S12          S15           2        TS         TS untested
236         S12          S16           9        TS         TS untested
237         S12          S17          12        TS         TS untested
238         S12          S18           7        TS         TS untested
239         S12          S19          12        TS         TS untested
240         S12          S20           7        TS         TS untested
254         S13          S14           1        TS         TS untested
255         S13          S15           5        TS         TS untested
256         S13          S16           5        TS         TS untested
257         S13          S17           8        TS         TS untested
[ reached 'max' / getOption("max.print") -- omitted 24 rows ]

rm(comment)

## export the dataframe into csv files
write.csv(dfl,file="PairwiseDifferencesDataframe.csv", row.names=FALSE, quote=FALSE)

####################################################################################
################################ Non-parametric test ###############################
####################################################################################

# dataframe splitting
## positive controls (PC)
dflPC <- subset(dfl,dfl$Status %in% c("related"))
dim (dflPC) # 45   6
## tested samples (TS)
dflTS <- subset(dfl,dfl$Status %in% c("tested"))
dim (dflTS) # 100   6

# retrieve tested samples in a variable
dflTS$Tested <- 
  ifelse((dflTS$FirstType == "TS"), dflTS$FirstSample,
         ifelse((dflTS$SecondType == "TS"), dflTS$SecondSample,
                "error"))

# transform long dataframes into lists
## positive controls (PC)
lsTS <- split.data.frame(dflTS,dflTS$Tested)
class(lsTS) # [1] "list"
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

## check outcomes
dfpvalues
comment <- scan(what="character", quiet = TRUE)
TestedSample  pKolmogorov     pKruskal    pWilcoxon
1           S11 9.986380e-01 4.761064e-01 4.830255e-01
2           S12 9.999975e-01 8.235065e-01 8.321964e-01
3           S13 8.989582e-01 2.654719e-01 2.702811e-01
4           S14 9.582372e-01 5.249224e-01 5.322192e-01
5           S15 9.999975e-01 8.583854e-01 8.671514e-01
6           S16 5.666991e-02 3.405835e-03 3.529358e-03
7           S17 1.330169e-04 6.196417e-06 6.528859e-06
8           S18 1.330169e-04 4.580405e-06 4.829255e-06
9           S19 1.330169e-04 4.586002e-06 4.835130e-06
10          S20 1.330169e-04 4.580405e-06 4.829255e-06

rm(comment)

## export the dataframe into csv files
write.csv(dfpvalues,file="PairwiseDifferencesTests.csv", row.names=FALSE, quote=FALSE)

# add a message
cat("Developped by Nicolas Radomski on February 16 (2022) with the R version 4.1.2 (2021-11-01)","\n")
