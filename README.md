# Usage
The R scripts PairwiseDifferences2NonParametric.R (detailed algorithm with Rstudio) and PairwiseDifferencesToNonParametric.R (automatic algorithm with Rscript) aim at deriving profiles of microbial mutations (e.g. cg/wgMLST, genes, SNPs, InDels, kmers) (Profiles.csv) into a matrix (PairwiseMatrix.csv) and a dataframe of pairwise differences (PairwiseDataframe.csv) with the objective to perform non-parametric tests Kolmogorov-Smirnov, Kruskal-Wallis and Wilcoxon-Mann-Whitney (PairwiseDifferencesTests.csv) to distinguish between samples related and unrelated (i.e. tested samples TS in Types.csv) to a studied outbreak (positive controls PC in Types.csv).
# Input
## Profiles of microbial mutations (i.e. Profiles.csv)
- S stands for sample
- L stands for locus
- G stands for genotype
- 0 stands for missing data
```
sample  L1  L2  L3  L4  L5  L6  L7  L8  L9 L10 L11 L12 L13 L14 L15 L16 L17 L18 L19 L20 L21 L22 L23 L24 L25 L26 L27 L28 L29 L30
    S1 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G5 G86 G54 G47  G1  G2  G3  G4  G5  G6   0 G91 G92 G93 G98  G3 G67 G12  G1
    S2 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90   0 G92 G93 G98  G3 G67 G12  G1
    S3 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G86 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91   0 G93 G98  G3 G67 G12  G1
    S4 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92   0 G98  G3 G67 G12  G1
    S5 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G98 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
    S6 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G5 G86 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
    S7 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
    S8 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G86 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
    S9 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
   S10 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G98 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
   S11 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
   S12 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G87 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
   S13 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G5 G86 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
   S14 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G86 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
   S15 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G98 G54 G47  G1  G2  G3  G4  G5  G6 G90 G91 G92 G93 G98  G3 G67 G12  G1
   S16 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G5 G86 G54 G47  G1  G2  G2  G5  G6  G7 G90 G91 G92 G93 G98  G3 G67 G12  G1
   S17 G20 G15 G55 G12 G30 G11 G24 G66 G12 G55 G66  G2 G86 G54 G47  G2  G1  G2  G5  G6  G7 G90 G91 G92 G93 G98  G3 G67 G12  G1
   S18 G20 G31 G55 G30 G30 G11 G55 G66 G55 G55 G66  G5 G98 G54 G47  G2  G1  G2  G5  G6  G7 G90 G91 G92 G93 G98  G3 G67 G12  G1
   S19 G20 G15 G55 G12 G30 G11 G24 G22 G12 G55 G66  G2 G87 G54 G47  G2  G1  G2  G5  G6  G7 G90 G91 G92 G93 G98  G3 G67 G12  G1
   S20 G20 G31 G55 G30 G30 G11 G55 G22 G55 G55 G66  G5 G87 G54 G47  G2  G1  G2  G5  G6  G7 G90 G91 G92 G93 G98  G3 G67 G12  G1
```
## Positive controls (PC) and tested samples (TS) (i.e. Types.csv)
```
Sample Type
    S1   PC
    S2   PC
    S3   PC
    S4   PC
    S5   PC
    S6   PC
    S7   PC
    S8   PC
    S9   PC
   S10   PC
   S11   TS
   S12   TS
   S13   TS
   S14   TS
   S15   TS
   S16   TS
   S17   TS
   S18   TS
   S19   TS
   S20   TS
```

# Output
## Matrix of pairwise differences of microbial mutations (i.e. PairwiseMatrix.csv)
```
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
```
## Dataframe of pairwise differences of microbial mutations (i.e. PairwiseDataframe.csv)
```
FirstSample SecondSample Differences FirstType SecondType   Status
         S1           S2           2        PC         PC  related
         S1           S3           1        PC         PC  related
         S1           S4           5        PC         PC  related
         S1           S5           5        PC         PC  related
         S1           S6           1        PC         PC  related
         S1           S7           3        PC         PC  related
         S1           S8           2        PC         PC  related
         S1           S9           6        PC         PC  related
         S1          S10           6        PC         PC  related
         S1          S11           2        PC         TS   tested
         S1          S12           5        PC         TS   tested
         S1          S13           1        PC         TS   tested
         S1          S14           2        PC         TS   tested
         S1          S15           6        PC         TS   tested
         S1          S16           4        PC         TS   tested
         S1          S17           7        PC         TS   tested
         S1          S18          11        PC         TS   tested
         S1          S19           9        PC         TS   tested
         S1          S20          12        PC         TS   tested
...
```
## Non-parametric tests (i.e. PairwiseDifferencesTests.csv)
```
TestedSample  pKolmogorov     pKruskal    pWilcoxon
         S11 9.986380e-01 4.761064e-01 4.830255e-01
         S12 9.999975e-01 8.235065e-01 8.321964e-01
         S13 8.989582e-01 2.654719e-01 2.702811e-01
         S14 9.582372e-01 5.249224e-01 5.322192e-01
         S15 9.999975e-01 8.583854e-01 8.671514e-01
         S16 5.666991e-02 3.405835e-03 3.529358e-03
         S17 1.330169e-04 6.196417e-06 6.528859e-06
         S18 1.330169e-04 4.580405e-06 4.829255e-06
         S19 1.330169e-04 4.586002e-06 4.835130e-06
         S20 1.330169e-04 4.580405e-06 4.829255e-06
```

# Install R and Rstudio
## 1/ Install R (from configured sources)
```
sudo apt update
sudo apt install r-base
R --version
```
## 2/ Install RStudio (from dowloaded rstudio-2021.09.1-372-amd64.deb)
```
sudo apt install gdebi-core
sudo gdebi /home/Downloads/rstudio-2021.09.1-372-amd64.deb
rstudio --version
```
# Update R version from 3.x to 4.x
## 1/ Check the current R version
```
R --version
```
## 2/ Add key and application repository
```
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
```
## 3/ Install the R version 4.x
```
sudo apt install r-base
```
## 4/ Check the current R version
```
R --version
```
# Update R 3.x packages to R 4.x
## 1/ Update installed packages from the R console
```
update.packages(ask = FALSE, checkBuilt = TRUE)
```
## 2/ Update missing packages from the R console
```
old_packages <- installed.packages(lib.loc = "/home/IZSNT/n.radomski/R/x86_64-pc-linux-gnu-library/3.6/")
new_packages <- installed.packages()
missing_packages <- as.data.frame(old_packages[!old_packages[, "Package"] %in% new_packages[, "Package"],])
install.packages(missing_packages$Package)
```
# Install R dependencies and launch with R
## 1/ Install dependencies from the R console
The R scripts PairwiseDifferences2NonParametric.R (detailed algorithm with Rstudio) and PairwiseDifferencesToNonParametric.R (automatic algorithm with Rscript) were prepared and tested with R version 4.1.2 and RStudio 2021.09.1.
```
install.packages("benchmarkme")
install.packages("data.table")
install.packages("spaa")
install.packages("scales")
```
## 2/ Launch each command from Rstudio (i.e. PairwiseDifferences2NonParametric.R detailed algorithm)
```
git clone https://github.com/Nicolas-Radomski/PairwiseDifferencesToNonParametric.git
cd PairwiseDifferencesToNonParametric
rstudio PairwiseDifferences2NonParametric.R
```
## 3/ Launch the whole script from Rscript (i.e. PairwiseDifferencesToNonParametric.R automatic algorithm)
```
git clone https://github.com/Nicolas-Radomski/PairwiseDifferencesToNonParametric.git
cd PairwiseDifferencesToNonParametric
Rscript PairwiseDifferencesToNonParametric.R Profiles.csv Types.csv
```
# Install Docker image and launch with Docker
## 1/ Pull Docker image from Docker Hub
```
docker pull nicolasradomski/pairwisedifferencestononparametric
```
## 2/ Launch with Docker and different paired-trees
```
docker run --name nicolas --rm -v /home/data:/data -v /home/output:/output nicolasradomski/pairwisedifferencestononparametric:latest sh -c 'Rscript code/PairwiseDifferencesToNonParametric.R data/Profiles.csv data/Types.csv' > output/std.log 2>&1
```
# Illustration
![Non-parametric figure](https://github.com/Nicolas-Radomski/PairwiseDifferencesToNonParametric/blob/main/illustration.png)
# References
- Jinlong Zhang J. Package 'spaa' - The Comprehensive R Archive Network. 2016, cran.r-project.org, Version 0.2.2
- First version: Radomski N., Cadel-Six S., Cherchame E., Felten A., Barbet P., Palma F., Mallet L., Le Hello S., Weill F.X., Guillier L. and Mistou M.Y. 1 A Simple and Robust Statistical Method to Define Genetic Relatedness of Samples Related to Outbreaks at the Genomic Scale - Application to Retrospective Salmonella Foodborne Outbreak Investigations. 2019, Frontiers in Microbiology, 10(2413): 1-13
# Acknowledgment
Adriano Di Pasquale for our discussions about algorithmic approaches
# Author
Nicolas Radomski
