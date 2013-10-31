#Set some data directories
tempdir="/Users/ogriffit/temp/"
outdir="/Users/ogriffit/Dropbox/BioStars/MachineLearning"

#install the core bioconductor packages, if not already installed
source("http://bioconductor.org/biocLite.R")
biocLite()

# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")
biocLite("affy")
biocLite("gcrma")
biocLite("org.Hs.eg.db")

# Install custom CDF packages
biocLite("hgu133ahsentrezgcdf")
biocLite("hgu133ahsentrezgprobe")
biocLite("hgu133ahsentrezg.db")
# Note, if the above do not install successfully with the biocLite method:
# Obtain current source files from:
# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
# Go to latest version, choose link to ENTREZG for appropriate version of R, find HGU95A row and download three source packages ('C', 'P' and 'A')
# In a terminal, run R CMD INSTALL <pckgfilename> where , <pckgfilename> is the name of the package file you downloaded

# Load all necessary libraries
library(GEOquery)
library(affy)
library(gcrma)
library(hgu133ahsentrezgcdf) #cdfname="HGU133A_HS_ENTREZG"
library(hgu133ahsentrezgprobe)
library(hgu133ahsentrezg.db)

# Set a data directory, download a GEO dataset, unpack and gunzip, and create a list of files for processing 
setwd(tempdir)
getGEOSuppFiles("GSE2990")
setwd(paste(tempdir,"GSE2990", sep=""))
untar("GSE2990_RAW.tar", exdir="data")
cels = list.files("data/", pattern = "cel")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "cel")

# Create AffyBatch object
setwd(paste(tempdir,"GSE2990/data", sep=""))
raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="HGU133A_HS_ENTREZG") 

# Perform GCRMA normalization
data.gcrma.norm=gcrma(raw.data)

# Extract expression values
gcrma=exprs(data.gcrma.norm)

# Remove Affy control probes - in this custom CDF, these are found in rows 12031 to 12098
# Check for yourself and look for probenames starting with "AFFX"
gcrma=gcrma[1:12030,] 

probes=row.names(gcrma)
symbol = unlist(mget(probes, hgu133ahsentrezgSYMBOL))
ID = unlist(mget(probes, hgu133ahsentrezgENTREZID))

#Combine gene annotations with raw data
gcrma=cbind(probes,ID,symbol,gcrma)

#Write GCRMA-normalized, mapped data to file
setwd(outdir)
write.table(gcrma, file = "trainset_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Get clinical details for this dataset
GSE2990_clindata=read.table("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2990/suppl/GSE2990_suppl_info.txt", header=TRUE) 
write.table(GSE2990_clindata, "trainset_clindetails.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
