#Set some data directories
tempdir="/Users/obigriffith/temp"
outdir="/Users/obigriffith/git/biostar-tutorials/MachineLearning"

#install the necessary bioconductor packages, if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
BiocManager::install("GEOquery")
BiocManager::install("affy")
BiocManager::install("gcrma")
BiocManager::install("org.Hs.eg.db")

# Install custom CDF packages
#BiocManager::install("hgu133ahsentrezgcdf")
#BiocManager::install("hgu133ahsentrezgprobe")
#BiocManager::install("hgu133ahsentrezg.db")

# Note, if the above do not install successfully you can alternately obtain current source files from:
# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
# Go to latest version, choose link to ENTREZG for appropriate version of R, find HGU133A row and copy download urls for three source packages ('C', 'P' and 'A')
# In a terminal, wget the urls and then run R CMD INSTALL <pckgfilename> where , <pckgfilename> is the name of the package file you downloaded
# wget http://mbni.org/customcdf/25.0.0/entrezg.download/hgu133ahsentrezgcdf_25.0.0.tar.gz
# wget http://mbni.org/customcdf/25.0.0/entrezg.download/hgu133ahsentrezgprobe_25.0.0.tar.gz
# wget http://mbni.org/customcdf/25.0.0/entrezg.download/hgu133ahsentrezg.db_25.0.0.tar.gz
# R CMD INSTALL hgu133ahsentrezgcdf_25.0.0.tar.gz 
# R CMD INSTALL hgu133ahsentrezgprobe_25.0.0.tar.gz
# R CMD INSTALL hgu133ahsentrezg.db_25.0.0.tar.gz

# Load all necessary libraries
library(GEOquery)
library(affy)
library(gcrma)
library(hgu133ahsentrezgcdf) #cdfname="HGU133A_HS_ENTREZG"
library(hgu133ahsentrezgprobe)
library(hgu133ahsentrezg.db)

# Set a data directory, download a GEO dataset, unpack and gunzip, and create a list of files for processing 
setwd(tempdir)
getGEOSuppFiles("GSE2034")
setwd(paste(tempdir,"GSE2034", sep="/"))
untar("GSE2034_RAW.tar", exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")

# Create AffyBatch object
setwd(paste(tempdir,"GSE2034/data", sep="/"))
raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="HGU133A_HS_ENTREZG") 

# Perform GCRMA normalization
data.gcrma.norm=gcrma(raw.data)

# Extract expression values
gcrma=exprs(data.gcrma.norm)

# Remove Affy control probes - look for probenames starting with "AFFX"
gcrma=gcrma[which(!grepl("AFFX", rownames(gcrma))),]

probes=row.names(gcrma)
symbol = unlist(mget(probes, hgu133ahsentrezgSYMBOL))
ID = unlist(mget(probes, hgu133ahsentrezgENTREZID))

#Combine gene annotations with raw data
gcrma=cbind(probes,ID,symbol,gcrma)

#Write GCRMA-normalized, mapped data to file
setwd(outdir)
write.table(gcrma, file = "trainset_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Get clinical details for this dataset
GSE2034_clindata=getGSEDataTables("GSE2034")[[2]][1:286,]
write.table(GSE2034_clindata, "trainset_clindetails.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
