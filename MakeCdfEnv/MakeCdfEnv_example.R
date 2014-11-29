#install the core bioconductor packages, if not already installed
source("http://bioconductor.org/biocLite.R")
biocLite()

# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")
biocLite("affy")
biocLite("gcrma")
biocLite("hugene10stv1cdf")
biocLite("hugene10stv1probe")
biocLite("hugene10stprobeset.db")
biocLite("hugene10sttranscriptcluster.db")

#Load the necessary libraries
library(GEOquery)
library(affy)
library(gcrma)
library(hugene10stv1cdf)
library(hugene10stv1probe)
library(hugene10stprobeset.db)
library(hugene10sttranscriptcluster.db)

#Set working directory for download
setwd("/Users/ogriffit/Dropbox/BioStars")

#Download the CEL file package for this dataset (by GSE - Geo series id)
getGEOSuppFiles("GSE27447")

#Unpack the CEL files
setwd("/Users/ogriffit/Dropbox/BioStars/GSE27447")
untar("GSE27447_RAW.tar", exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")

setwd("/Users/ogriffit/Dropbox/BioStars/GSE27447/data")
raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="hugene10stv1") #From bioconductor

#perform RMA normalization (I would normally use GCRMA but it did not work with this chip)
data.rma.norm=rma(raw.data)

#Get the important stuff out of the data - the expression estimates for each array
rma=exprs(data.rma.norm)

#Format values to 5 decimal places
rma=format(rma, digits=5)

#Map probe sets to gene symbols or other annotations
#To see all available mappings for this platform
ls("package:hugene10stprobeset.db") #Annotations at the exon probeset level
ls("package:hugene10sttranscriptcluster.db") #Annotations at the transcript-cluster level (more gene-centric view)

#Extract probe ids, entrez symbols, and entrez ids
probes=row.names(rma)
Symbols = unlist(mget(probes, hugene10sttranscriptclusterSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hugene10sttranscriptclusterENTREZID, ifnotfound=NA))

#Combine gene annotations with raw data
rma=cbind(probes,Symbols,Entrez_IDs,rma)

#Write RMA-normalized, mapped data to file
write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



#Create a package for downloaded CDF
#This assumes you have already downloaded or created a CDF and you just want to use it in Bioconductor with (for example) ReadAffy
#For example suppose you wanted to use a specific CDF from the Affy web site:
#http://www.affymetrix.com/support/technical/byproduct.affx?product=hugene-1_0-st-v1

#Note, this example is purely for illustration. There is no real need to create a package for HuGene-1_0-st-v1
#This chip is already supported by BioConductor and can be loaded with 
#library(hugene10stv1cdf)   #cdfname="hugene10stv1"

#Install package for making cdf packages
source("http://bioconductor.org/biocLite.R")
biocLite("makecdfenv")
library(makecdfenv)

#Create CDF package in temporary directory 
pkgpath <- tempdir()
make.cdf.package("HuGene-1_0-st-v1.r3.cdf", cdf.path="/Users/ogriffit/Downloads/HuGene-1_0-st-v1.r3.unsupported-cdf", compress=FALSE, species = "Homo_sapiens", package.path = pkgpath)
dir(pkgpath)

#Install that package at a terminal using 'pkgpath' from above
#R CMD INSTALL /var/folders/8j/bqry255x52q6_dhyw22w6sq80000gn/T//Rtmppd9YEK/hugene10stv1.r3cdf

#Then, load it for use here
library(hugene10stv1.r3cdf)

#Download a CEL file package for testing purposes
getGEOSuppFiles("GSE27447")

#Unpack the CEL files
untar("GSE27447/GSE27447_RAW.tar", exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")

setwd("/Users/ogriffit/data")
raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="hugene10stv1.r3cdf") #Custom installed

#You can now go on to whatever normalizing and analysis you wish with the data using your custom CDF package
#perform RMA normalization 
data.rma.norm=rma(raw.data)

#Get the important stuff out of the data - the expression estimates for each array
rma=exprs(data.rma.norm)

