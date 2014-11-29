# Load the library
library(biomaRt)

# Define biomart object
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Gives a list of all possible annotations; Currently there are 1668 listed
listAttributes(mart)

# Gives a list of all filters or criteria to search by; Currently there are 333 listed
# I chose to filter by: chromosome_name, start, end
listFilters(mart)

# Read in tab-delimited file with three columns: chromosome number, start position and end position
data = read.table("/Users/ogriffit/Desktop/test.txt")

# Extract HGNC gene symbol for a specific chr, start, stop range
#results <- getBM(attributes = c("hgnc_symbol", "chromosome_name","start_position", "end_position"), filters = c("chromosome_name","start", "end"), values = list(positions[,1], positions[,2],positions[,3]), mart = mart)

# Extract HGNC gene symbol, ENSG, etc for a specific chr, start, stop range for only KNOWN protein coding genes
#results <- getBM(attributes = c("hgnc_symbol", "chromosome_name","start_position", "end_position", "ensembl_gene_id", "gene_biotype", "status"), filters = c("chromosome_name","start", "end","biotype","status"), values = list(positions[,1], positions[,2], positions[,3], "protein_coding","KNOWN"), mart = mart)

getGeneSymbols = function(segment){
	chr=segment[1]
	start=segment[2]
	stop=segment[3]
	results <- getBM(attributes = c("hgnc_symbol", "chromosome_name","start_position", "end_position", "ensembl_gene_id", "gene_biotype", "status"), filters = c("chromosome_name","start", "end","biotype","status"), values = list(chr, start, stop, "protein_coding","KNOWN"), mart = mart)
	genes=results[,"ensembl_gene_id"]
	gene_string=paste(genes, collapse=",")
	return(gene_string)
}

#Could for loop through matrix, but for loops in R are very slow
#for (i in 1:3){
#	genes=getGeneSymbols(positions[i,])
#	print(genes)
#}

#Instead, use apply
results=apply(data,1,getGeneSymbols)
results_out=cbind(data,as.vector(results))
colnames(results_out)=c("chr","start","stop","ID","genes")

write.table(results_out, file="results.txt", quote=FALSE, row.names=FALSE, sep="\t")
