# Read a fasta file as a vector. You must have install bio3d
# install.packages ("bio3d")
library (bio3d)
fastaFile = read.fasta ("NC_001416-1.fasta")
vectorFasta = as.vector (fastaFile$ali)
