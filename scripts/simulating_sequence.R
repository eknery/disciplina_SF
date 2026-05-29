### library
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("phylotools")) install.packages("phylotools"); library("phylotools")

### output directory
dir_out = "2_sequences/"

### plot nodes
plot(truetree, cex= 0.4)
edgelabels(
  text = truetree$edge[,2],
  cex = 0.4
)
### tree to relax clock rates 
truetree_rel = truetree
### exponential clock
scalers = rexp(n = nrow(truetree$edge), rate = 0.5)
### modify tree
truetree_rel$edge.length = truetree_rel$edge.length * scalers

### plotting
plot(truetree_rel, cex= 0.4)
write.tree(truetree_rel, paste0("1_tree/truetree_rel.nwk") )

### rate matrix
Q = matrix(
  c(0, 2, 1, 1,
    2, 0, 1, 1,
    1, 1, 0, 2,
    1, 1, 2, 0),
  nrow = 4, 
  byrow = TRUE,
  dimnames = list(c("A","T","C","G"), c("A","T","C","G"))
)

################################# simulating locusA ############################

### sequences for each block
locusA = simSeq(
  x = truetree_rel, 
  l = 500, 
  Q = Q, 
  rate = 0.001,
  bf = c(0.25,0.25,0.25,0.25),
  ancestral = F
)
## mark extinct
exA = grepl("ex", names(locusA))
## remove extinct
locusA = locusA[!exA]
## conver to binary
binA = as.DNAbin(locusA)
### export locus1
write.dna(
  x = bin1,
  format = "fasta",
  file = paste0(dir_out, "locusA.fasta")
)

################################# simulating locusB ############################

### simulating locusB
locusB = simSeq(
  x = truetree_rel, 
  l = 500, 
  Q = Q, 
  rate = 0.003,
  bf = c(0.25,0.25,0.25,0.25),
  ancestral = F
)
## mark extinct
exB = grepl("ex", names(locusB))
## remove extinct
locusB = locusB[!exB]
## conver to binary
binB = as.DNAbin(locusB)
### export locus
write.dna(
  x = binB,
  format = "fasta",
  file = paste0(dir_out, "locusB.fasta")
)

################################# simulating locusC ############################

### simulating locusC
locusCa = simSeq(
  x = truetree_rel, 
  l = 100, 
  Q = Q, 
  rate = 0.001,
  bf = c(0.25,0.25,0.25,0.25),
  ancestral = F
)
locusCb = simSeq(
  x = truetree_rel, 
  l = 300, 
  Q = Q, 
  rate = 0.003,
  bf = c(0.25,0.25,0.25,0.25),
  ancestral = F
)
locusCc = simSeq(
  x = truetree_rel, 
  l = 100, 
  Q = Q, 
  rate = 0.001,
  bf = c(0.25,0.25,0.25,0.25),
  ancestral = F
)
### join categories
locusC = c(locusCa,locusCb,locusCc)
## mark extinct
exC = grepl("ex", names(locusC))
## remove extinct
locusC = locusC[!exC]
## conver to binary
binC = as.DNAbin(locusC)
### export locus
write.dna(
  x = binC,
  format = "fasta",
  file = paste0(dir_out, "locusC.fasta")
)
