### library
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### get tree
truetree = read.tree("calibration_simulation/1_trees/truetree.nwk")

### output directory
dir_out = "calibration_simulation/2_sequences/"

################################ relaxing tree  ################################

### plotting
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(truetree, cex= 0.4)
nodelabels(
  cex = 0.4
)
axisPhylo()

### tree to relax clock rates 
truetree_rel = truetree

### slow clade
clade1 = getDescendants(truetree_rel, node = 53)
truetree_rel$edge.length[truetree_rel$edge[, 2] %in% clade1] = truetree_rel$edge.length[truetree_rel$edge[, 2] %in% clade1] * 1.2

### fast clade
clade2 = getDescendants(truetree_rel, node = 36)
truetree_rel$edge.length[truetree_rel$edge[, 2] %in% clade2] = truetree_rel$edge.length[truetree_rel$edge[, 2] %in% clade2] *0.8

### plotting
plot(truetree_rel, cex= 0.4)
nodelabels(
  cex = 0.4
)

### export
write.tree(truetree_rel, paste0("calibration_simulation/1_trees/truetree_rel.nwk") )

############################## simulating locusA ##############################

### get tree
truetree_rel = read.tree("calibration_simulation/1_trees/truetree_rel.nwk")

### rate matrix
Q = matrix(
  c(0, 1, 1, 1,
    1, 0, 1, 1,
    1, 1, 0, 1,
    1, 1, 1, 0),
  nrow = 4, 
  byrow = TRUE,
  dimnames = list(c("a","c","g","t"), c("a","c","g","t"))
)

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
  x = binA,
  format = "fasta",
  file = paste0(dir_out, "locusA.fasta")
)

################################# simulating locusB ############################

### simulating locusB
locusB = simSeq(
  x = truetree_rel, 
  l = 1000, 
  Q = Q, 
  rate = 0.002,
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
  l = 250, 
  Q = Q, 
  rate = 0.001,
  bf = c(0.25,0.25,0.25,0.25),
  ancestral = F
)
locusCb = simSeq(
  x = truetree_rel, 
  l = 1000, 
  Q = Q, 
  rate = 0.002,
  bf = c(0.25,0.25,0.25,0.25),
  ancestral = F
)
locusCc = simSeq(
  x = truetree_rel, 
  l = 250, 
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

