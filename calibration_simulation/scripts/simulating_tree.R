### only once!
devtools::install_github("dosreislab/simclock")
### libraries
if(!require("diversitree")) install.packages("diversitree"); library("diversitree")
if(!require("geiger")) install.packages("geiger"); library("geiger")
if(!require("simclock")) install.packages("simclock"); library("simclock")

### output directory
dir_out = "calibration_simulation/1_trees/"

### simulate true tree
simulation = trees(
  pars = c(b= 0.1, 0.01), 
  type ="bd", 
  n = 1, 
  max.taxa = Inf, 
  max.t = 30,
  include.extinct = T
)
### picking tree
truetree = simulation[[1]]
### resolving politomies
truetree = ladderize(multi2di(truetree))
### plot
plot(truetree, cex= 0.4)
axisPhylo()
### export true tree
write.tree(truetree, paste0(dir_out,"truetree.nwk") )

### pruning out extinct lineages
livingtree = drop.extinct(truetree)
### plot
plot(livingtree, cex= 0.4)
axisPhylo()
### export true tree
write.tree(livingtree, paste0(dir_out,"livingtree.nwk") )
