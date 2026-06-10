### library
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### get estimate
mltree = read.tree("calibration_simulation/3_ml_inference/locusA.nwk")

### plot nodes
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(mltree, cex= 0.4)
nodelabels(
  cex = 0.4
)
### GENETIC DIVERGENCE ESTIMATE
## Clade (sp 15, sp16)
## age: 2.8 - 3.2 mya

### FOSSIL RECORD 
## All clades (MRCA) : 
## age: 25 - 30 mya

############################## Fitting clock models #######################

### calibration node
node = #
### age range
age = c("min" = 2.8, "max" = 3.2)

### Strict clock
chrono.strict = chronos(
  phy = mltree, 
  model = "discrete", 
  calibration = makeChronosCalib(
   phy = mltree, 
   node = node, 
   age.min = age["min"], 
   age.max = age["max"]
   ),
  control = chronos.control(
    nb.rate.cat = 1,
    dual.iter.max = 200
  )
)

### Relaxed uncorrelated
chrono.relunc = chronos(
  phy = mltree, 
  model = "relaxed", 
  calibration = makeChronosCalib(
    phy = mltree, 
    node = node, 
    age.min = age["min"], 
    age.max = age["max"]
  ),
  control = chronos.control(
    dual.iter.max = 200
  )
)

### Relaxed correlated
chrono.relcor = chronos(
  phy = mltree, 
  model = "correlated", 
  calibration = makeChronosCalib(
    phy = mltree, 
    node = node, 
    age.min = age["min"], 
    age.max = age["max"]
  ),
  control = chronos.control(
    dual.iter.max = 200
  )
)

############################# plotting chronograms ############################

### get true tree
livingtree = read.tree("calibration_simulation/1_trees/livingtree.nwk")

### visual parameters
par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))

### true tree
plot(
  ladderize(livingtree), 
  main = "Verdade", 
  cex = 0.6,
)
axisPhylo()

###
plot(
  ladderize(chrono.strict), 
  main = "Estrito", 
  cex = 0.6,
  )
axisPhylo()

###
plot(
  ladderize(chrono.relunc), 
  main = "Relaxado", 
  cex = 0.6,
  )
axisPhylo()

###
plot(
  ladderize(chrono.relcor), 
  main = "Relaxado correlacionado", 
  cex = 0.6,
  )
axisPhylo()

### total time
total_time = c(
  "Verdade" = max(branching.times(livingtree)),
  "Estrito" = max(branching.times(chrono.strict)),
  "Relaxado" = max(branching.times(chrono.relunc)),
  "Relaxado corr" = max(branching.times(chrono.relcor))
)
### 
total_time
