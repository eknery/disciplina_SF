### load libraries
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")

### input diretory
dir_input = "calibration_simulation/2_sequences/"
file_names = list.files(dir_input)

### loading data
fasta_list = list()
for(i in 1:length(file_names) ){
  fasta_list[[i]] = read.phyDat(paste0(dir_input, file_names[i]),
                                format = "fasta",
                                type = "DNA"
  )
  names(fasta_list)[i] = gsub(".fasta", "", file_names[i])
}

############################## INFERING ML TREE ###############################

## select output dir
dir_out = "calibration_simulation/3_ml_inference/"

### loci names
all_loci = names(fasta_list)

### loop
for(locus_name in all_loci){
  ### pick one aligment
  one_fasta = fasta_list[[locus_name]]
  ### ML fits
  ml_fit = optim.pml(pml(NJ(dist.ml(one_fasta)), data = one_fasta), model = "JC")
  ### ML trees
  ml_tree = ml_fit$tree
  ## export trees
  write.tree(phy =  ml_tree,
             file = paste0(dir_out,locus_name,".nwk")
  )
  print(paste0("ML tree done: ", locus_name))
}

