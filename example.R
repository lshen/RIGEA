
###########################################################
#     Run this example to see how to implement RIGEA      #  
# ------------------------------------------------------- #
# Author: Xiaohui Yao, Xiaohui.Yao@pennmedicine.upenn.edu #
# @University of Pennsylvania Perelman School of Medicine #
########################################################### 

library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('RIEGA.R')

# load data
load('example.rda') # wm.example; y.example.classification; y.example.regression; gene

# set parameters, can try different threshold to check the enrichment result
# As an example, we fix it
threshold = 1e-5 # nominal significant SNP-voxel association

result.RIGEA <- ROI.p.RIGEA(p.voxel, threshold) %>% 
  arrange(p.RIGEA)
write.csv(result.RIGEA, file.path(dirname(rstudioapi::getActiveDocumentContext()$path), 'example_result_RIGEA.csv'), 
          row.names = F)
