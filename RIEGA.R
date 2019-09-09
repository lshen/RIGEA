
####################################################################################
#   Regional Imaging Genetic Enrichment Analysis
#   Mining regional imaging genetic associations via voxelwise enrichment analysis
# --------------------------------------------------------------------
#   Input:
#   - p.voxel: voxel-wise SNP-iQT association p values
#   - parameters: p_threshold: for predefining significant voxel-wise SNP-iQT assciation
#     
#                 
# 
#   Output:
#   - p.enrich: enrichment p values for each brain ROI
#
# ------------------------------------------
#   # Author: Xiaohui Yao, Xiaohui.Yao@pennmedicine.upenn.edu #
#   Date created: Sep-01-2019
#   Date updated: --2019
#   @University of Pennsylvania Perelman School of Medicine
####################################################################################

ROI.p.RIGEA <- function(p.voxel, threshold)
{
  nROI = length(p.voxel) # number of ROIs
  nVoxel <- length(unlist(p.voxel)) # number of voxels across brain

  ### hypergeometric distribution function phyper(x, m, n, k), where: 
  #   - x: # of significant voxel-SNP associations located in the ROI; 
  #   - m: # of voxels in of the ROI
  #   - n: # of voxels out of the ROI
  #   - k: # of total significant voxel-SNP associations across the brain
  para.RIGEA <- matrix(NA, nrow = nROI, ncol = 4) # for each ROI, save parameters of phyper function phyper(x, m_, n ,k)
  colnames(para.RIGEA) <- c('hx', 'hm', 'hn', 'hk')
  
  para.RIGEA[, 1] <- sapply(p.voxel, function(x) length(which(x < threshold))) # x
  para.RIGEA[, 2] <- sapply(p.voxel, length) # m
  para.RIGEA[, 3] <- nVoxel - para.RIGEA[, 2] # n
  para.RIGEA[, 4] <- sum(para.RIGEA[, 1]) # k
  
  p.RIGEA <- phyper(para.RIGEA[, 1], para.RIGEA[, 2], para.RIGEA[, 3], para.RIGEA[, 4], lower.tail = F)
  RIGEA.result <- data.frame(ROI = paste0('ROI_', formatC(1:nROI, flag = '0', width = 3)), para.RIGEA, p.RIGEA = p.RIGEA)

  return(RIGEA.result)
}








