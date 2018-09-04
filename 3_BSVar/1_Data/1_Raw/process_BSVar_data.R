####################
#### TITLE:     Preprocess the raw data from the between-study heterogeneity assessment into an R object
#### Contents:
####
#### Source Files:
#### First Modified: 04/09/2018
#### Notes:
#################



##
###############
### Notes
###############
##

# Here we read in the t-values and save them to one **R** object.
# Note: we will read in the studies in random order for privacy reasons.


##
###############
### Preparation
###############
##

# Working directory should be source location

# Libraries
library(oro.nifti)
library(neurobase)

# Number of studies
nstud <- 33
database <- read.csv2('../../database.csv', header = TRUE, stringsAsFactors = FALSE)

# Read in ROI mask 
ROI <- readNIfTI2('../../pain_pAgF_z_FDR_0.01_forward_mask')[,,]
# Dimension in 3D
DIM3D <- dim(ROI)
# Switch to array
ROI <- array(ROI, dim = prod(DIM3D))
# MNI standard (2mm)
MNI <- readNIfTI2('../../MNI152_T1_2mm_brain')[,,]
# MNI mask
MNImask <- readNIfTI2('../../MNI152_T1_2mm_brain_mask')[,,]

# Scramble the sequence to read in studies
# -- done for privacy reasons.
scramStud <- sample(1:nstud, size = nstud, replace = FALSE)

# Small check
if(length(unique(scramStud)) != nstud){
  stop('Something went wrong with scrambling the read in sequence!')
}

##
###############
### Data processing
###############
##

# Empty data matrix
allStud <- matrix(NA, nrow = prod(DIM3D), ncol = nstud)

# For loop over the studies
for(i in scramStud){
  # Name of dataset is first column of database
  studDat <- readNIfTI2(paste(getwd(),'/', database[i,'img'], '.nii', sep = ''))[,,]
  # Check if dimensions match
  if(all(dim(studDat) != DIM3D)) stop(paste0('Dimension of image ', i,
                                             'does not match MNI space'))
  # Switch to vector
  studDat <- array(studDat, dim = prod(DIM3D))
  # Values outside MNI mask to NA
  studDat[MNImask == 0] <- NA
  # Bind to data matrix
  allStud[,i] <- studDat
  # Reset
  rm(studDat)
}

# Small checks
summary(allStud)
hist(allStud)


##
###############
### Write R object
###############
##

saveRDS(object = allStud, file = paste('../2_Processed/allStud.rds', sep = ''),
        compress = TRUE)




