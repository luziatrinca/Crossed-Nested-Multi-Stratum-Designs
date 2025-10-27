# Reproducing Results for Example 3

# REQUIRES PACKAGE digest
# library(digest)

# REQUIRES PACKAGE Matrix
library(Matrix)

# REQUIRES doFuture
library(doFuture)

# set the path
setwd("C:/CrossedPaper/Review1/Review2/Designs_with_Code")

# LOAD FUNCTIONS
source("FunctionsMSPE.R", encoding = 'UTF-8')

# Run the search 
# Design MSS_D_S
source("Code_for_designMSSD.R")

# Design MSS_DP_S
source("Code_for_designMSSDP.R")

# Design MSS_CP_S
source("Code_for_designMSSCP.R")

# Now run, in Rmarkdown, the Code_Table3_TableF_Figure6.RMD to gererate the 
# results in the paper