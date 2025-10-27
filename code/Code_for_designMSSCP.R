# REQUIRES PACKAGE Matrix
# library(Matrix)

# IF RUNNING PARALLEL, LOAD PACKAGE doFuture
# library(doFuture)

# SHOULD CHANGE THE PATH BELOW ACCORDINGLY
# SET PATH
# setwd("C:/CrossedPaper/Review1/CODE")

# LOAD FUNCTIONS
# source("FunctionsMSPE.R", encoding = 'UTF-8')

# ORDER STRATUM BEING DESIGNED (CURRENT STRATUM)
Ns <- 1   

# NUMBER OF FACTORS IN CURRENT STRATUM              
K <- 2   

# NUMBER OF LEVELS OF EACH FACTOR IN THIS STRATUM          
Nlev <- c(3, 3)

# MODEL FORMULA
model <- formula(~ X1+X2+I(X1^2)+I(X2^2)+X1:X2)

# NUMBER OF MODEL PARAMETERS
Npar <- length(attr(terms(model), "term.labels")) 

# WEIGHTS FOR PARAMETERS IN THIS STRATUM (NOT USED FOR "D" CRITERION)
# 1 FOR LINEAR AND INTERACTIONS, 1/4 FOR QUADRATIC
# WEIGTHS IN THE SAME ORDER AS MODEL TERMS IN THE FORMULA
W <- rep(1,Npar)
W[(K+1):(2*K)] <- 1/4
W <- matrix(W/sum(W),nr=1) 

# FACTORS NAMES
indprefix <- 1:(K)
prefix <- "X"
Fnames <- paste(prefix, indprefix, sep="")

# DESIGN REGION: FOR CUBIC ENTER 'Y'       
Cubic <- 'Y'

# NUMBER OF UNITS (THE DESIGN IN THIS STRATUM IS COMPLETELY RANDOMIZED 
# IN 10 MAIN PLOTS)
N <- 10        

# COMPOUND CRITERIA
KapD <- 0                            # weight for D
KapDP <- 1/3                          # weight for DP
KapLP <- 0                            # weight for CIs (AP)
KapL <- 1/3                          # weight for point estimation (A)
KapDF <- 1/3                          # weight for df efficiency (lof)
Kappa <- c(KapD,KapDP,KapLP,KapL,KapDF)

# CONFIDENCE COEFFICIENT FOR DP
# (1-alpha)
probDP <- 0.95 

# CONFIDENCE COEFFICIENT FOR LP
# (1-alpha)
probLP <- 0.95 

# MULTIPLE COMPARISONS CORRECTION OF ALPHA?                                      
MC <- 'N'      # 'Y'=yes or 'N'=no  correction for multiple comparisons?

# SEARCHING A DESIGN STRATUM-BY-STRATUM

################################################################################
################################################################################
# POINT EXCHANGE ALGORITHM

# FORM LIST OF LEVELS OF THE K FACTORS (COMPLETE UNTIL 1:Nlev[K])

# SET candExternal <- NULL FOR THE FULL FACTORIAL AS CANDIDATE SET, OTHERWISE, CREAT 
# A MATRIX NAMED candExternal WITH THE REQUIRED TREATMENT SET
candExternal <- NULL

# TO RUN THE SEARCH
# NUMBER OF 'TRIES' FOR THE POINT EXCHANGE ALGORITHM
Ntries <- 100  
seed <- 5
set.seed(seed)

parallelTries <- TRUE         # FOR PARALLEL RUNNING
update.info <- TRUE           # FOR SHOWING PROGRESS OF THE SEARCH
plan(multisession)            # REQUIRED FOR PARALLEL RUNNING
# SEARCH DESIGN
CPphase1 <- SearchTreat(N, Npar, W, Kappa, Ntries, K, Nlev, Levels=NULL, 
                        candExternal=NULL, model, Fnames, probDP, probLP, MC,
                        parallelTries, update.info)
plan(sequential)              # REQUIRED FOR PARALLEL RUNNING

# save(CPphase1, file="CPphase1.RData")

################################################################################
# DESIGNING THE SECOND PHASE
# ORDER STRATUM BEING DESIGNED (CURRENT STRATUM)
Ns <- 2   

# NUMBER OF FACTORS IN ALL PREVIOUS STRATA (add all)             
K1 <- 2

# NUMBER OF FACTORS IN CURRENT STRATUM              
K2 <- 3   

# NUMBER OF LEVELS OF EACH FACTOR IN THIS STRATUM          
Nlev <- c(3, 3, 3)       

# MODEL FORMULA
model <- formula(~ X3+X4+X5+I(X3^2)+I(X4^2)+I(X5^2)+X1:X3+X1:X4+X1:X5+
                   X2:X3+X2:X4+X2:X5+X3:X4+X3:X5+X4:X5)

Npar <- length(attr(terms(model), "term.labels")) 

# WEIGHTS FOR PARAMETERS IN THIS STRATUM (NOT USED FOR "D" CRITERION)
# 1 FOR LINEAR AND INTERACTIONS, 1/4 FOR QUADRATIC
# WEIGTHS IN THE SAME ORDER AS MODEL TERMS IN THE FORMULA
W <- rep(1,Npar)

W[(K2+1):(2*K2)] <- 1/4
W <- matrix(W/sum(W),nr=1) 

# FACTORS NAMES
indprefix <- 1:(K1+K2)
prefix <- "X"
Fnames <- paste(prefix, indprefix, sep="")

# NUMBER OF 'BLOCKS' (UNITS IN PREVIOUS STRATUM, HERE IS 10 OVENS TIMES 3 BATCHS)
Nbloc <- 30     

# SIZE OF 'BLOCKS'       
Bsize <- 2

# nUMBER OF UNITS
N <- Nbloc*Bsize             

# THE DESIGN IN THE PREVIOUS STRATA 
# ENTER OR LOAD Xwhole: THE MATRIX OF TREATMENTS IN THE PREVIOUS STRATA
# NUMBER OF COLUMNS = K1
# NUMBER OF ROWS = Nbloc

# load("CPphase1.RData")  # IF NOT AVAILABLE IN WORKSPACE

# Because each row of previous design must be replicated more than
# Bsize (actually there are Nblocp "superblocks"), the row are pre-replicated
# here
Nblocp <- 3
Xwhole <- matrix(apply(CPphase1$Xopt,1,Expand,Nblocp),nc=K1,byrow=TRUE)

# COMPOUND CRITERIA
KapD <- 0                            # weight for D
KapDP <- 1/3                          # weight for DP
KapLP <- 0                            # weight for CIs (AP)
KapL <- 1/3                          # weight for point estimation (A)
KapDF <- 1/3                          # weight for df efficiency (lof)
Kappa <- c(KapD,KapDP,KapLP,KapL,KapDF)

# CONFIDENCE COEFFICIENT FOR DP
# (1-alpha)
probDP <- 0.95

# CONFIDENCE COEFFICIENT FOR LP
# (1-alpha)
probLP <- 0.95

# MULTIPLE COMPARISONS CORRECTION OF ALPHA?                                      
MC <- 'N'      # 'Y'=yes or 'N'=no  correction for multiple comparisons?

# SEARCHING A DESIGN STRATUM-BY-STRATUM

# SET candExternal <- NULL FOR THE FULL FACTORIAL AS CANDIDATE SET
# Construct candExternal otherwise

# TO RUN THE SEARCH
# POINT EXCHANGE
# NUMBER OF 'TRIES' FOR THE POINT EXCHANGE ALGORITHM
Ntries <- 100

algorithm <- "Point"  # or anything else for Coordinate (not good for PE criterion)
Levels <- NULL
inverse.updating <- "Yes"

seed <- 1738467814
set.seed(seed)

parallelTries <- TRUE         # FOR PARALLEL RUNNING
update.info <- TRUE           # FOR SHOWING PROGRESS OF THE SEARCH
plan(multisession)
CPphase2 <- SS_Search(algorithm, inverse.updating, Kappa, Xwhole,Nbloc,Bsize,
                         K1,K2,W,Npar,Nlev,Levels=NULL, Ntries,
                         probDP, probLP, MC,
                         model, Fnames,  candExternal=NULL, treat.restriction=NULL, 
                         parallelTries=TRUE,
                         update.info = TRUE)
plan(sequential)

# save(CPphase2, file="CPphase2.RData")

# THIRD PHASE - EX3
# INTERCHANGE GROUP UNITS 
set.seed(1738467814)
NCol <- 3               # Number of "superblocks" in previous strata                
Colsize <- 20           # Superblock size in the previous stratum  
NRow <- 10
Rowsize <- 6
PLOTsize <- 2
K1 <- 2                 # Number of factors in previous stratum VARYING IN SUPERBLOCKS
K2 <- 3

KapD <- 0               # weight for D
KapDP <- 1/3            # weight for DP
KapLP <- 0              # weight for CIs (AP)
KapL <- 1/3             # weight for point estimation (A)
KapDF <- 1/3            # weight for df efficiency (lof)   
Kappa <- c(KapD, KapDP, KapLP, KapL, KapDF)

# CONFIDENCE COEFFICIENT FOR DP
# (1-alpha)
probDP <- 0.95

# CONFIDENCE COEFFICIENT FOR LP
# (1-alpha)
probLP <- 0.95

# MULTIPLE COMPARISONS CORRECTION OF ALPHA?                                      
MC <- 'N'      # 'Y'=yes or 'N'=no  correction for multiple comparisons?

# FACTORS NAMES
indprefix <- 1:(K1+K2)
prefix <- "X"
Fnames <- paste(prefix, indprefix, sep="")

# TERMS IN THE MODEL TO ESTIMATE FROM SUPERBLOCKS (INDICATORS)                                                                          #
# INDICATORS OF MODEL TERMS WHEN THE DESIGN IS SEEN AS            
# BLOCKED BY THE SUPERBLOCKS                                                          
model <- formula(~ X3+X4+X5+I(X3^2)+I(X4^2)+I(X5^2)+
                   X1:X3+X1:X4+X1:X5+
                   X2:X3+X2:X4+X2:X5+
                   X3:X4+X3:X5+
                   X4:X5)

Npar <- length(attr(terms(model), "term.labels"))                                                          #

# WEIGHTS FOR PARAMETERS                                #
# 1 FOR LINEAR AND INTERACTIONS, 1/4 FOR QUADRATIC  

W <- rep(1,Npar)
W[(K2+1):(2*(K2))] <- .25
W <- W/sum(W)

# load("CPphase2.RData")  # IF NOT AVAILABLE IN THE WORKSPACE

Xfull <- CPphase2$Xopt
N <- nrow(Xfull)
Batches <- rep(rep(1:NCol,rep(PLOTsize,NCol)),NRow)
Ovens <- rep(1:NRow, rep(Rowsize,NRow))

XF <- cbind(Ovens, Batches, Xfull)
Xfull <- XF[order(XF[,2], XF[,1]),]
N <- nrow(Xfull)
row.names(Xfull) <- 1:N
MSS_CP <- FinalAdjust(Xfull, N, NRow, Rowsize, NCol, Colsize, PLOTsize, K1, K2, model, 
                           Npar, W, Fnames, Kappa, probDP, probLP)
colnames(MSS_CP$Xopt) <- c("Ovens", "Batches", Fnames)

save(MSS_CP, file="MSS_CP.RData")

