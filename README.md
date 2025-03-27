These are R codes to reproduce all the Modified Stratum-by-Stratum designs, for Example 3, presented in the manuscript. 

Packages Matrix, doFuture, progressr and rmarkdown are required as well as the functions file FunctionsMSPE.R provided. 

The experiment, a strip-split-plot layout, has four strata (Ovens*Batches/Runs): 1. Ovens, 2. Batches, 3. Ovens*Batches, 4. Runs. Level combinations of factors X1 and X2 are applied to Oven units (randomized design); no factors are applied to Batch units, and Ovens*Batches units and level combinations of factors X3 to X5 are applied to (Ovens*Batches/Runs). 

Following the algorithm presented in the paper, for each criterion, the design is built in three phases. The designs in phases 2 and 3 are conditional on the designs built in previous phases.

To reproduce the designs in Supplemental Table E (except D*, which is obtained from JMP, which can be read from EX3_JMP.CSV):
1) Save all files provided in your directory
2) Run the Main_code.R (this generates the designs)
3) Run Code_Table3_TableF_Figure6.RMD (this reproduces Table 3, Figure 6 and Table F).
