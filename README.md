These are files of the Supplemental Material of "Response Surface Designs for General Crossed and Nested Multi-Stratum Structures".

Supp.pdf: describes the algorithms used to construct designs, gives the formulae for DF calculations, gives a brief sensitivity study to explore the dependence on weights and gives the designs obtained for the examples in the paper and their efficiencies.

designs: designs for all examples in CSV.

code: contains R code to reproduce all the Modified Stratum-by-Stratum designs and their properties for Example 3.
Packages Matrix, doFuture, progressr, digest, and rmarkdown are required, as well as the functions file FunctionsMSPE.R provided. 

Following the algorithm presented in the paper, for each criterion, the design is built in three phases. The designs in phases 2 and 3 are conditional on the designs built in previous phases.

To reproduce the designs in Supplemental Table S5 (except D*, which is obtained from JMP, which can be read from EX3_JMP.CSV), run, sequentially, the files: 
1) Code_for_designMSSD.R; 
2) Code_for_designMSSDP.R and 
3) Code_for_designMSSCP.R. 
The code assumes that, for each criterion, the final design is stored externally in the same directory where the codes are stored.

To reproduce Table 3, Table S5 (Supplemental) and Figure 6, run the code Code_Table3_TableS5_Figure6.RMD (rmarkdown).

For examples 1, 2 and 4, if desired, Tables and Figures can be reproduced by running the respective RMD code, which reads the designs stored in the CSV files. 
