These are files of the Supplemental Material of "Response Surface Designs for General Crossed and Nested Multi-Stratum Structures".

Supp.pdf: describes the algorithms used to construct designs, gives the formulae for DF calculations, gives a brief sensitivity study to explore the dependence on weights and gives the designs obtained for the examples in the paper and their efficiencies.

designs: designs for all examples in CSV.

code: contains R code to reproduce all the Modified Stratum-by-Stratum designs and their properties for Example 3.
Packages Matrix, doFuture, progressr, digest, and rmarkdown are required, as well as the functions file FunctionsMSPE.R provided. 

To reproduce the designs in Supplemental Table S5 (except D*, which is obtained from JMP, which can be read from EX3_JMP.CSV):
1) Download all files in the code to your computer
2) Run the Main_code.R (this generates the designs)
3) Run Code_Ex3_Table3_TableS6_Figure7.RMD (this reproduces Table 3, Figure 7 and Table S6).
