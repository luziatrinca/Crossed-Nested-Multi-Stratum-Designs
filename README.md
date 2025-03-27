These are files of the Supplemental Material.

Supp.pdf: describes the algorithms used to construct designs, gives the designs obtained for the examples in the paper and their efficiencies and gives a brief sensitivity study to explore the dependence on weights.

designs.zip: designs for all examples in CSV.

code.zip: contains R code to reproduce all the Modified Stratum-by-Stratum designs and their properties for Example 3.
Packages Matrix, doFuture, progressr and rmarkdown are required as well as the functions file FunctionsMSPE.R provided. 

To reproduce the designs in Supplemental Table E (except D*, which is obtained from JMP, which can be read from EX3_JMP.CSV):
1) Download code.zip to your computer and unzip it.
2) Run the Main_code.R (this generates the designs)
3) Run Code_Table3_TableF_Figure6.RMD (this reproduces Table 3, Figure 6 and Table F).
