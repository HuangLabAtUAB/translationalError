# translationalError
The R codes for the translational error SAAV project.

Below is the annotation for the code files.

1. getMousePeptide.R - the R script for obtaining the theoretical mouse peptides.
2. backgroundDataFrames.R and savStatistics.R - to summarize the patterns for SAAVs from PDX or from TCGA somatic mutation.
3. CollectGS.R - collect gold-standard PSMs from mouseDb-informed search
4. evaFpOpen.R and fpOpenEvaFuc.R - parse and evalute MSFragger open search results.
5. dpSearch.R and fpEvaFun.R - parse and evalute MaxQuant dependent search results.
6. pFind_eva_all.R - parse and evaluate Open-pFind search results.
7. plotFixedSearch.R, variableSearchMqParse.R, and VariableModi.R - parse and evaluate fixed search results with predefine SAAVs.
8. ParseMaf.R – translate somatic mutations into customized databases.


Other files:
1. HOM_MouseHumanSequence.rpt - the data table containing all the human and mouse homologous gene information
2. aaDelta.rds, aaMass.rds, hmPairs.rds - the data frames for the data analysis

Citation:
[Song Y, Zhou, Q, and Huang C. J Proteome Res. 2025](https://pubmed.ncbi.nlm.nih.gov/41328725/)
