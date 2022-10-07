# Epigenetics_Manuscript

This repository houses the analyses that appear in Tyler, Spruce, et al. (in prep): 
Variation in histone configurations correlates with gene expression across nine 
inbred strains of mice. 

The directories contain the following things:
* Code - mostly small functions that are called repeatedly.
* Documents - multiple directories containing different stages of the analysis and the manuscript.
  The analysis files include the following:
  
  * Run ChromHMM: 1.ChromHMM.Rmd
  ChromHMM (Ernst and Kellis, 2012) uses a hidden Markov model (HMM) 
  to identify states representing combinations of epigenetic marks. 
  For example, state 1 might be the co-occurrence of H3K4me1 and H3K4me3. 
  This step in the pipeline runs ChromHMM on the epigenetic data from all
  mouse strains across a series of states (4 to 2^n, where n is the number 
  of measured chromatin modifications.). This workflow starts with bam files, 
  generates bed files, binarizes the bed files using ChromHMM, and fits models 
  over a series of states. Results for each model are written to a separate 
  directory named for the number of states in the model. 

  * Prepare Gene Expression Data: 1.Expression.Rmd
  This workflow can be run in parallel to the above, hence the same number.
  It reads in gene expression data from all mouse strains and normalizes is 
  for further analysis.

  * Compare Models: 1.3.Compare_ChromHMM_States.Rmd
  This workflow combines the data from the previous step so we can look at all 
   ChromHMM models together. We select the model with the best corelation with 
   expression for follow-up analyses.
   
  * Analyze selected model: 1.4_Chromatin_States_and_Expression.Rmd
  This workflow is the primary analysis for the manuscript. It analyzes
  the selected model from the above step to investigate the relationship
  between chromatin states and gene expression.
  
  * Other analyses: 
  Other Rmarkdown files include additional analyses that were performed 
  for this study.
   
* Data and Results folders are not stored on GitHub.
* Data are available on NCBI Gene Expression Omnibus 
(GEO; https://www.ncbi.nlm.nih.gov/geo/) under accession 
number GSE213968.
* Results are generated locally as needed.



