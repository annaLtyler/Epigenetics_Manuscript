This document describes individual steps in this analysis.

1. Run ChromHMM: 1.ChromHMM.Rmd
ChromHMM (Ernst and Kellis, 2012) uses a hidden Markov model (HMM) to identify states representing combinations of epigenetic marks. For example, state 1 might be the co-occurrence of H3K4me1 and H3K4me3. This step in the pipeline runs ChromHMM on the epigenetic data from all mouse strains across a series of states (4 to 2^n, where n is the number of measured chromatin modifications.). This workflow starts with bam files, generates bed files, binarizes the bed files using ChromHMM, and fits models over a series of states. Results for each model are written to a separate directory named for the number of states in the model. 

1. Prepare Gene Expression Data: 1.Expression.Rmd
This workflow can be run in parallel to the above, hence the same number.
It reads in gene expression data from all mouse strains and normalizes is for further analysis.

2. Correlate ChromHMM Models and Gene expression: 2.Chromatin_and_Expression.Rmd
This workflow correlates chromatin state with gene expression for each model from the previous step.


3. Compare Models: 3.Compare_ChromHMM_States.Rmd
This workflow combines the data from the previous step so we can look at all ChromHMM models together. We select the model with the best corelation with expression for follow-up analyses.

