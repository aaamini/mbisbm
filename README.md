# mbisbm


The code here implements the Matched BIpartite SBM (with node covariates/metadata) algorithms 
described in ... 

The main fitting file is `common/fit_mbisbm.m' which implements v.7 of the algorithm.
This is the defulat version and should be used if you have good initialization. 
See `example_wiki_cities.m` and `example_wiki_topart.m` for how to use this code.

The function `common/biSpec.m' implements the matched bipartite clustering described in the paper.
This can be used as the initialization for the mbisbm or as a standalone algorithm.

There is also `common/fit_mbisbm_v2.m' which implements v.2 of the algorithm. This is used for example
in creating the some of the early figures in the paper. You can sun `simulations/compare_nmi_test.m'
to see it in action. This version is more robust to the initial labeling 
(and implements less features that v.7). Use this if you have unreliable/bad initial labels.
