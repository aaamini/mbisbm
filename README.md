# mbisbm


The code here implements algorithms for Matched BIpartite SBM (with node covariates/metadata) as
described in [this paper](https://arxiv.org/abs/1703.04943). 

The main file is `common/fit_mbiSBM.m` which implements the latest version of the algorithm.
This is the defulat version and should be used if you have a good initialization. 
See `example_wiki_cities.m` and `example_wiki_topart.m` for how to use this code.

The function `common/biSpecClust.m` implements the matched bipartite clustering described in the paper.
This can be used as the initialization for the mbisbm or as a standalone algorithm.

There is also `common/fit_mbiSBM_v2.m` which implements v.2 of the algorithm. This is used for example
in creating the some of the early figures in the paper. You can sun `simulations/compare_nmi_sims.m`
to see it in action. This version is more robust to the initial labeling 
(and implements less features that v.7). Use this if you have unreliable/bad initial labels.
