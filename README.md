# mbisbm

`mbisbm` stands for Matched BIpartite Stochastic Block Model.
The code here implements algorithms in the paper entitled
[Matched bipartite block model with covariates](https://arxiv.org/abs/1703.04943). 

The main function is `common/fit_mbiSBM.m` which implements the latest version of the algorithm.
This is the defulat version and should be used if you have a good initialization. 
See `example_wiki_cities.m` and `example_wiki_topart.m` for how to use this function.

`common/biSpecClust.m` implements the matched bipartite clustering described in the paper.
This can be used as the initialization for the mbisbm or as a standalone algorithm.

There is also `common/fit_mbiSBM_v2.m` which implements v.2 of the algorithm. This is used for example
in creating some of the early figures in the paper. You can run `simulations/compare_nmi_sims.m`
to see it in action. This version is more robust to the initial labeling 
(and implements less features that v.7). Use it if you have unreliable/bad initial labels.
