data_name = 'wiki_TopArt.mat';
load(['data/' data_name])
%%
%Xs{2} = zscore(Xs{2});

DC = 0;
diag_rest = 0;

fname = sprintf('%s_DC_%d_diag_%d', data_name, DC, diag_rest); 

N = size(A); 
K = 3;
alpha = [0 1 10]; % choose one of these

% random initialization with Bipartite Spectral Clustering
[l1_sc, l2_sc, Z_2] = biSpecClust2(A,K,1:K,'perturb',1,'pert_geom',1,'alpha',alpha(2));

[tau_1, tau_2, theta, Psi, opt_str, Sigt, mut, Sig, mu, sig2] =  ...
        fit_mbiSBM_v7(A, Xs, K, l1_sc, l2_sc, 'ignore_theta', ~DC, 'mode','poi','gen_sbm', 0, 'diag_rest', diag_rest);


%%
tic
[l1_sc, l2_sc, Z_2] = biSpecClust(A, ... % the adjacency matrix (n1 x n2) 
    K, ... % number of communities
    'perturb',true, ... % use perturbation
    'pert_geom', true, ... % use geometric scaling for the perturbation
    'alpha',alpha(2));
dt1 = toc;


%% mbiSBM, without degree correction (try the degree correction as well use 'ignore_theta', false)
tic
[tau_1, tau_2, theta, Psi, opt_str, Sigt, mut, Sig, mu, sig2] =  ...
        fit_mbiSBM(A, ... % the adjacency matrix (n1 x n2) 
            X, ... % covariates a cell array with two elements { [n1 x d1], [n2 x d2] } 
            K, ... % number of communities
            l1_sc, l2_sc, ... % initial labels for the two sides
            'ignore_theta', true, ... % ignore theta -- i.e., no degree correction
            'mode','poi', ... % use Poisson likelihood
            'gen_sbm', false, ... % do not use general SBM (using the simple p versus q model)
            'diag_rest', false); % do not use diagonal restriction for covariate covariances 
dt2 = toc;


%%
[~,estim_lab2] = max(tau_2'); 
%%
nmi = compute_mutual_info(estim_lab2,l2);



    