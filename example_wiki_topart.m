data_name = 'wiki_TopArt.mat';
load(['data/wikipedia_networks/' data_name])

addpath(genpath('common/'))

K = 3; % number of communtities to fit
alpha = [0 1 10]; % try one of these; controls perturbation level in spectral clustering

%% Initialization using (matched) Bipartite Spectral Clustering
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


%% We have true labels, hence we compute the NMI w.r.t. true labels
% compute_mutual_info can computed NMI between any combination of
% (hard/soft, harf/soft) pair of input labels

% In this example we only have labels for side 2, so only compute the NMI on that side
nmi_sc = compute_mutual_info(l2_sc,l2);
nmi_mbisbm = compute_mutual_info(tau_2,l2);
fprintf('\n----- NMIs -----\n')
fprintf('%7s = %3.3f\n', 'SC', nmi_sc)
fprintf('%7s = %3.3f\n', 'mbiSBM', nmi_mbisbm)







    
