function [tau_1,tau_2, theta, p, q, str] = ...
                fit_mbiSBM_v2(A, X, K, tau_1_init, tau_2_init, varargin)
% X is a cell array

p_temp = mean(A(:));
q_temp = 0.1*p_temp;

parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'VERBOSE',1)
addOptional(parser,'ignore_theta',0)
addOptional(parser,'first_tau_update',1)
%addOptional(parser,'p_init',0.1)
%addOptional(parser,'q_init',0.01)
addOptional(parser,'p_init', p_temp)
addOptional(parser,'q_init', q_temp)
addOptional(parser,'epsil',1e-3)
addOptional(parser,'model','poi')
addOptional(parser,'init_p_q_pr',false)

parse(parser, varargin{:});
VERBOSE = parser.Results.VERBOSE;
ignore_theta = parser.Results.ignore_theta;
first_tau_update = parser.Results.first_tau_update;
p_init = parser.Results.p_init;
q_init = parser.Results.q_init;
epsil = parser.Results.epsil;
model = parser.Results.model;
init_p_q_pr = parser.Results.init_p_q_pr;

if size(tau_1_init,2) == 1
    tau_1_init = label_vec2mat(tau_1_init,K);
end

if size(tau_2_init,2) == 1
    tau_2_init = label_vec2mat(tau_2_init,K);
end

N = size(A);
tau_1 = struct('val', tau_1_init, 'sum', sum(tau_1_init,1));
tau_2 = struct('val', tau_2_init, 'sum', sum(tau_2_init,1));


tau_norm = @(dtau) norm(dtau,'inf');
Tmax = 500;
%tol = 1e-4;
tol = .5*(1/K);
    
dim = [size(X{1},2) size(X{2},2)];
X_idx = { 1:dim(1); dim(1) + (1:dim(2)) };
D = sum(dim);
COVAR = D > 0;

%if COVAR 
Sig = eye(D);
mu = zeros(1,D);
sig2 = 10*[1,1];
Sigt = zeros(D,D,K);
for k=1:K
    Sigt(:,:,k) = eye(D);
end
mut = zeros(K,D);
%end

degs{1} = sum(A,2);
degs{2} = sum(A,1)';

%if VERBOSE, fprintf('%3s, %7s, %7s \n','t', 'delta_1', 'delta_2'), end
if VERBOSE > 0
    
    covar_str = 'unknown X'; % this shouldn't appear
    switch bin2dec(num2str(dim > 0))
        case 0, covar_str = 'no X';
        case 1, covar_str = 'X2';
        case 2, covar_str = 'X1';
        case 3, covar_str = 'X12';
    end
    
    %if COVAR, str = [str sprintf('X, ')]; else  str = [str sprintf('no X, ')]; end
    if ignore_theta, dc_str = sprintf('no DC'); else dc_str = sprintf('DC'); end
    if strcmpi(model,'poi')
        poi_model = 1;
        mod_str = sprintf('Poi');
    else
        poi_model = 0;
        mod_str = sprintf('Ber');
    end
    str = sprintf('mbiSBM: %5s, %5s, %3s', covar_str, dc_str, mod_str);
    fprintf('%30s',str)
    
end

% p = p_init;
% q = q_init;
% [pi_1, pi_2] = deal( ones(1,K)/K );

if ignore_theta
    theta{1} = ones(N(1),1); 
    theta{2} = ones(N(2),1);
end

%% Initialize p, q, pr by fixed values or based on initial labels
if init_p_q_pr
    p = p_init; 
    q = q_init; 
    [pi_1, pi_2] = deal(ones(1,K)/K);
    %fprintf('%f, %f \n',p,q)
else
    pi_1 = update_pi(tau_1);
    pi_2 = update_pi(tau_2);
    [p,q] = update_Q(A, tau_1, tau_2, N(1), N(2));
end

for t = 1:Tmax

    tau_1_old = tau_1;
    tau_2_old = tau_2;
    
    if ~ignore_theta
        theta{1} = dcbm_theta_update_spingarn(degs{1}, tau_1.val, varargin{:});
        theta{2} = dcbm_theta_update_spingarn(degs{2}, tau_2.val, varargin{:});
        %theta{1} = dcbm_theta_update_cham_pock(degs{1}, tau_1.val, varargin{:});
        %theta{2} = dcbm_theta_update_cham_pock(degs{2}, tau_2.val, varargin{:});
    end
    
    if poi_model  
        phi_1 = log(p+epsil) - log(q+epsil);
        phi_0 = q-p;
    else
        phi_1 = log(p*(1-q)+epsil) - log(q*(1-p)+epsil);
        phi_0 = log(1-p+epsil) - log(1-q+epsil);
    end
      
    [beta, betap] = update_beta(X, X_idx, N, K, Sigt, mut, sig2);
    
    if first_tau_update == 1
        tau_1 = update_tau(A ,  tau_2, theta{1}, beta{1}, pi_1, phi_1, phi_0, epsil, varargin{:});
        tau_2 = update_tau(A',  tau_1, theta{2}, beta{2}, pi_2, phi_1, phi_0, epsil, varargin{:});
    else
        tau_2 = update_tau(A',  tau_1, theta{2}, beta{2}, pi_2, phi_1, phi_0, epsil, varargin{:});
        tau_1 = update_tau(A ,  tau_2, theta{1}, beta{1}, pi_1, phi_1, phi_0, epsil, varargin{:});
    end
    
    delta_1 = tau_norm(tau_1.val - tau_1_old.val);
    delta_2 = tau_norm(tau_2.val - tau_2_old.val);
    
    if COVAR
        [Sigt, mut, Sig, mu, sig2] = update_params(X, dim, N, K, betap, tau_1, tau_2,  Sigt, mut, Sig, mu, sig2);
    end
    
    pi_1 = update_pi(tau_1);
    pi_2 = update_pi(tau_2);
    [p,q] = update_Q(A, tau_1, tau_2, N(1), N(2));


    if max([delta_1,delta_2]) < tol
       break
    end
end
%if VERBOSE, fprintf('Total # of itr = %3d\n',t), end
if VERBOSE > 0,  fprintf('%3d, %3.3e, %3.3e \n',t, delta_1, delta_2);  end
tau_1 = tau_1.val;
tau_2 = tau_2.val;
end


function out = row_normalize(X) 
    out = bsxfun(@times, X, 1./sum(X,2));
end

function out = safe_exp(X) 
    out = exp( bsxfun(@plus, X , -max(X,[],2)) );
end

function tau_out = update_tau(A, tau, theta, beta, pr, phi_1, phi_0, epsil, varargin)
    
    N = size(A,1);
    %tau = row_normalize( safe_exp(phi_1*A*t + phi_0*repmat(tsum,N,1) * diag(pr)) );
    
    a =  phi_1*A*tau.val + phi_0*repmat(tau.sum + log(pr),N,1) + beta;
    tval = dcbm_label_update_augL(a,theta, varargin{:});
    
    %tval = row_normalize( safe_exp( phi_1*A*tau.val + phi_0*repmat(tau.sum,N,1)+beta) * diag(pr) );
    
    tsum = max( sum(tval,1), epsil );
    %tsum = sum(tval,1);
    
    tau_out = struct('val', tval, 'sum',tsum);
end

function [p,q] = update_Q(A, tau_1, tau_2, N_1, N_2)
    temp = sum(tau_1.sum .* tau_2.sum);
  
    p = trace( tau_1.val' * (A * tau_2.val) ) / temp;
    r = N_1*N_2 / temp;
    q = (r*full(mean(A(:))) - p)/(r-1);
end

function pr = update_pi(tau)
    pr = tau.sum / sum(tau.sum);
end

function [beta, betap] = beta_func(X,St,mt,sig2)
   betap = trace(St) + sum( bsxfun(@minus, X, mt).^2,2);
   beta = - betap /(2*sig2);
end

function [beta, betap] = update_beta(X, idx, N, K, Sigt, mut, sig2)
    beta = {[];[]};
    betap = {[];[]};
    for r = 1:2
        beta{r} = zeros(N(r),K);
    
        if isempty(idx{r}), continue, end
    
        for k = 1:K
            [beta{r}(:,k), betap{r}(:,k)] = beta_func(X{r}, Sigt(idx{r},idx{r},k), mut(k,idx{r}), sig2(r));
        end
    end    
end

function [Sigt, mut, Sig, mu, sig2] = update_params(X, d, N, K, betap, tau_1, tau_2,  Sigt, mut, Sig, mu, sig2)
    tau{1} = tau_1.val;
    tau{2} = tau_2.val;
    tau_sum{1} = tau_1.sum;
    tau_sum{2} = tau_2.sum;

    Xch = {[],[]};
    much = {[],[]};
    for r = 1:2
        if d(r) > 0
            Xch{r} = tau{r}'*X{r}; % Dimension: K x d(r)
            much{r} = bsxfun(@times, Xch{r}, 1./(tau_sum{r}')); % Dimension: K x d(r)
        end
    end
    much = [much{:}]; % horizontal concatination

    %Sigt and mut update
    for k = 1:K
        D_k_inv = blkdiag( tau_sum{1}(k)/sig2(1) * eye(d(1)), tau_sum{2}(k)/sig2(2) * eye(d(2)));
        Sigt(:,:,k) = inv( D_k_inv + inv(Sig) );
        mut(k,:) = Sigt(:,:,k)*( D_k_inv*much(k,:)' + Sig\mu');
    end

    %mu update
    mu = mean(mut);

    %Sig update
    diff = bsxfun(@minus, mut, mu);
    Sig = mean(Sigt,3)+ diff'*diff/K;

    %sig2 update
    sig2 = zeros(2,1);
    for r = 1:2
        if d(r) > 0, sig2(r) = trace(tau{r}'*betap{r})/(N(r)*d(r)); end
    end
end

function out = flatness(pr)
    out = min(pr)/max(pr);
end
    