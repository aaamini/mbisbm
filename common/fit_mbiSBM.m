function [tau_1,tau_2, theta, Psi, str, Sigt, mut, Sig, mu, sig2] = ...
                fit_mbiSBM(A, X, K, tau_1_init, tau_2_init, varargin)
% X is a cell array

parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'VERBOSE',true)
addOptional(parser,'gen_sbm',false)
addOptional(parser,'diag_rest',false)
addOptional(parser,'ignore_theta',false)
addOptional(parser,'first_tau_update',1)
addOptional(parser,'manual_init',false)
%addOptional(parser,'p_init',0.1)
%addOptional(parser,'q_init',0.01)
addOptional(parser,'epsil',1e-3)
addOptional(parser,'model','poi')
%addOptional(parser,'tau_init',{})
%addOptional(parser,'perturb_init','None')   % 'None', 'Dir', 'hard'
%addOptional(parser,'perturb_init_t',0.1)   % add 10\% noise
%addOptional(parser,'perturb_init_alpha',0.5) 

parse(parser, varargin{:});
VERBOSE = parser.Results.VERBOSE;
ignore_theta = parser.Results.ignore_theta;
first_tau_update = parser.Results.first_tau_update;
%p_init = parser.Results.p_init;
%q_init = parser.Results.q_init;
epsil = parser.Results.epsil;
model = parser.Results.model;
gen_sbm = parser.Results.gen_sbm;
diag_rest = parser.Results.diag_rest;
manual_init = parser.Results.manual_init;
% perturb_init = parser.Results.perturb_init;
% perturb_init_t = parser.Results.perturb_init_t;
% perturb_init_alpha = parser.Results.perturb_init_alpha;

if manual_init,
    % override gen_sbm option, we will use planted partition (pp) model
    gen_sbm = false;
end

% perturb_str = '';
% switch lower(perturb_init)
%     case 'dir'
%         [tau_1_init, tau_2_init] = dirchlet_perturb(tau_1_init, tau_2_init, perturb_init_t, perturb_init_alpha);
%         perturb_str = '~ dir';
%     case 'hard'
% end

%% The following function is now external
% random_init = @(N,K) mnrnd(1,ones(1,K)/K,N);  

N = size(A);

if isempty(tau_1_init) 
    tau_1_init = random_init(N(1),K);
end

if isempty(tau_2_init) 
    tau_2_init = random_init(N(2),K);
end


if size(tau_1_init,2) == 1
    tau_1_init = label_vec2mat(tau_1_init,K);
end

if size(tau_2_init,2) == 1
    tau_2_init = label_vec2mat(tau_2_init,K);
end


tau_1 = struct('val', tau_1_init, 'sum', sum(tau_1_init,1));
tau_2 = struct('val', tau_2_init, 'sum', sum(tau_2_init,1));


tau_norm = @(dtau) norm(dtau,'inf');
Tmax = 500;
%tol = 1e-4;
tol = .5*(1/K);

if isempty(X)
    dim = [0 0];
else
    try
        dim = [size(X{1},2) size(X{2},2)];    
    catch
        error('X should be either {}, {[],[]}, {X_1,[]}, {[],X_2}, or {X_1,X_2}.')
    end
end
X_idx = { 1:dim(1); dim(1) + (1:dim(2)) };
D = sum(dim);
COVAR = D > 0;

if ~COVAR
    diag_rest = 1;
    Sigt = ones(D,K);
    %Sigt = zeros(D,D,K);
    %for k=1:K
    %    Sigt(:,:,k) = eye(D);
    %end
    mut = zeros(K,D);
end

if ~diag_rest
    Sig = eye(D);
else
    Sig = ones(D,1);
end
mu = zeros(1,D);
sig2 = 10*[1,1];

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
    if ignore_theta, dc_str = 'no DC'; else dc_str = 'DC'; end
    if gen_sbm, gen_str = 'gen'; else gen_str = 'pp'; end
    if diag_rest, diag_str = 'diag_cov'; else diag_str = 'full_cov'; end
    if strcmpi(model,'poi')
        poi_model = 1;
        mod_str = sprintf('Poi');
    else
        poi_model = 0;
        mod_str = sprintf('Ber');
    end
    str = sprintf('mbiSBM (v.7): %4s, %4s, %3s, %3s, %3s, %3s', covar_str, dc_str, mod_str, gen_str, diag_str);
    fprintf('%30s',str)
    
end

% p = p_init;
% q = q_init;
% [pi_1, pi_2] = deal( ones(1,K)/K );

if ignore_theta
    theta{1} = ones(N(1),1); 
    theta{2} = ones(N(2),1);
end


for t = 1:Tmax   
    if manual_init && (t==1)
        [pi_1, pi_2] = deal(ones(1,K)/K);
        p = mean(A(:));
        q = 0.1*p;
        if COVAR
            Sig = eye(D);
            mu = zeros(1,D);
            sig2 = 10*[1,1];
            Sigt = zeros(D,D,K);
            for k=1:K
                Sigt(:,:,k) = eye(D);
            end
            mut = zeros(K,D);
        end
    else % initialize everything based on the initial labels
        pi_1 = update_pi(tau_1);
        pi_2 = update_pi(tau_2);
    
        if gen_sbm
            Psi = update_Psi(A, tau_1, tau_2);
        else
           [p,q] = update_Q(A, tau_1, tau_2, N(1), N(2));
        end
    
    
        if COVAR
            [Sigt, mut] = update_Gamma(X, dim, K, tau_1, tau_2, Sig, mu, sig2, diag_rest);
            [Sig, mu] = update_Sig_mu(Sigt, mut, K, diag_rest);
        end

        [beta, betap] = update_beta(X, X_idx, N, K, Sigt, mut, sig2, diag_rest);

        sig2 = update_sig2(dim, N, betap, tau_1, tau_2);
    end

    tau_1_old = tau_1;
    tau_2_old = tau_2;
    
    if ~ignore_theta
        theta{1} = dcbm_theta_update_spingarn(degs{1}, tau_1.val, varargin{:});
        theta{2} = dcbm_theta_update_spingarn(degs{2}, tau_2.val, varargin{:});
        %theta{1} = dcbm_theta_update_cham_pock(degs{1}, tau_1.val, varargin{:});
        %theta{2} = dcbm_theta_update_cham_pock(degs{2}, tau_2.val, varargin{:});
    end
    
    if gen_sbm
        if poi_model  
            Phi_1 = log(Psi+epsil);
            Phi_0 = -Psi;
        else
            temp = log(1-Psi+epsil);
            Phi_1 = log(Psi+epsil) - temp;
            Phi_0 = temp;
        end
    else
        if poi_model  
            phi_1 = log(p+epsil) - log(q+epsil);
            phi_0 = q-p;
        else
            phi_1 = log(p*(1-q)+epsil) - log(q*(1-p)+epsil);
            phi_0 = log(1-p+epsil) - log(1-q+epsil);
        end
    end
    
      
    % [beta, betap] = update_beta(X, X_idx, N, K, Sigt, mut, sig2);
    
    if gen_sbm % general SBM
        if first_tau_update == 1
            tau_1 = update_tau_gen(A ,  tau_2, theta{1}, beta{1}, pi_1, Phi_1, Phi_0, epsil, varargin{:});
            tau_2 = update_tau_gen(A',  tau_1, theta{2}, beta{2}, pi_2, Phi_1', Phi_0', epsil, varargin{:});
        else
            tau_2 = update_tau_gen(A',  tau_1, theta{2}, beta{2}, pi_2, Phi_1', Phi_0', epsil, varargin{:});
            tau_1 = update_tau_gen(A ,  tau_2, theta{1}, beta{1}, pi_1, Phi_1, Phi_0, epsil, varargin{:});
        end
    else % planted partion model
        if first_tau_update == 1
            tau_1 = update_tau(A ,  tau_2, theta{1}, beta{1}, pi_1, phi_1, phi_0, epsil, varargin{:});
            tau_2 = update_tau(A',  tau_1, theta{2}, beta{2}, pi_2, phi_1, phi_0, epsil, varargin{:});
        else
            tau_2 = update_tau(A',  tau_1, theta{2}, beta{2}, pi_2, phi_1, phi_0, epsil, varargin{:});
            tau_1 = update_tau(A ,  tau_2, theta{1}, beta{1}, pi_1, phi_1, phi_0, epsil, varargin{:});
        end
    end
    
    
    delta_1 = tau_norm(tau_1.val - tau_1_old.val);
    delta_2 = tau_norm(tau_2.val - tau_2_old.val);
    
%     if COVAR
%         [Sigt, mut, Sig, mu, sig2] = update_params(X, dim, N, K, betap, tau_1, tau_2,  Sigt, mut, Sig, mu, sig2);
%     end
%     
%     pi_1 = update_pi(tau_1);
%     pi_2 = update_pi(tau_2);
%     [p,q] = update_Q(A, tau_1, tau_2, N(1), N(2));


    if max([delta_1,delta_2]) < tol
       break
    end
end
%if VERBOSE, fprintf('Total # of itr = %3d\n',t), end
if VERBOSE > 0,  fprintf('%3d, %3.3e, %3.3e \n',t, delta_1, delta_2);  end
tau_1 = tau_1.val;
tau_2 = tau_2.val;

if ~gen_sbm
    Psi = q*ones(K) + (p-q)*eye(K);
end
end


% function out = row_normalize(X) 
%     out = bsxfun(@times, X, 1./sum(X,2));
% end
% 
% function out = safe_exp(X) 
%     out = exp( bsxfun(@plus, X , -max(X,[],2)) );
% end

function tau_out = update_tau(A, tau, theta, beta, pr, phi_1, phi_0, epsil, varargin)
    
    N = size(A,1);
    %tau = row_normalize( safe_exp(phi_1*A*t + phi_0*repmat(tsum,N,1) * diag(pr)) );
    
    a =  phi_1*A*tau.val + phi_0*repmat(tau.sum + log(pr),N,1) + beta;
    tval = dcbm_label_update_augL(a,theta, varargin{:});
    
    %tval = row_normalize( safe_exp( phi_1*A*tau.val + phi_0*repmat(tau.sum,N,1)+beta) * diag(pr) );
    
    tsum = max( sum(tval,1), epsil ); % tsum is a 1xK row vector 
    %tsum = sum(tval,1);
    
    tau_out = struct('val', tval, 'sum',tsum);
end

function tau_out = update_tau_gen(A, tau, theta, beta, pr, Phi_1, Phi_0, epsil, varargin)
    
    N = size(A,1);
    %tau = row_normalize( safe_exp(phi_1*A*t + phi_0*repmat(tsum,N,1) * diag(pr)) );
    
    a =  A*tau.val*(Phi_1') + repmat(tau.sum*(Phi_0') + log(pr),N,1) + beta;
    tval = dcbm_label_update_augL(a,theta, varargin{:});
    
    %tval = row_normalize( safe_exp( phi_1*A*tau.val + phi_0*repmat(tau.sum,N,1)+beta) * diag(pr) );
    
    tsum = max( sum(tval,1), epsil ); % tsum is a 1xK row vector 
    %tsum = sum(tval,1);
    
    tau_out = struct('val', tval, 'sum',tsum);
end


function [p,q] = update_Q(A, tau_1, tau_2, N_1, N_2)
    temp = sum(tau_1.sum .* tau_2.sum);
  
    p = trace( tau_1.val' * (A * tau_2.val) ) / temp;
    r = N_1*N_2 / temp;
    q = (r*full(mean(A(:))) - p)/(r-1);
end

function Psi = update_Psi(A, tau_1, tau_2)

    Psi =  tau_1.val' * (A * tau_2.val)  ./ (tau_1.sum(:) * tau_2.sum(:)');
    
end

function pr = update_pi(tau)
    pr = tau.sum / sum(tau.sum);
end

% function [beta, betap] = beta_func(X,St,mt,sig2) 
%    betap = trace(St) + sum( bsxfun(@minus, X, mt).^2,2);
%    beta = - betap /(2*sig2);
% end

function [beta, betap] = update_beta(X, idx, N, K, Sigt, mut, sig2, diag_rest)
    beta = {[];[]};
    betap = {[];[]};
    for r = 1:2
        beta{r} = zeros(N(r),K);
    
        if isempty(idx{r}), continue, end
    
        for k = 1:K
            %[beta{r}(:,k), betap{r}(:,k)] = beta_func(X{r}, Sigt(idx{r},idx{r},k), mut(k,idx{r}), sig2(r));     
            if ~diag_rest
                temp = trace(Sigt(idx{r},idx{r},k));
            else % diagonal restriction, Sigt(:,k) is a vector
                temp = sum(Sigt(idx{r},k));
            end
            betap{r}(:,k) = temp + sum( bsxfun(@minus, X{r},  mut(k,idx{r})).^2, 2);
            beta{r}(:,k) = - betap{r}(:,k) /(2*sig2(r));
        end
    end    
end

function [Sigt, mut] = update_Gamma(X, d, K, tau_1, tau_2, Sig, mu, sig2, diag_rest)
    %Gamma = (Sigt, mut) update
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
    
    ds = sum(d);
    if ~diag_rest
        Sigt = zeros(ds,ds,K);
        mut = zeros(K,ds);
        %disp(Sig), disp(mu)
        Sig_inv_mu = Sig\mu';
        for k = 1:K
            D_k_inv = blkdiag( tau_sum{1}(k)/sig2(1) * eye(d(1)), tau_sum{2}(k)/sig2(2) * eye(d(2)));
            Sigt(:,:,k) = inv( D_k_inv + inv(Sig) );
            mut(k,:) = Sigt(:,:,k)*( D_k_inv*much(k,:)' + Sig_inv_mu);
        end
    else % diagonal restriction, Sig is a vector
        Sigt = zeros(ds,K);
        mut = zeros(K,ds);
        %disp(Sig), disp(mu)
        Sig_inv_mu = Sig.*mu'; % both column vectors
        for k = 1:K
            D_k_inv =  [tau_sum{1}(k)/sig2(1) * ones(d(1),1); tau_sum{2}(k)/sig2(2) * ones(d(2),1)];
            Sigt(:,k) = 1./( D_k_inv + 1./Sig );
            mut(k,:) = Sigt(:,k).*( D_k_inv.*much(k,:)' + Sig_inv_mu );
        end
    end
end

function [Sig, mu] = update_Sig_mu(Sigt, mut, K, diag_rest)
    %mu update
    mu = mean(mut);

    %Sig update
    diff = bsxfun(@minus, mut, mu);
    if ~diag_rest
        Sig = mean(Sigt,3)+ diff'*diff/K;
    else % diagonal restriction, Sig is a vector, so is Sigt(:,k)
        Sig = mean(Sigt,2) + diag(diff'*diff/K);
    end
end

function sig2 = update_sig2(d, N, betap, tau_1, tau_2)
    tau{1} = tau_1.val;
    tau{2} = tau_2.val;
    
    sig2 = zeros(2,1);
    for r = 1:2
        if d(r) > 0, sig2(r) = trace(tau{r}'*betap{r})/(N(r)*d(r)); end
    end
end
% 
% function [Sigt, mut, Sig, mu, sig2] = update_params(X, d, N, K, betap, tau_1, tau_2,  Sigt, mut, Sig, mu, sig2)
%     tau{1} = tau_1.val;
%     tau{2} = tau_2.val;
%     tau_sum{1} = tau_1.sum;
%     tau_sum{2} = tau_2.sum;
% 
%     Xch = {[],[]};
%     much = {[],[]};
%     for r = 1:2
%         if d(r) > 0
%             Xch{r} = tau{r}'*X{r}; % Dimension: K x d(r)
%             much{r} = bsxfun(@times, Xch{r}, 1./(tau_sum{r}')); % Dimension: K x d(r)
%         end
%     end
%     much = [much{:}]; % horizontal concatination
% 
%     %Sigt and mut update
%     for k = 1:K
%         D_k_inv = blkdiag( tau_sum{1}(k)/sig2(1) * eye(d(1)), tau_sum{2}(k)/sig2(2) * eye(d(2)));
%         Sigt(:,:,k) = inv( D_k_inv + inv(Sig) );
%         mut(k,:) = Sigt(:,:,k)*( D_k_inv*much(k,:)' + Sig\mu');
%     end
% 
%     %mu update
%     mu = mean(mut);
% 
%     %Sig update
%     diff = bsxfun(@minus, mut, mu);
%     Sig = mean(Sigt,3)+ diff'*diff/K;
% 
%     %sig2 update
%     sig2 = zeros(2,1);
%     for r = 1:2
%         if d(r) > 0, sig2(r) = trace(tau{r}'*betap{r})/(N(r)*d(r)); end
%     end
% end