function [label_1 ,label_2, Z_2]= biSpecClust(A,K,varargin)
parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'normalize',1)
addOptional(parser,'perturb',0)
addOptional(parser,'alpha',0.1) % default perturbation parameter
addOptional(parser,'pert_geom',1) % geometric versus arithmatic scaling
addOptional(parser,'idx',1:K) % indices of eigenvector to keep for low-d representation.


parse(parser, varargin{:});
normalize = parser.Results.normalize;
perturb = parser.Results.perturb;
alpha = parser.Results.alpha;
pert_geom = parser.Results.pert_geom;
idx = parser.Results.idx;

N = size(A);

D_1 = sum(A,2);
D_2 = sum(A,1);

if pert_geom % scaling of the perturbation
    alphas = alpha/sqrt(prod(N)); %geometric mean 
else
    alphas = alpha/(sum(N)/2); %arithmatic mean 
end

if ~perturb
    g1 = safe_diag_pwr(D_1);
    g2 = safe_diag_pwr(D_2);
    A_n = diag(g1)*A*diag(g2);
    [U,~,V] = svds(A_n,K+1);
else
    g1 = safe_diag_pwr(D_1(:) + alphas*N(2)*ones(N(1),1));
    g2 = safe_diag_pwr(D_2(:) + alphas*N(1)*ones(N(2),1));
    
    %[U,~,V] = svds(@(x,tflag) Afun(x,tflag,A,g1,g2,alphas), N, K+1);
    Af = struct;
    Af.size = N;
    Af.times = @(x,varargin) bsxfun(@times,g1,A*(bsxfun(@times,g2,x)))+ alphas*g1*(g2'*x); %Afun(x,'notransp',A,g1,g2,alphas);
    Af.trans = @(x,varargin) bsxfun(@times,g2,A'*bsxfun(@times,g1,x)) + alphas*g2*(g1'*x); %Afun(x,'transp',A,g1,g2,alphas);
    Af.param = [];
    opts.tol = 1e-8;
    [U,~,V,] = lmsvd(Af,K+1,opts);
    %[U,~,V] = lansvd(@(x) Atfun(x,A*1.,g1,g2,alphas), @(x) Atfun(x,A'*1.,g2,g1,alphas), N(1), N(2), K+1);
end

U_2 = U(:,idx);
V_2 = V(:,idx);
    
    %row_l2_norms = @(X) sqrt(sum(X.^2,2));
    %row_l2_normalize = @(X) bsxfun(@times, X, 1./row_l2_norms(X));
if normalize
    U_2=row_l2_normalize(U_2);
    V_2=row_l2_normalize(V_2);
end

Z_2 = [U_2;V_2];
    
%[e, C] = kmeans(Z_2,K,'Replicates',50, 'Start','sample');
[e, C] = kmeans(Z_2,K,'Replicates',100);
label_1 = e(1:N(1));
label_2 = e(N(1) + (1:N(2)));

%figure(1), clf, scatter3(Z_2(:,1),Z_2(:,2),Z_2(:,3),'.'), hold on, scatter3(C(:,1), C(:,2), C(:,3),'ro')
        
end

function g = safe_diag_pwr(d)
    idx = d ~= 0;
    g = zeros(size(d));
    g(idx) = d(idx).^(-.5);
    %G = diag(dinv);
end

function y = Afun(x,tflag,A,g1,g2,alphas)
    if strcmp(tflag,'notransp')
        y = g1.*(A*(g2.*x)) + alphas*g1*(g2'*x);
    else
        y = g2.*(A'*(g1.*x)) + alphas*g2*(g1'*x);
    end
end


% function y = Ax(x,A,g1,g2,alphas)
%      y = g1.*(A*(g2.*x)) + alphas*g1*(g2'*x);
% end
