function [A, label, theta, emp_avg_deg, Q, pr, expected_A] = bip_dcsbm_pareto(N, K, beta, lambda, varargin)

parser = inputParser;
addOptional(parser,'alpha',inf)
addOptional(parser,'pr_max_min_ratio',1.7)
addOptional(parser,'poisson',false)

parse(parser, varargin{:});
alpha = parser.Results.alpha;
pr_max_min_ratio = parser.Results.pr_max_min_ratio;
poisson = parser.Results.poisson;

if isinf(alpha)
    DEG_VAR = false;
else
    DEG_VAR = true;
    simple_pareto = @(theta,alpha) makedist('GeneralizedPareto','k',1/alpha,'sigma',theta/alpha,'theta',theta) ;
    pd = simple_pareto((alpha-1)/alpha,alpha);
end

pr = {[];[]};
for r = 1:2
    pr{r} = (pr_max_min_ratio-1)*rand(K,1)+1;
    pr{r} = pr{r}/sum(pr{r});
end


compute_avg_deg = @(N,p,q,pri1,pri2) (2*prod(N)/sum(N)) * (q + (p-q)*sum(pri1(:).*pri2(:)));

deg_scale_factor = compute_avg_deg(N, 1, beta, pr{1}, pr{2});
p = (lambda/deg_scale_factor)*1;
q = p*beta;

Q = [p,q];

label = {[];[]};
theta = {[];[]};

for r = 1:2
    label{r} = gen_labels(N(r),pr{r}, K);
    
    if DEG_VAR, theta{r} = random(pd, N(r),1); end
end


prob_matrix = q*ones(N(1),N(2));


for k=1:K
    idx = {[];[]};
    for r = 1:2
        idx{r} = label{r}==k;
        if DEG_VAR
            factor = sum(theta{r}(idx{r})) / sum(idx{r});
            theta{r}(idx{r}) = theta{r}(idx{r})/factor;
        end
    end
    
    prob_matrix(idx{1},idx{2}) = p; 
    
end

if DEG_VAR
    expected_A = (theta{1}*theta{2}').*prob_matrix;
else
    expected_A = prob_matrix;
end

above_one = sum(expected_A(:) > 1) / (prod(N));
if  above_one > 0, warning('%3.3e fraction of edge probabilities are larger than 1.', above_one), end

if poisson
    A = poissrnd(expected_A,N(1),N(2));
else
    A = rand(N(1),N(2)) <= expected_A;
end

emp_avg_deg = 2*sum(A(:))/sum(N);
end


function [labels, sizes] = gen_labels(n,pr,K)
%K should devise n 
%sizes = repmat(n/K,1,K);

sizes = mnrnd(n,pr);
labels = [];
for i=1:K
  labels = [labels repmat(i,1,sizes(i))];
end

end


