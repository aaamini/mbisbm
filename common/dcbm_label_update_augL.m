function [X,lambda] = dcbm_label_update_augL(a, theta, varargin)
%function [X,lambda] = dcbm_label_update_augL(a,theta,stepsize)

parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'mu',0.01)
addOptional(parser,'Tmax',100)
addOptional(parser,'tol',1e-4)
addOptional(parser,'damping',0.95)

parse(parser, varargin{:});
mu = parser.Results.mu;
Tmax = parser.Results.Tmax;
tol = parser.Results.tol;
damping = parser.Results.damping;

% if nargin < 3
%     mu = .01;
% else
%     mu = stepsize;
% end

[n,K] = size(a);

theta = theta(:);
if length(theta) ~= n
    error('Length of theta should match size(a,1).')
end

%tol = 1e-4;
%T = 100;

lambda = zeros(1,K);

tht = theta - ones(n,1);

for t = 1:Tmax
    X = row_softmax(a + tht*lambda);
    lambda_old = lambda;
    lambda = lambda - mu*(tht'*X);
    
    mu = mu*damping;
    
    delta = norm(lambda-lambda_old)/max(1,norm(lambda_old));
    if delta < tol
        break
    end
   
end
%fprintf('tau-update in %3d iters\n',t)