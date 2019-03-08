function [theta, delta] = dcbm_theta_update_cham_pock(d,tau,varargin)
%function [theta, delta] = dcbm_theta_update_cham_pock(d,tau,alpha,r,sig,VERBOSE)

parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'alpha',1)
addOptional(parser,'r',1)
addOptional(parser,'sig',1)
addOptional(parser,'VERBOSE',0)

parse(parser, varargin{:});
alpha = parser.Results.alpha;
r = parser.Results.r;
sig = parser.Results.sig;
VERBOSE = parser.Results.VERBOSE;

[n,K] = size(tau);

% if nargin < 3
%     alpha = 1;
% end
% 
% if nargin < 4
%     r = 1;
% end
% 
% if nargin < 5
%     sig = 1;
% end
% 
% if nargin < 6 
%     VERBOSE = 0;
% end
% 

tol = 1e-4;
T = 100;
    function out = err(x,xo)
        out = norm(x-xo,'inf')/max(norm(xo,'inf'),1);
    end
    
    

    function out = prox_theta(th)
        out = 0.5*th + sqrt( 0.25*th.^2 + alpha*d);
    end
 

    xi = zeros(K,1);
    theta = zeros(n,1);
    thb = zeros(n,1);
    delta = zeros(T,1);
    for t = 1:T
        
        theta_old = theta;
        
        xi = xi + sig*tau'*(thb - 1);
        theta = prox_theta(theta_old - alpha*tau*xi);
        thb = theta + r*(theta - theta_old);
        
        
        
        delta(t) = err(theta,theta_old);
        if delta(t) < tol
            delta(t+1:end) = [];
            break
        end
        if VERBOSE > 1, fprintf('cham_pock: (alpha,r,sig) = %3.2f, %3.2f, %3.2f ----> %3d, delta = %3.3e\n',alpha,r,sig,t,delta(t)), end
    end
    %fprintf('norm (th-1)''tau = %5.3f\n', norm(theta-ones(n,1))'*tau);

end

