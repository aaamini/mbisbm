function [theta, delta] = dcbm_theta_update_spingarn(d,tau,varargin)
%function [theta, delta] = dcbm_theta_update_spingarn(d,tau,alpha,VERBOSE)

parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'alpha',1)
addOptional(parser,'VERBOSE',0)

parse(parser, varargin{:});
alpha = parser.Results.alpha;
VERBOSE = parser.Results.VERBOSE;

[n,K] = size(tau);

%H = tau* ((tau'*tau)\tau');
%H = tau* (pinv(tau'*tau)*tau');
%TAU = pinv(tau'*tau);
TAU = (tau'*tau) + (1e-7)*eye(K);

tol = 1e-4;
T = 1000;
    function out = err(x,xo)
        out = norm(x-xo,'inf')/max(norm(xo,'inf'),1);
    end
    
    function out =  proj_V(theta)
        %out = theta - H*(theta - 1);
        %out = theta - tau*lsqr(tau,theta - 1);
        %out = theta - tau*(TAU*(tau'*(theta - 1)));
        out = theta - tau*(TAU\(tau'*(theta - 1)));
    end

    function out = prox_theta(th)
        out = 0.5*th + sqrt( 0.25*th.^2 + alpha*d);
    end
 

    xi = zeros(n,1);
    theta = zeros(n,1);
    delta = zeros(T,1);
    for t = 1:T
        
        theta_old = theta;
        
        theta = prox_theta(xi);
        xi = xi + proj_V(2*theta - xi) - theta;
        
        
        delta(t) = err(theta,theta_old);
        if delta(t) < tol
            delta(t+1:end) = [];
            break
        end
        if VERBOSE > 1, fprintf('Spingarn: alpha = %3.2f ----> %3d, delta = %3.3e\n',alpha,t,delta(t)), end
    end
    %fprintf('norm (th-1)''tau = %5.3f\n', norm(theta-ones(n,1))'*tau);

end

