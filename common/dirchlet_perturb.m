function [z1_p,z2_p] = dirchlet_perturb(z1,z2,t,alpha)

if nargin < 4
    alpha = .5;
end

[N1,K] = size(z1);
[N2,Kp] = size(z2);
if Kp ~= K
 error('2nd dimension of z1 and z2 should have the same size.')
end

 
da = alpha*ones(1,K); 
z1_p = (1-t)*z1 + t*drchrnd(da,N1);
z2_p = (1-t)*z2 + t*drchrnd(da,N2);