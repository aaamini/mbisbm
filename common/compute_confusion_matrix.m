function M = compute_confusion_matrix(c,e,varargin)
% Compute the confusion matrix between labels "c" and "e"
%
% c,e Two sets of labels
% K   number of labels in both "c" and "e"

if nargin > 2
    K = varargin{1};
else
    K = max(max(c),max(e));
end

M = label_vec2mat(c,K)'*label_vec2mat(e,K);
    
% M = zeros(K);
% for k = 1:K
%     for r = 1:K
%         M(k,r) = sum( (c(:) == k) .* ( e(:) == r ) );
%     end
% end