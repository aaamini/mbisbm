function out = row_l2_normalize(X) 
out = zeros(size(X));
l2norms = row_l2_norms(X);
nzidx = l2norms ~= 0;
out(nzidx,:) = bsxfun(@times, X(nzidx,:), 1./l2norms(nzidx,:));
%out = bsxfun(@times, X, 1./row_l2_norms(X));