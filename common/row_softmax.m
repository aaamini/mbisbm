function out = row_softmax(X)

%out = row_normalize_ell1( exp(X) );
out = row_normalize_ell1( exp(bsxfun(@minus, X, max(X,[],2))) );
