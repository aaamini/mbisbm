function out = row_normalize_ell1(X)

out = bsxfun(@rdivide, X, sum(abs(X),2));