function out = row_l2_norms(X)
out = sqrt(sum(X.^2,2));
