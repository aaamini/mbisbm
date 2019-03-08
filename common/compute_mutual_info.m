function normMuI = compute_mutual_info(c,e,varargin)
% normMUI Computes the normalized mutual information between two clusters
% Labels should be either vectors or n x k matrices

c = turn_into_column_ifvec(c);
e = turn_into_column_ifvec(e);

if size(c,2) > 1
    c = label_mat2vec(c);
end

if size(e,2) > 1
    e = label_mat2vec(e);
end


CM = compute_confusion_matrix(c,e);


N = sum(CM(:));
normCM = CM/N; % normalized confusion matrix

IDX = CM == 0; % index of zero elements of CM so that we can avoid them

jointEnt = - sum( (normCM(~IDX)).*log(normCM(~IDX)) );

indpt = sum(normCM,2) * sum(normCM,1);
muInfo = sum(normCM(~IDX) .* log(normCM(~IDX) ./ indpt(~IDX)) );

normMuI = muInfo / jointEnt;
