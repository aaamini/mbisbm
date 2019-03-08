function [ ac ] = compute_acc( c, e )
%ACC Summary of this function goes here
%   Detailed explanation goes here


c = turn_into_column_ifvec(c);
e = turn_into_column_ifvec(e);


if size(c,2) > 1
    c = label_mat2vec(c);
end

if size(e,2) > 1
    e = label_mat2vec(e);
end

n = length(c);
if length(e) ~= n, error('c and e should have the same length.'), end

CM = compute_confusion_matrix(c,e);

[assignment,cost] = munkres(-CM/n);
ac = -cost;

% fprintf('new acc')

% K = size(y,2);
% L = size(z,2);
% ac=0;
% if K==L
%     n = size(y,1);
%     R = y'*z./n;
%     for k = 1:K
%         [~,ind]=max(R);
%         [v2,ind1]=max(max(R));
%         ac = ac+v2;
%         if k<K
%             i = ind(ind1);
%             j = ind1;
%             R(i,:)=[];
%             R(:,j)=[];
%         end
%     end
% end
% end

