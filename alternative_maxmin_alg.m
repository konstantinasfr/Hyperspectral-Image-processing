function [W,threshold]=alternative_maxmin_alg(X,num_clusters,sed)

%INPUT ARGUMENTS:
%   X:              lxN matrix,  each column of which corresponds to
%           an l-dimensional data vector.
%   threshold:      the distance in which new clusters will stop be
%                   created.
%           
% OUTPUT ARGUMENTS:
%   w:      lxm matrix, whose columns are the selected "most
%           distant" vectors.
%
% W: points that define clusters
% Y=X-w



[W]=distant_init(X,num_clusters,sed);

[l,N]=size(W);

minimum=10000;
for i = 1:N-1
    for j = i+1:N
        if distan(W(:,i),W(:,j))<minimum
            minimum = distan(W(:,i),W(:,j));
        end
    end
end
threshold = minimum;
