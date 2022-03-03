function [W,num_clusters]=maxmin_alg(X,threshold)

%INPUT ARGUMENTS:
%   X:              Nxl matrix, whose rows are the data vectors.
%   threshold:      the distance in which new clusters will stop be
%                   created.
%           
% OUTPUT ARGUMENTS:
%   w:      lxm matrix, whose columns are the selected "most
%           distant" vectors.
%
% W: points that define clusters
% Y=X-w

X=X';
[l,N]=size(X);
[W]=distant_init(X,2,1);
ind1 = find_index(X,W(:,1));
ind2 = find_index(X,W(:,2));
Y = X;
Y(:,ind1) = [];
Y(:,ind2) = [];

clusters_construction = 1;
while clusters_construction == 1
    [~,n]=size(Y);
    total_mins = zeros(1,n);
    for i = 1:n
        [~,clusters]=size(W);
        minimum = 100000;
        for j = 1:clusters
            if distan(Y(:,i),W(:,j))<minimum
                minimum = distan(Y(:,i),W(:,j));
%                 break;
            end
            total_mins(i) = minimum;
        end
    end
    [maximum,ind] = max(total_mins);
%     disp(maximum);
    if maximum>threshold
        W = [W Y(:,ind)];
        Y(:,ind) = [];
    else
        clusters_construction = 0;
    end
    [~,num_clusters] = size(W);
end
    
    
    
