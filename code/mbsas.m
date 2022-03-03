function [m,C]=mbsas(X,threshold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION (auxiliary)

%
% INPUT ARGUMENTS:
%   X:      lxN matrix, whose columns are the data vectors.

%
% OUTPUT ARGUMENTS:
%   C:      lxm matrix
%
% (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 1;
[C]=distant_init(X,m,1);
Nc = (1);
[l,N]=size(X);

ind = find_index(X,C);
X(:,ind) = [];

for i = 1:N-1
   minimum = 10000;
   for j = 1:m
       if distan(X(:,i),C(:,j))<minimum 
          minimum = distan(X(:,i),C(:,j));
          closest_repr = j;
          break;
       end
   end
   if minimum>threshold
   	m = m + 1;
   	C = [C X(:,i)];
   	Nc = [Nc 1];
   else
   	C(:,closest_repr) = (Nc(closest_repr)*C(:,closest_repr)+X(:,i))/(Nc(closest_repr)+1);
   	Nc(closest_repr) = Nc(closest_repr) + 1;
   end  
end



