function [index] = find_index(X,w)
[l,N]=size(X);
for i = 1:N
    if X(:,i) == w
        index = i;
    end
end