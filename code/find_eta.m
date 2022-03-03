function [dist]=find_eta(X)

[N,l]=size(X); % px= no. of rows (pixels) and nx=no. of columns (bands)


M = mean(X);
dist = zeros(N,1);
for i = 1:N
      
     dist(i)=   sqrt(sum((X(i,:) - M) .^ 2));
end



% C = cov(dataset);