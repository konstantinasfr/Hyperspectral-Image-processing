% pkg load statistics
load data_country
Z = zscore(Countrydata);
subdata=Z(:,[1,2,3,5,6]);
data=subdata.';

J_total = [];
threshold_total = [];
for i = 2:10
[bel,cost,w,a]=k_medoids(data,i,30);
J_total(i) = cost;
threshold_total(i) = i;
end
plot(threshold_total,J_total) , axis([2 10 0 1200])

