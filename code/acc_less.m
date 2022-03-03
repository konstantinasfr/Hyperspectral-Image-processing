function []=acc_less(conf_matrix)

[classes,clusters] = size(conf_matrix);

correct_matches = zeros(clusters,1);
index = zeros(clusters,1);

i=1;
j=1;
for k = i:clusters
    
    for l = j:classes
    
        if conf_matrix(l,k) ~= 0
          regr(i+1,j+1)
            
        end
    end
    [m,i]= max(conf_matrix(k,:));
    if m > correct_matches(i)
        correct_matches(i) = m;
        index(i) = k;
    end 
    
    
    
end

