function [max_perm,max_acc]=my_accuracy(n,true_labels,pred_labels)

v = (1:n);
P = perms(v);
disp(P(1,:));


[perm,kk]=size(P);
[pixels,zz] = size(true_labels);

max_acc = 0;
for i = 1:perm
    pred_labels_perm = interp1(v,P(i,:),pred_labels);
    same = 0;
    for j = 1:pixels
%         disp(true_labels);
        if true_labels(j)== pred_labels_perm(j)
           same = same+1;
        end
    end
    temp_acc = same/pixels;
    if temp_acc>max_acc
        max_acc = temp_acc;
        max_perm = P(i,:);
    end
    
end
    