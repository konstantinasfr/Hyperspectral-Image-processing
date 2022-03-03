%%%%%% random initialization %%%%%%%% 1 24
J_total = zeros(clusters,1);
threshold_total = zeros(clusters,1);
for i = 1:clusters
[init_theta]= rand_init(dataset,i,1);
[theta,bel,J]= k_means(dataset,init_theta);
if i == 8
          pred_labels = bel';
          visualize(pred_labels,Salinas_Labels,Salinas_Image,"random initialization - k-means");
          [max_match_perm,max_acc_perm]=my_accuracy(8,true_labels,pred_labels);
          %Perform mapping operations
          [NewLabel,max_match_hung]=BestMapping(true_labels,pred_labels);
          %Y:Real Tags Idx:Cluster label NewLabel:Map rearranged labels
          %ACC
          max_acc_hung=Acc(true_labels,NewLabel);
          fprintf('Accuracy of clustering is: %f\n',max_acc_hung); %Show classification results
          C = confusionmat(true_labels,pred_labels);
          confusionchart(C);
end
J_total(i) = J;
threshold_total(i) = i;
end
figure('Name',"random initialization"), plot(clusters_total_maxmin,J_total);


%%%%%%%%%%%%%%%%%%%%%%% mbsas inilization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[w]=distant_init(dataset,2,1);
[max_distance]=distan(w(:,1),w(:,2));
% for k = 1:1
% cols = size(dataset,2);
% P = randperm(cols);
% dataset_suffled = dataset(:,P);
% dataset = dataset_suffled;
% end
dataset_suffled = dataset;

representatives = zeros(15);
window = 100;
step = max_distance*(1/window);
clusters_mbsas_total = zeros(window-1,1);
threshold_total = zeros(window-1,1);
J_total_mbsas = [];
clusters_total_mbsas = [];
for i = 1:window-3
    threshold = max_distance - i*step;
    [clusters_mbsas,C]=mbsas(dataset_suffled,threshold);
     [~,a] = size(C);
     if representatives(a) == 0
         %do k-means
         [theta,bel,J]= k_means(dataset,C);
        if i == 6
          pred_labels = bel';
          visualize(pred_labels,Salinas_Labels,Salinas_Image,"random initialization - k-means");
          [max_match_perm,max_acc_perm]=my_accuracy(8,true_labels,pred_labels);
          %Perform mapping operations
          [NewLabel,max_match_hung]=BestMapping(true_labels,pred_labels);
          %Y:Real Tags Idx:Cluster label NewLabel:Map rearranged labels
          %ACC
          max_acc_hung=Acc(true_labels,NewLabel);
          fprintf('Accuracy of clustering is: %f\n',max_acc_hung); %Show classification results
          C = confusionmat(true_labels,pred_labels);
          confusionchart(C);
        end
        % fprintf("J and clusters %d %d \n",J,i);
         J_total_mbsas = [J_total_mbsas, J];
         clusters_total_mbsas =[clusters_total_mbsas, a];
         representatives(a)=1;
     end
    clusters_mbsas_total(i) = clusters_mbsas;
    threshold_total(i) = threshold;
end
 figure('Name',"mbsas"), histogram(clusters_mbsas_total);
 figure('Name',"mbsas1"), plot(threshold_total,clusters_mbsas_total);

 figure('Name',"mbsas initialization"), plot(clusters_total_mbsas,J_total_mbsas);
