clear
format compact
close all

load Salinas_Data

[p,n,l]=size(Salinas_Image); % Size of the Salinas cube

disp([p,n,l]);
% disp(Salinas_Image(6,1,:));

% %Depicting the bands of the Salinas cube
%  for i=1:l
%      figure(1); imagesc(Salinas_Image(:,:,i))
% %      pause(0.1)
%  end

 % Making a two dimensional array whose rows correspond to the pixels and
% the columns to the bands, containing only the pixels with nonzero label.
X_total=reshape(Salinas_Image, p*n,l);
L=reshape(Salinas_Labels,p*n,1);
existed_L=(L>0);   %This contains 1 in the positions corresponding to pixels with known class label
X=X_total(existed_L,:);
true_labels = L(existed_L,:);
[px,nx]=size(X); % px= no. of rows (pixels) and nx=no. of columns (bands)


N = normalize(X);
[eigenval,eigenvec,explain,Y,mean_vec]=pca_fun(N',7);


%{
 s=Y';
 x=s(:,1);
 y=s(:,2);
% %  y2=s(:,3);
%  plot(x,y,".");
%  
% Plot the clusters
figure(2), hold on
figure(2), plot(x,y,'g.')
% figure(2), plot(mean_vec,'rd')
figure(2), axis equal


%}

cummulative_explain = zeros(15);

cummulative_explain(1) = explain(1);
for i = 2:16
    cummulative_explain(i) = explain(i) + cummulative_explain(i-1);
end
%  figure('Name',"variance error"), plot(cummulative_explain(1:15));

%  figure('Name',"eigenval"), plot(eigenval);
 
 
 
  Y = Y';
 minimum = min(Y);
maximum = max(Y);
range = maximum - minimum;
disp(maximum);
disp(minimum);
disp(range);
%  for i = 1:4
%       cl_label_tot=zeros(p*n,1);
%       cl_label_tot(existed_L)=Y(:,i);
%       im_cl_label=reshape(cl_label_tot,p,n);
%       figure(i), imagesc(im_cl_label), axis equal
%  end
%   figure(11), imagesc(Salinas_Labels), axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% k-means %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dataset = Y';
clusters = 15;
%%%%%% maxmin algorithm initialization %%%%%%%%

%{
J_total = zeros(clusters,1);
clusters_total_maxmin = zeros(clusters,1);
for i = 1:clusters
    [W,threshold1]=alternative_maxmin_alg(dataset,i,173);
%     [init_theta]= rand_init(dataset,i,95);
    [theta,bel,J]= k_means(dataset,W);
      if i == 9
          pred_labels = bel';
          visualize(pred_labels,Salinas_Labels,Salinas_Image,"maxmin algorithm initialization - k-means");
          [max_match_perm,max_acc_perm]=my_accuracy(8,true_labels,pred_labels);
          %Perform mapping operations
          [NewLabel,max_match_hung]=BestMapping(true_labels,pred_labels);
          %Y:Real Tags Idx:Cluster label NewLabel:Map rearranged labels
          %ACC
          max_acc_hung=Acc(true_labels,NewLabel);
%           fprintf('Accuracy of clustering is: %f\n',max_acc_hung); %Show classification results
          C = confusionmat(true_labels,pred_labels);
          confusionchart(C);
          Suuuum = sum(C);
          disp(sum(Suuuum));
      end
    % fprintf("J and clusters %d %d \n",J,i);
    J_total(i) = J;
    clusters_total_maxmin(i) = i;
end
%    figure('Name',"maxmin algorithm initialization"), plot(clusters_total_maxmin,J_total);
 
 %}
 


%%%%%% random initialization %%%%%%%% 1 24

%{
J_total = zeros(clusters,1);
threshold_total = zeros(clusters,1);
for i = 1:clusters
[init_theta]= rand_init(dataset,i,24);
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
figure('Name',"random initialization"), plot(threshold_total,J_total);

%}

%{
q = 1.4;
clusters = 15;
options = [q NaN NaN 0];

[~,total_pixels] = size(dataset);
objFun_total = zeros(clusters,1);
i_total = zeros(clusters,1);
max_bel = zeros(1,total_pixels);
for i = 2:clusters
[centers,U,objFun] = fcm(dataset',i, options);
    for j = 1:total_pixels
         [maximum_U,maximum_cluster] = max(U(:,j));
          max_bel(j) = maximum_cluster;
    end
    if i == 8
              pred_labels = max_bel';
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
        
    objFun_total(i) = objFun(end);
    i_total(i) = i;
end
figure('Name',"fuzzy"), plot(i_total,objFun_total);
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Possibilistic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
X_total=reshape(Salinas_Image, p*n,l);
L=reshape(Salinas_Labels,p*n,1);
existed_1=(L==1);   %This contains 1 in the positions corresponding to pixels with known class label
X1=X_total(existed_1,:);
% true_labels = L(existed_L,:);

N = normalize(X1);
[eigenval,eigenvec,explain,Y1,mean_vec]=pca_fun(N',2);
Y1=Y1';
[A1,C1]=find_eta(X1);

existed_2=(L==2);
X2=X_total(existed_2,:);
N = normalize(X2);
[eigenval,eigenvec,explain,Y2,mean_vec]=pca_fun(N',2);
Y2=Y2';
[A2,C2]=find_eta(X2);

existed_3=(L==3);
X3=X_total(existed_3,:);
N = normalize(X3);
[eigenval,eigenvec,explain,Y3,mean_vec]=pca_fun(N',2);
Y3=Y3';
[A3,C3]=find_eta(X3);

existed_4=(L==4);
X4=X_total(existed_4,:);
N = normalize(X4);
[eigenval,eigenvec,explain,Y4,mean_vec]=pca_fun(N',2);
Y4=Y4';
[A4,C4]=find_eta(X4);

existed_5=(L==5);
X5=X_total(existed_5,:);
N = normalize(X5);
[eigenval,eigenvec,explain,Y5,mean_vec]=pca_fun(N',2);
Y5=Y5';
[A5,C5]=find_eta(X5);

existed_6=(L==6);
X6=X_total(existed_6,:);
N = normalize(X6);
[eigenval,eigenvec,explain,Y6,mean_vec]=pca_fun(N',2);
Y6=Y6';
[A6,C6]=find_eta(X6);

existed_7=(L==7);
X7=X_total(existed_7,:);
N = normalize(X7);
[eigenval,eigenvec,explain,Y7,mean_vec]=pca_fun(N',2);
Y7=Y7';
[A7,C7]=find_eta(X7);

existed_8=(L==8);
X8=X_total(existed_8,:);
N = normalize(X8);
[eigenval,eigenvec,explain,Y8,mean_vec]=pca_fun(N',2);
Y8=Y8';
[A8,C8]=find_eta(X8);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_total=reshape(Salinas_Image, p*n,l);
L=reshape(Salinas_Labels,p*n,1);
existed_L=(L>0);   %This contains 1 in the positions corresponding to pixels with known class label
X=X_total(existed_L,:);
true_labels = L(existed_L,:);
[px,nx]=size(X); % px= no. of rows (pixels) and nx=no. of columns (bands)


N = normalize(X);
[eigenval,eigenvec,explain,Y,mean_vec]=pca_fun(N',2);
Y=Y';

existed_1=(true_labels==1);   %This contains 1 in the positions corresponding to pixels with known class label
Y1=Y(existed_1,:);

existed_2=(true_labels==2);   %This contains 1 in the positions corresponding to pixels with known class label
Y2=Y(existed_2,:);

existed_3=(true_labels==3);   %This contains 1 in the positions corresponding to pixels with known class label
Y3=Y(existed_3,:);

existed_4=(true_labels==4);   %This contains 1 in the positions corresponding to pixels with known class label
Y4=Y(existed_4,:);

existed_5=(true_labels==5);   %This contains 1 in the positions corresponding to pixels with known class label
Y5=Y(existed_5,:);

existed_6=(true_labels==6);   %This contains 1 in the positions corresponding to pixels with known class label
Y6=Y(existed_6,:);

existed_7=(true_labels==7);   %This contains 1 in the positions corresponding to pixels with known class label
Y7=Y(existed_7,:);

existed_8=(true_labels==8);   %This contains 1 in the positions corresponding to pixels with known class label
Y8=Y(existed_8,:);

%}


figure(2), hold on
figure(2), plot(Y1(:,1),Y1(:,2),'g.',...
Y2(:,1),Y2(:,2),'y.',Y3(:,1),Y3(:,2),'m.',...
Y4(:,1),Y4(:,2),'c.',Y5(:,1),Y5(:,2),'r.',Y6(:,1),Y6(:,2),'k.',Y7(:,1),Y7(:,2),'.',Y8(:,1),Y8(:,2),'b.')
% figure(2), plot(Y1(:,1),Y1(:,2),'mp',...
% Y2(:,1),Y2(:,2),'c*',Y3(:,1),Y3(:,2),'go',...
% Y4(:,1),Y4(:,2),'yx',Y5(:,1),Y5(:,2),'k+',Y6(:,1),Y6(:,2),'rd',Y7(:,1),Y7(:,2),'bs',Y8(:,1),Y8(:,2),'.')
% figure(2), plot(theta(1,:),theta(2,:),'k+',init_theta(1,:),init_theta(2,:),'rd',m(:,1),m(:,2),'bs')
figure(2), axis equal


%{

% GMModel = fitgmdist(dataset',8);
[pred_labels,c1,p1,JDF1,iter1]=PDclust(dataset',7);
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
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Agglomerative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clusters=8;
C=100000000000000;
D = 10000;
% Y = pdist(dataset');
% squareform(Y);
Y=dataset';
Z = linkage(Y,"median");
% T = cluster(Z,'Cutoff',C);
%     T = cluster(Z,'Cutoff',140,'Criterion','distance');

    T = cluster(Z,'maxclust',8);
%  T = cluster(Z,'Cutoff',C,'Depth',D);
% cutoff = median([Z(end-7,clusters) Z(end-6,clusters) Z(end-5,clusters) Z(end-4,clusters) Z(end-3,clusters) Z(end-2,clusters) Z(end-1,clusters)]);
% cutoff = median([Z(end-2,clusters) Z(end-1,clusters)]);
 dendrogram(Z)
pred_labels = T;
visualize(pred_labels,Salinas_Labels,Salinas_Image,"agglomerative");
%               [max_match_perm,max_acc_perm]=my_accuracy(8,true_labels,pred_labels);
              %Perform mapping operations
              [NewLabel,max_match_hung]=BestMapping(true_labels,pred_labels);
              %Y:Real Tags Idx:Cluster label NewLabel:Map rearranged labels
              %ACC
              max_acc_hung=Acc(true_labels,NewLabel);
              fprintf('Accuracy of clustering is: %f\n',max_acc_hung); %Show classification results
              C = confusionmat(true_labels,pred_labels);
              confusionchart(C);
Cr= crosstab(T,true_labels);





