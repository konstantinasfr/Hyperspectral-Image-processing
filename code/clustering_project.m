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
[eigenval,eigenvec,explain,Y,mean_vec]=pca_fun(N',9);


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
%   figure('Name',"variance error"), plot(cummulative_explain(1:15));

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
      if i == 10
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
%   figure('Name',"maxmin algorithm initialization"), plot(clusters_total_maxmin,J_total);
 
 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fuzzy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
q = 1.5;
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
    if i == 10
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
[eigenval,eigenvec,explain,Y,mean_vec]=pca_fun(N',9);
Y=Y';

existed_1=(true_labels==1);   %This contains 1 in the positions corresponding to pixels with known class label
Y1=Y(existed_1,:);
[D1]=find_eta(Y1);
eta1 = mean(D1);

existed_2=(true_labels==2);   %This contains 1 in the positions corresponding to pixels with known class label
Y2=Y(existed_2,:);
[D2]=find_eta(Y2);
eta2 = mean(D2);

existed_3=(true_labels==3);   %This contains 1 in the positions corresponding to pixels with known class label
Y3=Y(existed_3,:);
[D3]=find_eta(Y3);
eta3 = mean(D3);

existed_4=(true_labels==4);   %This contains 1 in the positions corresponding to pixels with known class label
Y4=Y(existed_4,:);
[D4]=find_eta(Y4);
eta4 = mean(D4);

existed_5=(true_labels==5);   %This contains 1 in the positions corresponding to pixels with known class label
Y5=Y(existed_5,:);
[D5]=find_eta(Y5);
eta5 = mean(D5);

existed_6=(true_labels==6);   %This contains 1 in the positions corresponding to pixels with known class label
Y6=Y(existed_6,:);
[D6]=find_eta(Y6);
eta6 = mean(D6);

existed_7=(true_labels==7);   %This contains 1 in the positions corresponding to pixels with known class label
Y7=Y(existed_7,:);
[D7]=find_eta(Y7);
eta7 = mean(D7);

existed_8=(true_labels==8);   %This contains 1 in the positions corresponding to pixels with known class label
Y8=Y(existed_8,:);
[D8]=find_eta(Y8);
eta8 = mean(D8);


%{
figure(2), hold on
figure(2), plot(Y1(:,1),Y1(:,2),'m.',...
Y2(:,1),Y2(:,2),'c.',Y3(:,1),Y3(:,2),'g.',...
Y4(:,1),Y4(:,2),'y.',Y5(:,1),Y5(:,2),'k.',Y6(:,1),Y6(:,2),'r.',Y7(:,1),Y7(:,2),'b.',Y8(:,1),Y8(:,2),'.')
% figure(2), plot(Y1(:,1),Y1(:,2),'mp',...
% Y2(:,1),Y2(:,2),'c*',Y3(:,1),Y3(:,2),'go',...
% Y4(:,1),Y4(:,2),'yx',Y5(:,1),Y5(:,2),'k+',Y6(:,1),Y6(:,2),'rd',Y7(:,1),Y7(:,2),'bs',Y8(:,1),Y8(:,2),'.')
% figure(2), plot(theta(1,:),theta(2,:),'k+',init_theta(1,:),init_theta(2,:),'rd',m(:,1),m(:,2),'bs')
figure(2), axis equal
%}

%%%%%%%%%%%%%% Probabilistic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% GMModel = fitgmdist(dataset',8);
[pred_labels,c1,p1,JDF1,iter1]=PDclust(dataset',9);
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
Z = linkage(Y,"ward");
% T = cluster(Z,'Cutoff',C);
%     T = cluster(Z,'Cutoff',145,'Criterion','distance');

    T = cluster(Z,'maxclust',10);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Possibilistic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
clusters = 10;
[l,N]=size(dataset);
objFun_total = zeros(clusters,1);
i_total = zeros(clusters,1);
max_bel = zeros(1,167);
% eta = [ 10 1 1 1 1 1 1 1];
%  eta = [ eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8];
%  eta = [2.44,2.38,3.91,7.67,2.82,2.23,1.97,3.24];
% eta = [ 3 1 2 5 6 7 8 4];
eta = [2.44,1,3,1.67,2.82,2.23,1.97,1.24];
eta = [2.7,1,1.6,2.1,0.4,3.641,1.3,1.35,0.9];
eta = [2.7,1,1.6,1.8,0.35,3.641,1.7,1.35,0.9,1.2]; %0.753020
eta = [2.7,1,1.6,1.8,0.33,3.641,1.6,1.35,0.9,1.2];
eta = [2.7,1,1.6,1.7,0.33,3.641,1.6,1.35,0.9,1.15];
eta = [2.7,1,1.6,1.7,0.31,3.641,1.5,1.35,0.9,1.15];
eta = [2.7,1,1.6,1.2,0.31,3.641,1.5,1.35,0.9,1.15];
eta = [2.7,1,1.6,1.2,0.31,3.4,1.3,1.35,0.7,1.1];
q=1.5;
lll = clusters;
for i = lll:lll
[U,theta]=possibi(dataset,i,eta,q,1,3,0.00001);
    for j = 1:13908
         [maximum_U,maximum_cluster] = max(U(j,:));
          max_bel(j) = maximum_cluster;
    end
          if i == 10
              pred_labels = max_bel';
              visualize(pred_labels,Salinas_Labels,Salinas_Image,"possibilistic");
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
        
%     objFun_total(i) = objFun(end);
%     i_total(i) = i;
end


%}


