% This is for case study for the 10th Lecture in the "Fundamentals of
% Remote Sensing" course (see the scenario in the slides)

clear
format compact
close all

load Salinas_cube
load Salinas_ground_truth
[p,n,l]=size(Xinit);

X_all=reshape(Xinit,p*n,l)'; % Taking the lx(p*n) matrix with the vectors of the pixels in its rows

L_all=reshape(Labels,p*n,1); % Taking a p*n column array of the pixel labels


for i=1:l
    figure(1), imagesc(Xinit(:,:,i)); axis off; axis image
    pause(0.2)
end

%Choose the data points that have available ground truth
Q_ground=find(L_all>0);
X=X_all(:,Q_ground);
L=L_all(Q_ground);

%Perform PCA 
[eigenval,eigenvec,explain,Y,mean_vec]=pca_fun(X,3);

%Plot the 3-dimensional points
figure; plot3(Y(1,:),Y(2,:),Y(3,:),'.b')

% Creating the first three principal components
Y=Y-min(min(Y)); % Making Y>=0
PCmat=zeros(3,p*n);
PCmat(1,Q_ground)=Y(1,:);
PCmat(2,Q_ground)=Y(2,:);
PCmat(3,Q_ground)=Y(3,:);
PCmat=PCmat';
PCmat=reshape(PCmat,p,n,3);
PCmat=(PCmat-min(min(min(PCmat))) )/( max(max(max(PCmat)))-min(min(min(PCmat))) );

% Depicting the first PC's of the salinas image cube
for i=1:3
    figure(i+11), imagesc(PCmat(:,:,i)); colorbar; axis off; axis image
    pause(1)
end

%This depicts a 3-channel image where each channel corresponds to one of
%the first 3 PCs
figure(10); imshow(PCmat)

% Saving the data 
save data X Y L Q_ground p n