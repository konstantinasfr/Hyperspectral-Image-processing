function []=visualize(cl_label,Salinas_Labels,Salinas_Image,label_name)

[p,n,l]=size(Salinas_Image); % Size of the Salinas cube
X_total=reshape(Salinas_Image, p*n,l);
L=reshape(Salinas_Labels,p*n,1);
existed_L=(L>0);  

%This contains 1 in the positions corresponding to pixels with known class label
% The following code can be used after the execution of an algorithm 
% Let "cl_label" be the px-dimensional vector, whose i-th element is the label
% of the class where the i-th vector in X has been assigned
% The code below, helps in depicting the results as an image again

cl_label_tot=zeros(p*n,1);
cl_label_tot(existed_L)=cl_label;
disp(size(cl_label_tot));
im_cl_label=reshape(cl_label_tot,p,n);
disp(size(im_cl_label));
figure('Name',label_name), imagesc(im_cl_label), axis equal
figure('Name',"Correct clustering"), imagesc(Salinas_Labels), axis equal
