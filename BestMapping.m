function [NewLabel,max_match_hung] = BestMapping(La1,La2)

%Real label: La1 clustering result label: La2 mapped label: NewLabel

Label1=unique(La1');
L1=length(Label1);
Label2=unique(La2');
L2=length(Label2);

%Construct a matrix G that calculates the repeatability of two classification labels
G = zeros(max(L1,L2),max(L1,L2));
for i=1:L1
    index1= La1==Label1(1,i);
    for j=1:L2
        index2= La2==Label2(1,j);
        G(i,j)=sum(index1.*index2);
    end
end

%Use Hungarian algorithm to calculate the matrix after rearrangement
[index]=munkres(-G);
[l,l] = size(index);
max_match_hung = zeros(l,1);
for i = 1:l
   for j = 1:l
      if index(j,i) ~= 0
           max_match_hung(i) = j;
      end
   end
end
%Convert the map reordering result to a row vector that stores the label order after map reordering
[temp]=MarkReplace(index);
%Generate label NewLabel after rearrangement
NewLabel=zeros(size(La2));
for i=1:L2
    NewLabel(La2==Label2(i))=temp(i);
end

end
