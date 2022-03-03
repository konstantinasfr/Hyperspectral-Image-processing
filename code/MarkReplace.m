%Convert space matrix storing label order to a row vector
function [assignment] = MarkReplace(MarkMat)

[rows,cols]=size(MarkMat);

assignment=zeros(1,cols);

for i=1:rows
    for j=1:cols
        if MarkMat(i,j)==1
            assignment(1,j)=i;
        end
    end
end

end
