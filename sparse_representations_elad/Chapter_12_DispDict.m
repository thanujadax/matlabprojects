function I=Chapter_12_DispDict(D,numRows,numCols,X,Y,sortVarFlag)

% This function displays the dictionary D which consists of image patches.
% The dicitonary has numRows*numCols atoms to show, and those are presented
% as 2D array of patches of size X-by-Y pixels. The sortVarFlag flag
% enables a sorting of the atoims by their variance

borderSize=1;

% Preparing the image 
sizeForEachImage =sqrt(size(D,1))+borderSize;
I=zeros(sizeForEachImage*numRows+borderSize,...
            sizeForEachImage*numCols+borderSize,3);
I(:,:,1)=0;
I(:,:,2)=0; 
I(:,:,3)=1; 

% Stretching the atoms
for counter=1:size(D,2)
    D(:,counter)=D(:,counter)-min(D(:,counter));
    if (max(D(:,counter)))
        D(:,counter)=D(:,counter)./max(D(:,counter));
    end
end

if (sortVarFlag)
    vars=var(D);
    [V,indices]=sort(vars');
    indices=fliplr(indices);
    D=[D(:,1:sortVarFlag-1),D(:,indices+sortVarFlag-1)];
    signs=sign(D(1,:));
    signs(find(signs==0))=1;
    D=D.*repmat(signs,size(D,1),1);
    D=D(:,1:numRows*numCols);
end

% Fill the image with the atoms
counter=1;
for j = 1:numRows
    for i = 1:numCols
        I(borderSize+(j-1)*sizeForEachImage+1:j*sizeForEachImage,...
           borderSize+(i-1)*sizeForEachImage+1:i*sizeForEachImage,1)...
           =reshape(D(:,counter),X,Y);
        I(borderSize+(j-1)*sizeForEachImage+1:j*sizeForEachImage,...
           borderSize+(i-1)*sizeForEachImage+1:i*sizeForEachImage,2)...
           =reshape(D(:,counter),X,Y);
        I(borderSize+(j-1)*sizeForEachImage+1:j*sizeForEachImage,...
           borderSize+(i-1)*sizeForEachImage+1:i*sizeForEachImage,3)...
           =reshape(D(:,counter),X,Y);
        counter = counter+1;
    end
end

return;
