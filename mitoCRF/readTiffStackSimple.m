% Reads tiff image stack and returns 3D matrix
function FinalImage = readTiffStackSimple(FileTif)

%FileTif='ImageStack.tif';
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint8');
 
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   %TifLink.setTag('SubFileType', Tiff.SubFileType.Page);
   % tmp receives 3 copies of the same image as a 3D array
   tmp=TifLink.read();
   FinalImage(:,:,i)=tmp(:,:,1);
end
TifLink.close();