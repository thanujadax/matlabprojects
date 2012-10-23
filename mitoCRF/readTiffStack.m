% Read Tiff Stack and return 3D matrix
% Faster version
% uses tifflib.mexa64.MEX file which is to be copied in the path
% REF: http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/

function FinalImage = readTiffStack(FileTif)

%FileTif='../data/mitoData/stem.tiff';
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint8');
FileID = tifflib('open',FileTif,'r');
rps = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);
 
for i=1:NumberImages
   tifflib('setDirectory',FileID,i);
   % Go through each strip of data.
   rps = min(rps,mImage);
   for r = 1:rps:mImage
      row_inds = r:min(mImage,r+rps-1);
      stripNum = tifflib('computeStrip',FileID,r);
      tmp = tifflib('readEncodedStrip',FileID,stripNum);
      % tmp receives 3 copies of the same image as a 3D array
      FinalImage(row_inds,:,i) = tmp(:,:,1); 
   end
end

tifflib('close',FileID);

