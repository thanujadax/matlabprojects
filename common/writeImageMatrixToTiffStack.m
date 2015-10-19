function writeImageMatrixToTiffStack(imgdata,outputFileName)
% remove outputfile if already exists
if exist(outputFileName,'file') == 2
    delete(outputFileName);
end

    t = Tiff(outputFileName,'w');
    tagstruct.ImageLength = size(imgdata,1);
    tagstruct.ImageWidth = size(imgdata,2);
    tagstruct.Photometric = Tiff.Photometric.RGB;
    tagstruct.BitsPerSample = 8;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    t.setTag(tagstruct);
    t.write(imgdata);
    
% for K=1:length(imgdata(1, 1, :))
%     imwrite(imgdata(:, :, K), outputFileName, 'WriteMode', 'append',  'Compression','none');
% 
% end

    t.close();