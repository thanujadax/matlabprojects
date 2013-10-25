function c_wsIDsInCell = getCells2WSregions(labelImg_indexed,ws,numLabels)

c_wsIDsInCell = cell(numLabels,1);

for i=1:numLabels
    pixelsForLabel_i = (labelImg_indexed==i);
    wsIDsForLabel_i = ws(pixelsForLabel_i);
    wsIDsForLabel_i = unique(wsIDsForLabel_i);
    wsIDsForLabel_i = wsIDsForLabel_i(wsIDsForLabel_i>0);
    c_wsIDsInCell{i} = wsIDsForLabel_i; 
end
