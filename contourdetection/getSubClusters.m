function subClusterIDs = getSubClusters(clusPixInds,sizeR,sizeC)

% Input:
%   pixInds of a clustered node

% Output:
%   subClusterIDs: nx2 matrix with first col the pixInds and the second
%   col the corresponding subcluster ids

numPixInds = numel(clusPixInds);

subClusterIDs = zeros(numPixInds,2);
subClusterIDs(:,1) = clusPixInds;
neighborhood = zeros(numPixInds,4);

if(numPixInds==0)
    subClusterIDs = 0;
else
    % get 4nh for each clus node pixel
    for i=1:numPixInds
        neighbors = getNeighbors(clusPixInds(i),sizeR,sizeC);
        neighbors = intersect(clusPixInds,neighbors);
        if(numel(neighbors)>0)
            neighborhood(i,1:numel(neighbors)) = neighbors;
        end
    end    
    
    % assign cluster ids
    k = 0;
    for i=1:numPixInds
        if(subClusterIDs(i,2)==0)
            % assign id
            k = k+1;
            subClusterIDs(i,2)= k;
            K=0;
        else
            % get its id and assign it to its neighbors
            K = subClusterIDs(i,2);            
        end
        
        % get neighbors and assign them this id
        % if a neighbor has a different id already, merge the two ids
        neighbors = neighborhood(i,:);
        neighbors = neighbors(neighbors>0);
        if(numel(neighbors)>0)
            [~,pixLIDs] = intersect(clusPixInds,neighbors);
            if(K==0)
                % assign k
                subClusterIDs(pixLIDs,2) = k;
            else
                % assign K
                subClusterIDs(pixLIDs,2) = K;
            end
        end
        
    end
end