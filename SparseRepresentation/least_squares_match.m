function Coefs = least_squares_match(Data, Dictionary)

    Coefs = sparse(size(Dictionary,2),size(Data,2));
    
    for i = 1:size(Data, 2)
        
        bestWord = -1;
        minDistance = 0;
        bestScale = 0;
        
        dataL2norm = Data(:,i)'*Data(:,i);
        
        % find the best dictionary word regarding the L2 norm without
        % consideration of the scale
        for j = 1:size(Dictionary, 2)
            
            dataWordDotProduct = Data(:,i)'*Dictionary(:,j);
            
            % scaling factor for the data-point to minimize the L2 norm to
            % the current dictionary word
            scale = dataWordDotProduct/dataL2norm;
            
            difference = scale*Data(:,i) - Dictionary(:,j);
            distance = norm(difference, 2);
            
            if (distance < minDistance || bestWord == -1)
                
                minDistance = distance;
                bestWord = j;
                bestScale = scale;
            end
        end
        
        % the coefficient is the inverse scale of the data-point
        Coefs(bestWord, i) = 1/bestScale;
    end        