function Coefs = least_squares_match(Data, Dictionary)

    Coefs = sparse(size(Dictionary,2),size(Data,2));
    
    for i = 1:size(Data, 2)
        
        bestWord = -1;
        minDistance = 0;
        
        dataL2norm = Data(:,i)'*Data(:,i);
        
        for j = 1:size(Dictionary, 2)
            
            dataWordDotProduct = Data(:,i)'*Dictionary(:,j);
            
            scale = dataWordDotProduct/dataL2norm;
            
            difference = scale*Data(:,i) - Dictionary(:,j);
            distance = norm(difference, 2);
            
            if (distance < minDistance || bestWord == -1)
                
                minDistance = distance;
                bestWord = j;
            end
        end
        
        Coefs(bestWord, i) = 1.0;
    end        