function constraints = writeConstraintsFile2(Aeq,b,senseArray,pathToSave)

% create constraints.txt
% # contains one linear constraints per row on the allowed labels in the following form:
% #
% #   <coef>*<var_num> [<coef>*<var_num> ... ] <rel> <value>
% #
% # where
% #   <coef>    ... a real number
% #   <var_num> ... the number of the variable, in accordance to label.txt and features.txt
% #   <rel>     ... the relation, one of "<=", "==", ">="
% #   <value>   ... a real number
% 
% 1*0 1*1 == 1 # y_0 + y_1 == 1

filename = 'constraints.txt';
% filename = fullfile(pathToSave,filename);
fileID = fopen(filename,'w');
[numConstraints,numVar] = size(Aeq);
[r,c,s] = find(Aeq);
% var_num starts from zero

for i=1:numConstraints
    rel1 = senseArray(i);
    switch rel1
        case '='
            rel = ' == ';
        case '>'
            rel = ' >= ';
        case '<'
            rel = ' <= ';
    end
    
    % [ r,c,coeff] = find(Aeq(i,:));
    constr_i_logicalInd = (r==i);
    var_num = c(constr_i_logicalInd) -1;
    % var_num = c - 1; % var_num starts with zero
    coeff = s(constr_i_logicalInd);
    
    term = '%4.4f*%d ';
    numTerms = numel(coeff);
    bstr = int2str(b(i));
%     formatSpec = repmat(term,1,numTerms);
%     
%     
%     formatSpec = strcat(formatSpec,rel);
%    
%     formatSpec = strcat(formatSpec,bstr);
%     
%     sprintf(formatSpec,coeff',var_num')
    str = [];
    for j=1:numTerms
       str1 = sprintf(term,coeff(j),var_num(j));
       %str1 = horzcat(' ',str1);
       str = horzcat(str,str1);
    end
    
    str = horzcat(str,rel);
    str = horzcat(str,bstr);
    str = horzcat(str,'\n');
    fprintf(fileID,str);
       
end
fclose(fileID);
constraints = 0;