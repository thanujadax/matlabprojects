function constraints = writeConstraintsFile(Aeq,senseArray)

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

[numConstraints,numVar] = size(Aeq);

for i=1:numConstraints
    rel1 = senseArray(i);
    switch rel1
        case '='
            rel = '==';
        case '>'
            rel = '>=';
        case '<'
            rel = '<=';
    end
    
    
    
end
