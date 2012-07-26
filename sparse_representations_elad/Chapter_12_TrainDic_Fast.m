function [Dictionary,output] = Chapter_12_TrainDic_Fast(Data,param)
% ========================================================
% MOD and KSVD algorithms for dictionary learning
% This function uses the OMP package by Ron Rubinstein for speed-up 
% ========================================================
% INPUT ARGUMENTS:
% Data - an nXN matrix that contains N signals (Y), each of dimension n.
% param - a structure that contains all required parameters for the MOD execution.
%          Required fields are:
%               method - 'MOD' or 'KSVD'
%               K - the number of dictionary elements to train
%               numIteration - number of iterations to perform.
%               errorFlag - if =0, a fix number of coefficients is used for representation
%                               of each signal. If so, param.L must be specified as the number of
%                               representing atom.
%                               if =1, arbitrary number of atoms represent each signal, until a
%                               specific representation error is reached. If so, param.errorGoal
%                               must be specified as the allowed error.
%               L (optional, see errorFlag) - maximum coefficients to use in OMP
%               errorGoal (optional, see errorFlag) - allowed representation error
%               InitializationMethod - mehtod to initialize the dictionary, can
%                               be one of the following arguments:
%                                 * 'DataElements' (initialization by the signals themselves), or:
%                                 * 'GivenMatrix' (initialization by a given matrix param.initialDictionary).
%               initialDictionary (optional, see InitializationMethod) the matrix to initialize with
%               TrueDictionary (optional) - if specified, in each iteration the difference between
%                                  this dictionary and the trained one is measured and displayed.
%               displayProgress - if =1 progress information is displyed. If param.errorFlag==0,
%                                 the average repersentation error (RMSE) is displayed, while if
%                                 param.errorFlag==1, the average number of required coefficients for
%                                 representation of each signal is displayed.
% ========================================================
% OUTPUT ARGUMENTS:
%  Dictionary - The extracted dictionary of size nX(param.K).
%  output - A struct that contains information about the current run. It includes the fields:
%               CoefMatrix - The final coefficients matrix (it should hold that Data equals
%                                 approximately Dictionary*output.CoefMatrix.
%               ratio - If the true dictionary was defined (in synthetic experiments), this
%                                 parameter holds a vector of length param.numIteration that
%                                  includes the detection ratios in each iteration).
%               totalErr - The total representation error after each iteration (defined only
%                                 if param.displayProgress=1 and param.errorFlag = 0)
%               numCoef - A vector of length param.numIteration that includes the average
%                                 number of coefficients required for representation of each signal
%                                 (in each iteration) (defined only if param.displayProgress=1 and
%                                 param.errorFlag = 1)
% ========================================================

% initialization
counter=1;
if strcmp(param.InitializationMethod,'DataElements')
    Dictionary(:,1:param.K)=Data(:,1:param.K);
elseif strcmp(param.InitializationMethod,'GivenMatrix')
    Dictionary=param.initialDictionary;
end
Dictionary=Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));

% zeroing various result vectors
output.totalErr=zeros(1,2*param.numIteration+1);
output.numCoef=zeros(1,2*param.numIteration+1);
if (size(param.TrueDictionary)==size(Dictionary))
    displayErrorWithTrueDictionary=1;
    output.ratio=zeros(2*param.numIteration+1,1);
else
    displayErrorWithTrueDictionary=0;
end

% Sparse coding with the initial dictionary
if (param.errorFlag==0)
    % CoefMatrix=OMP(Dictionary,Data,param.L);
    CoefMatrix=omp(Dictionary'*Data,Dictionary'*Dictionary,param.L); 
else
    % CoefMatrix=OMPerr(Dictionary,Data, param.errorGoal);
    % ==> The next line uses Ron Rubinstein's very fast OMP package
    CoefMatrix=omp2(Dictionary'*Data,sum(Data.*Data),Dictionary'*Dictionary,n*param.errorGoal); 
end

% compute the errors
output.totalErr(counter)=sqrt(sum(sum((Data-Dictionary*CoefMatrix).^2))/...
                                        numel(Data));
output.numCoef(counter)=length(find(CoefMatrix))/size(Data,2);
disp(['Iteration ',num2str(0),':  Error=',num2str(output.totalErr(counter)), ...
    ', Average cardinality: ',num2str(output.numCoef(counter))]);
if (displayErrorWithTrueDictionary)
    output.ratio(counter)=DictionariesDistance(param.TrueDictionary,Dictionary);
    disp(strcat(['Iteration  ', num2str(0),' ratio of restored elements: ',...
                        num2str(output.ratio(counter))]));
end
counter=counter+1;

% Main Iterations
for iterNum=1:param.numIteration

    % Update the dictionary
    if strcmp(param.Method,'MOD')
        Dictionary=Data*CoefMatrix'*inv(CoefMatrix*CoefMatrix'+...
                         1e-7*speye(size(CoefMatrix,1)));
        sumDictElems=sum(abs(Dictionary));
        zerosIdx=find(sumDictElems<eps);
        Dictionary(:,zerosIdx) = randn(size(Dictionary,1),length(zerosIdx));
        Dictionary=Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));
        
    elseif strcmp(param.Method,'KSVD')
        for j=1:1:size(Dictionary,2)
            relevantDataIndices=find(CoefMatrix(j,:));
            if ~isempty(relevantDataIndices)
                tmpCoefMatrix=CoefMatrix(:,relevantDataIndices);
                tmpCoefMatrix(j,:)=0;
                errors=Data(:,relevantDataIndices)-Dictionary*tmpCoefMatrix;
                [betterDictionaryElement,singularValue,betaVector] = svds(errors,1);
                CoefMatrix(j,relevantDataIndices)=singularValue*betaVector';
                Dictionary(:,j)=betterDictionaryElement;
            end;
        end;
    end;
    
    % Compute the errors and display
    output.totalErr(counter)=sqrt(sum(sum((Data-Dictionary*CoefMatrix).^2))...
        /numel(Data));
    output.numCoef(counter)=length(find(CoefMatrix))/size(Data,2);
    disp(['Iteration ',num2str(iterNum),': Error=',...
        num2str(output.totalErr(counter)),', Average cardinality: ',...
        num2str(output.numCoef(counter))]);
    if (displayErrorWithTrueDictionary)
        output.ratio(counter)=DictionariesDistance(param.TrueDictionary,Dictionary);
        disp(strcat(['Iteration  ', num2str(iterNum),' ratio of restored elements: ',...
            num2str(output.ratio(counter))]));
    end
    counter=counter+1;
    
    % lean-up rarely used or too-close atoms
    T2=0.99; T1=3;
    Er=sum((Data-Dictionary*CoefMatrix).^2,1);
    G=Dictionary'*Dictionary;
    G=G-diag(diag(G));
    for j=1:1:size(Dictionary,2)
        if max(G(j,:))>T2
            alternativeAtom=find(G(j,:)==max(G(j,:)));
            [val,pos]=max(Er);
            Er(pos(1))=0;
            Dictionary(:,j)=Data(:,pos(1))/norm(Data(:,pos(1)));
            % CoefMatrix(alternativeAtom,:)=CoefMatrix(j,:)+CoefMatrix(alternativeAtom,:);
            % CoefMatrix(j,:)=0;
            % CoefMatrix(j,pos(1))=1;
            G=Dictionary'*Dictionary;
            G=G-diag(diag(G));            
        elseif length(find(abs(CoefMatrix(j,:))>1e-7))<=T1
            [val,pos]=max(Er);
            Er(pos(1))=0;
            Dictionary(:,j)=Data(:,pos(1))/norm(Data(:,pos(1)));
            % CoefMatrix(j,:)=0;
            % CoefMatrix(j,pos(1))=1;
            G=Dictionary'*Dictionary;
            G=G-diag(diag(G));
        end;
    end;
    
    % Sparse coding: find the coefficients
    if (param.errorFlag==0)
        % CoefMatrix=OMP(Dictionary,Data,param.L);
        CoefMatrix=omp(Dictionary'*Data,Dictionary'*Dictionary,param.L); 
    else
        % CoefMatrix=OMPerr(Dictionary,Data,param.errorGoal);
        CoefMatrix=omp2(Dictionary'*Data,sum(Data.*Data),Dictionary'*Dictionary,n*param.errorGoal); 
    end

    % Compute the errors and display
    output.totalErr(counter)=sqrt(sum(sum((Data-Dictionary*CoefMatrix).^2))...
        /numel(Data));
    output.numCoef(counter)=length(find(CoefMatrix))/size(Data,2);
    disp(['Iteration ',num2str(iterNum),': Error=',...
        num2str(output.totalErr(counter)),', Average cardinality: ',...
        num2str(output.numCoef(counter))]);
    if (displayErrorWithTrueDictionary)
        output.ratio(counter)=DictionariesDistance(param.TrueDictionary,Dictionary);
        disp(strcat(['Iteration  ', num2str(iterNum),' ratio of restored elements: ',...
            num2str(output.ratio(counter))]));
    end
    counter=counter+1;
    
end
output.CoefMatrix=CoefMatrix;

return;

% ========================================================
% ========================================================

function [ratio]=DictionariesDistance(original,new)

T1=0.01;
catchCounter=0;
distances=abs(original'*new);
for i=1:1:size(original,2)
    minValue=1-max(distances(i,:));
    catchCounter=catchCounter+(minValue<T1);
end;
ratio=100*catchCounter/size(original,2);
return; 

% ========================================================
% ========================================================

function [A]=OMPerr(D,X,errorGoal)
% ========================================================
% Sparse coding of a group of signals based on a given dictionary and specified representation
% error to get.
% input arguments: D - the dictionary
%                           X - the signals to represent
%                           errorGoal - the maximal allowed representation error
% output arguments: A - sparse coefficient matrix.
% ========================================================
[n,P]=size(X);
[n,K]=size(D);
E2 = errorGoal^2*n;
maxNumCoef = n/2;
A = sparse(size(D,2),size(X,2));
for k=1:1:P,
    a=[];
    x=X(:,k);
    residual=x;
    indx = [];
    a = [];
    currResNorm2 = sum(residual.^2);
    j = 0;
    while currResNorm2>E2 && j < maxNumCoef,
        j = j+1;
        proj=D'*residual;
        pos=find(abs(proj)==max(abs(proj)));
        pos=pos(1);
        indx(j)=pos;
        a=pinv(D(:,indx(1:j)))*x;
        residual=x-D(:,indx(1:j))*a;
        currResNorm2 = sum(residual.^2);
    end;
    if (~isempty(indx))
        A(indx,k)=a;
    end
end;
return;

% ========================================================
% ========================================================

function [A]=OMP(D,X,L)
% ========================================================
% Sparse coding of a group of signals based on a given dictionary and specified number of
% atoms to use.
% input arguments: D - the dictionary
%                           X - the signals to represent
%                           errorGoal - the maximal allowed representation error
% output arguments: A - sparse coefficient matrix.
% ========================================================
[n,P]=size(X);
[n,K]=size(D);
for k=1:1:P,
    a=[];
    x=X(:,k);
    residual=x;
    indx=zeros(L,1);
    for j=1:1:L,
        proj=D'*residual;
        [maxVal,pos]=max(abs(proj));
        pos=pos(1);
        indx(j)=pos;
        a=pinv(D(:,indx(1:j)))*x;
        residual=x-D(:,indx(1:j))*a;
        if sum(residual.^2) < 1e-6
            break;
        end
    end;
    temp=zeros(K,1);
    temp(indx(1:j))=a;
    A(:,k)=sparse(temp);
end;
return;

