%function model = classRF_train(X,Y,ntree,mtry, extra_options)
    %requires 2 arguments and the rest 2 are optional
    %X: data matrix
    %Y: target values
    %ntree (optional): number of trees (default is 500)
    %mtry (default is max(floor(D/3),1) D=number of features in X)
    %there are about 14 odd options for extra_options. 
    % Refer to tutorial_ClassRF.m to examine them

%function Y_hat = classRF_predict(X,model)
    %requires 2 arguments
    %X: data matrix
    %model: generated via classRF_train function
    
ntree = 500;

% load /home/thanuja/Dropbox/RESULTS/sparserepresentations/dict/SparseCoeffMat_reconst_01.mat
% load /home/thanuja/Dropbox/RESULTS/sparserepresentations/dict/labels_mem_01_256.mat

data = full(sparsecoeff);
data = data';
X = data;
Y = labelVec;

% train RF    
t1 = cputime
model = classRF_train(X,Y,ntree);
t2 = cputime

timetaken = t2-t1


% predict new image labels
% mem10


