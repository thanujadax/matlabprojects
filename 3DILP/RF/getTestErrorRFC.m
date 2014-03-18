function err = getTestErrorRFC(model,X_tst,Y_tst)
Y_hat = classRF_predict(X_tst,model);
err = length(find(Y_hat~=Y_tst))/length(Y_tst);
% fprintf('\nexample 1: error rate %f\n', err);