function f = getILPObj3Dtoy(edgeUnary,nodeAngleCosts,...
            regionUnary,T_unary,...
            w_off_e,w_on_e,w_off_n,w_on_n,w_off_r,w_on_r,w_on_T,w_off_T)
% Inputs:
% T_unary: unary probability for each T to be active. Learned from Rf

%% 2D costs
f = getILPObjectiveVectorParametric(edgeUnary,nodeAngleCosts,...
            regionUnary,w_off_e,w_on_e,w_off_n,w_on_n,w_off_r,w_on_r);

%% 3D costs
% T variables



