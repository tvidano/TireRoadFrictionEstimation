function [xk] = state_eqn(model_param,xprev,uprev)
% Calls state equation 
%
% INPUTS: 
%   model_param {struct}: A struct containing model parameters,
%   xprev {vector}: Previous A Posteriori State Vector [k-1|k-1],
%   u {vector}: Previous Input Vector [k-1]
%
% OUTPUTS: 
%   x_hat {vector}: Current A Priori State Vector [k|k-1]

%-------------------------------------------------------
% Unpack Model Parameters:
del_t = model_param.del_t;
Cbat = model_param.Cbat;

xk = xprev - (del_t/Cbat)*uprev;

end

