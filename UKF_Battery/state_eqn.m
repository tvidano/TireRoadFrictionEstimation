function [xk] = state_eqn(xprev,uprev)
% Calls state equation 

% INPUTS: x[k-1|k-1], u[k-1]

% OUTPUTS: x_hat[k|k-1]

%-------------------------------------------------------
global del_t Cbat

xk = xprev - (del_t/Cbat)*uprev;

end

