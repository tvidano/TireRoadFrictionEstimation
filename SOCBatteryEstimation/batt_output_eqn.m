function [yhatk] = batt_output_eqn(model_param,xbar,uprev)
% Calls OUTPUT EQUATION
%
% INPUTS: 
%   model_param {struct}: A struct containing model parameters, 
%   x {vector}: A Posteriori State Vector [k|k-1],
%   u {vector}: A Posteriori Input Vector [k-1],
%
% OUTPUTS: 
%   y_hat {vector}: Output vector [k]

% Unpack Model Parameters:
soc_intpts_OCV = model_param.soc_intpts_OCV;
OCV_intpts = model_param.OCV_intpts;
R0 = model_param.R0;

% Equations:
diff_OCV = xbar - soc_intpts_OCV;
[~,ind_OCV] = min(abs(diff_OCV));
Voc = OCV_intpts(ind_OCV);

yhatk = Voc - R0*uprev;

end

