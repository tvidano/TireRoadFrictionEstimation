function [yhatk] = output_eqn(xbar,uprev)
% Calls output equation

% INPUTS: whatever is necessary to compute it --> USUALLY x[k|k-1], u[k-1]
% FOR THE NONLINEAR BATTERY MODEL --> Voc[k-1], I[k-1]

% OUTPUTS: y_hat[k]

global soc_intpts_OCV OCV_intpts R0

diff_OCV = xbar - soc_intpts_OCV;
[~,ind_OCV] = min(abs(diff_OCV));
Voc = OCV_intpts(ind_OCV);

yhatk = Voc - R0*uprev;

end

